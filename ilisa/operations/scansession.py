import sys
import os
import time
import datetime
import shutil
import yaml
import argparse
import logging

import ilisa
import ilisa.operations
from ilisa.operations import USER_CACHE_DIR
import ilisa.operations.directions as directions
import ilisa.operations.modeparms as modeparms
import ilisa.operations.programs as programs
from ilisa.operations.stationdriver import StationDriver, waituntil
from ilisa.operations.scan import still_time_fun, LScan

OBSLOGFILE = "ilisa_cmds.log"  # Log of each ilisa_obs or ilisa_adm CLI run
_LOGGER = logging.getLogger(__name__)


def projid2meta(projectid):
    """\
    Get the project metadata for project with projectid

    Parameters
    ----------
    projectid : str
        The ID for project.

    Returns
    -------
    projectmeta : dict
        The projects metadata.
    accessfiles : dict
        The access files keyed on station ID.
    projectfile : str
        Path to the file for which the above information was take.
    proj_srcs : list
        List of source names and directions.
    """
    # Setup projectmeta:
    if projectid is not None:
        projectfile = "project_" + projectid + ".yml"
        projectfile = os.path.join(ilisa.operations.USER_CONF_DIR, projectfile)
        with open(projectfile) as projectfilep:
            projectprofile = yaml.safe_load(projectfilep)
        projectmeta = projectprofile['project']
        accessfiles = projectprofile['accessfiles']
        proj_srcs = projectprofile.get('sources')
    else:
        projectmeta = {'observer': None, 'name': None}
        accessfiles = None
        projectfile = None
        proj_srcs = []
    return projectmeta, accessfiles, projectfile, proj_srcs


def process_scansess(sesscans_in):
    """\
    Function for parsing a station session schedule

    Parameters
    ----------
    sesscans_in : dict
        Input session scans.

    Returns
    -------
    sessmeta : dict
        Processed metadata for sesscan_in.
    generate_scans() : generator
        scans settings.
    """
    # Set the session_id to something
    session_id = sesscans_in.get('session_id')
    mockrun = sesscans_in.get('mockrun', False)
    projectid = sesscans_in.get('projectid', '0')
    note = sesscans_in.get('note')
    stnid = sesscans_in.get('station')

    # Initialize processed station session schedule
    sessmeta = {'session_id': session_id,
                'projectid': projectid,
                'station': stnid,
                'note': note,
                'mockrun': mockrun,
                'cli_start': sesscans_in['cli_start'],
                'start': sesscans_in.get('start')
                # ,'scans': []
                }

    # Process the start time for the session
    # # cli_start overrides scan session start (if not None)
    if sessmeta['cli_start']:
        # Command Line Interface requested starttime takes priority
        sessmeta['start'] = modeparms.as_asapdatetime(sessmeta['cli_start'])
    if not sessmeta['start']:
        sessmeta['start'] = 'ASAP'

    def generate_scans():
        # Margin of time between two scans in seconds (fastest safe switching)
        margintime = datetime.timedelta(seconds=13.0)

        # Initialize "previous" values for scan loop
        scanstarttimeprev = modeparms.timestr2datetime(sessmeta['start'])
        duration_totprev = 0
        for scan in sesscans_in['scans']:
            # Get a title for this scan
            scan_id = scan.get('id')

            # Prepare observation arguments:

            # - Compute `starttime` of this scan
            #   `after` converted to `starttime` using formula:
            #      `starttime` = `starttimeprev` + `after`
            scanstarttime = scan.get('starttime')
            if scan == sesscans_in['scans'][0] and \
                    (not scanstarttime or scanstarttime == 'ASAP'):
                # First scan set to session start if not set already or ASAP
                scanstarttime = sessmeta['start']
            scanstarttime_guess = scanstarttime
            if not scanstarttime_guess or scanstarttime_guess == 'ASAP':
                dur_dprev = datetime.timedelta(seconds=duration_totprev)
                after = scan.get('after')
                if after:
                    # 'after' field overrides scanstarttime set to 'ASAP'
                    after_delta = modeparms.hmsstr2deltatime(after)
                    if after_delta < (dur_dprev + margintime):
                        time2nxtscan = (dur_dprev, after_delta, margintime)
                        _LOGGER.warning(
                            'No time for scan: after {} - dur {} < marg {}'
                            .format(*time2nxtscan))
                    scanstarttime = scanstarttimeprev + after_delta
                    scanstarttime_guess = scanstarttime
                else:
                    # Now we will have to predict scanstarttime
                    # and it should include margintime between scans
                    scanstarttime = 'ASAP'
                    scanstarttime_guess = (scanstarttimeprev + dur_dprev
                                           + margintime)

            # - Duration total
            duration_in = scan.get('duration')
            duration_tot = modeparms.hmsstr2deltatime(str(duration_in)
                                                      ).total_seconds()
            # duration_tot = eval(str(scan['duration']))
            file_dur = scan.get('file_dur')

            # - Beam
            beam = scan.get('beam', {})
            # -- Freq
            freqspec = beam.get('freqspec')
            # -- Pointing
            pointing = beam.get('pointing')
            # -- direction: alternative to pointing but can't be name
            direction = beam.get('direction')
            # -- Source name
            source = beam.get('source')

            # Postprocess beam to get direction
            if not direction:
                if pointing:
                    direction = directions.normalizebeamctldir(pointing)
                elif source:
                    direction = directions.std_pointings(source)
            if not directions.check_directionstr(direction):
                raise ValueError('direction {} is not correct format'
                                 .format(direction))
            if not source:
                if pointing:
                    source = pointing

            # - Record
            #     defaults
            rec = scan.get('rec', [])
            acc = False
            if 'acc' in rec:
                acc = True
            bfs = False
            if 'bfs' in rec:
                bfs = True
            bsx_stat = None
            if 'bst' in rec:
                bsx_stat = 'bst'
            elif 'sst' in rec:
                bsx_stat = 'sst'
            elif 'xst' in rec:
                bsx_stat = 'xst'
            # - Integration for rec bsx
            integration = scan.get('integration', 1.0)
            if integration > duration_tot:
                raise ValueError("Integration longer than duration ({}>{})."
                                 .format(integration, duration_tot))

            # If dedicated observation program chosen, set it up
            # otherwise run main obs program
            obsprog = scan.get('obsprog')

            # Collect observation parameters specified
            obsargs_in = {'beam':
                              {'freqspec': freqspec,
                               'pointing': pointing,
                               'direction': direction,
                               'source': source
                               },
                          'rec': rec,
                          'acc': acc,
                          'bfs': bfs,
                          'bsx_stat': bsx_stat,
                          'integration': integration,
                          'duration': duration_tot,
                          'file_dur': file_dur,
                          'starttime': modeparms.astimestr(scanstarttime),
                          'starttime_guess': \
                              modeparms.astimestr(scanstarttime_guess),
                          'id': scan_id
                          }
            obsargs_in.update({'obsprog': obsprog})
            yield obsargs_in
            # Next scan use current time as previous time and current time
            scanstarttimeprev = scanstarttime_guess
            duration_totprev = duration_tot

    return sessmeta, generate_scans()


def lcu_services(scan):
    """
    Find out what LCU services are needed

    Parameters
    ----------
    scan : dict
        The scan specification.

    Return
    ------
    service : str
        The service required.
    """
    beam = scan.get('beam')
    if beam.get('pointing'):
        return 'beam'
    elif beam.get('freqspec') and \
            modeparms.FreqSetup(beam['freqspec']).rcumodes[0] > 4:
        return 'tof'
    return None


def check_sess_start_passed(sessmeta):
    """\
    Check if Session can be started in future or if starttime is now in the past

    Parameters
    ----------
    sessmeta : dict
        The ScanSession metadata.
    """
    paststart = False
    if sessmeta['start'] != 'ASAP':
        utcnow = datetime.datetime.utcnow()
        if sessmeta['start'] < utcnow:
            paststart = True
    return paststart


class ScanSession(object):
    """Class that runs a session on a station."""

    def __init__(self, stndrv, session_id=None):
        """Initialize Session."""
        self.stndrv = stndrv
        if session_id is None:
            session_id = self.make_session_id()
        self.session_id = session_id
        latestdatafile = 'latestdatafiles_' + self.stndrv.get_stnid() + '.txt'
        self.latestdatafile = os.path.join(USER_CACHE_DIR, latestdatafile)
        # Create the file with 1st line noting ongoing session
        # This should be removed when session is finished.
        with open(self.latestdatafile, 'w') as f:
            f.writelines(['ONGOING\n'])

    def set_stn_session_id(self, parent_session_id):
        self.stn_sess_id = '{}_{}'.format(parent_session_id, self.stndrv.get_stnid())

    def get_stn_session_id(self):
        return self.stn_sess_id

    def make_session_id(self, ref_dattim=datetime.datetime.utcnow()):
        """Make a session ID based on time of creation.
        session_id has format 'sid<CT>' where <CT> is the datetime
        in the format '%Y%m%dT%H%M%S' of the time of creation.
        """
        session_id = "sid{}".format(
            modeparms.normalizetimestr(ref_dattim, dt_fformat='%Y%m%dT%H%M%S'))
        return session_id

    def get_datastorepath(self):
        return self.stndrv.dru_data_root

    def _projsubpath(self):
        projsubpath = os.path.join('Projects',
                                   'proj{}'.format(self.projectmeta['id']))
        return projsubpath

    def _sesssubpath(self):
        sesssubpath = os.path.join(self._projsubpath(),
                                   'sess_{}'.format(self.get_stn_session_id()))
        return sesssubpath

    def get_sesspath(self):
        sesspath = os.path.join(self.get_datastorepath(), self._sesssubpath())
        return sesspath

    def save_scansess(self, sched):
        sesspath = self.get_sesspath()
        ses_sched_file = os.path.join(sesspath, 'ScanSess.yml')
        with open(ses_sched_file, 'w') as f:
            yaml.dump(sched, f, explicit_start=True)

    def run_scansess(self, sesscans_in, session_id=None):
        """\
        Run local session given a stn_ses_schedule dict.

        Dispatches to the stationdriver to execute session operations
        """

        def _update_latestdatafile(rec_name, rec_state):
            """\
            Update the latestdata file for the given recording if needed
            """
            if rec_state and not scanrecpath[rec_name]:
                scanrecpath[rec_name] = lscan.scanresult[rec_name].scanrecpath
                if scanrecpath[rec_name]:
                    with open(self.latestdatafile, 'a') as f:
                        f.write(scanrecpath[rec_name])
                        f.write('\n')

        sessmeta, scans_iter = process_scansess(sesscans_in)
        if check_sess_start_passed(sessmeta):
            raise ValueError('Starttime in the past')
        # Starttime handling
        self.session_id = session_id
        if not self.session_id:
            self.session_id = self.make_session_id(sessmeta['start'])
        self.projectmeta, _, _, _ = projid2meta(sessmeta['projectid'])
        self.set_stn_session_id(self.session_id)
        # Set where ldata should be put after recording on LCU
        sesspath = self.get_sesspath()
        bfdsesdumpdir = self._sesssubpath()
        # Sess start delta time handling (=boot time plus beam)
        sesstart_dt = datetime.timedelta(minutes=1)
        # Wait until it is time to bootup
        waituntil(sessmeta['start'], sesstart_dt)
        self.stndrv.goto_observingstate()
        _LOGGER.info("Started {}".format(self.session_id))

        scans_done = []
        for scan in scans_iter:
            if scan['obsprog'] is not None:
                _mockrun = False
                if not _mockrun:
                    scanresult = programs.record_obsprog(self.stndrv, scan)
                else:
                    scanresult = {}
            else:
                # Only pointing used not source name but it's in scan metadata
                pointing_spec = {'pointing': scan['beam']['pointing'],
                                 'direction': scan['beam']['direction'],
                                 'source': scan['beam']['source']}
                acc = scan['acc']
                bfs = scan['bfs']
                bsx_stat = scan['bsx_stat']
                integration = scan['integration']
                freqspec = scan['beam']['freqspec']
                freqsetup = modeparms.FreqSetup(freqspec)
                scan_dur = scan['duration']

                # Calculate scan schedule fundamental timings
                starttime = modeparms.timestr2datetime(scan['starttime'])
                _lcu_services = lcu_services(scan)

                margin_scan_start = datetime.timedelta(seconds=20)
                startedtime = waituntil(starttime, margin_scan_start)
                _lag_delta = startedtime - starttime  # positive for proper lag
                if scan['starttime'] == 'ASAP':
                    _lag_delta = datetime.timedelta(0.0)
                if _lag_delta > datetime.timedelta(0.0):
                    _LOGGER.warning("Scan start is lagging by {}"
                                    .format(_lag_delta))
                    _dur_corr = _lag_delta.total_seconds()
                    # Prioritize start of next scan over this scan's duration
                    # by correcting `scan_dur` by lag:
                    _LOGGER.info("Correcting duration by {}s for lag"
                                 .format(_dur_corr))
                    scan_dur -= _dur_corr
                    if scan_dur <= 0.0:
                        scan['id'] = 'NoID' if not scan.get('id', None) else scan['id']
                        _LOGGER.warning('No time for scan,'
                                        ' so skipping this one.')
                        scan['id'] += '_SKIPPED'
                        scans_done.append(scan)
                        continue

                # Initialize LScan
                lscan = LScan(self.stndrv, bsx_stat, freqsetup, pointing_spec,
                              scan_dur, integration, starttime,
                              acc=acc, bfs=bfs, destpath=sesspath,
                              destpath_bfs=bfdsesdumpdir,
                              file_dur=scan['file_dur'], scan_id=scan['id'])
                _LOGGER.info("Started LScan:: {}"
                             .format(lscan.describe_scan()))
                subscan = iter(lscan)
                # Start the subscan
                scanrecpath = {'acc': None, 'bfs': None, 'bsx': None}
                next(subscan)

                # Compute stop time as now + scan_dur
                stoptime = datetime.datetime.utcnow() + datetime.timedelta(
                    seconds=scan_dur)
                _LOGGER.info('Will stop @ {}'.format(stoptime))
                stop_cond = still_time_fun(stoptime)
                stop_scan_cond = stop_cond()
                while stop_scan_cond:
                    try:
                        _LOGGER.debug('scansession: IN TOP SUBSCAN LOOP')
                        _ = subscan.send(stop_scan_cond)
                    except StopIteration:
                        stop_scan_cond = False
                    else:
                        stop_scan_cond = stop_cond()
                    list(map(_update_latestdatafile,
                             *zip(('acc', acc), ('bfs', bfs), ('bsx', bsx_stat)))
                         )
                _LOGGER.debug('scansession: END SUBSCAN LOOP')
                lscan.close()
                # If no recording then simply wait for duration_tot
                if not bfs and not modeparms._xtract_bsx(bsx_stat) and not acc:
                    _sleepfor = (datetime.timedelta(seconds=scan_dur)
                                 + margin_scan_start)
                    _LOGGER.info('Not recording for {}s'.format(_sleepfor))
                    time.sleep(_sleepfor)
                scanresult = lscan.scanresult
            scan['id'] = scanresult.get('scan_id', None)
            scanpath_scdat = scanresult.get('scanpath_scdat', None)
            _LOGGER.info("Saved scan here: {}".format(scanpath_scdat))
            scan_ended_at = datetime.datetime.utcnow()
            duration_actual = scan_ended_at - startedtime
            duration_req = datetime.timedelta(seconds=scan['duration'])
            _LOGGER.info("End LScan:: id: {}".format(scan['id']))
            _LOGGER.debug("Request dur={}, ".format(duration_req)
                          + "Actual dur={}".format(duration_actual))
            scans_done.append(scan)
        # Collate session metadata with scan meta data as scansession metadata
        sessmeta['scans'] = scans_done
        # and save it:
        self.save_scansess(sessmeta)
        # Update LATESTDATAFILE that session is no longer on-going:
        with open(self.latestdatafile, 'r') as f:
            filecontents = f.readlines()
            # 1st line containing ongoing status. Remove it.
            _ongoing = filecontents.pop(0)
        with open(self.latestdatafile, 'w') as f:
            f.writelines(filecontents)


def get_proj_stn_access_conf(projid, stnid=None):
    """\
    Get access conf for project projid on stnid

    Parameters
    ----------
    projid : str
        Project ID.
    stnid : str
        Station ID. If set to None, we try to get default station for projid.

    Returns
    -------
    ac : dict
        Access info for given projid and stnid.

    Raises
    ------
    RuntimeError
        In projid's project file: if there is no accessfile defined for stnid
        or if stnid is not given and there are no stnid's.
    """
    projectmeta, accessfiles, projectfile, _ = projid2meta(projid)
    if stnid is None:
        # Try to get station from accessfiles in project config file:
        #    1st key in accessfiles dict taken as default stnid
        try:
            stnid = list(accessfiles.keys()).pop(0)
        except:
            raise RuntimeError("No stations for project {}".format(projid))
    try:
        # See if station has an access config file
        acf_name = accessfiles[stnid]
    except:
        raise RuntimeError("Station {} not found for project {}"
                           " (Check project file {})"
                           .format(stnid, projid, projectfile))
    userilisadir = ilisa.operations.USER_CONF_DIR
    acf_path = os.path.join(userilisadir, acf_name)
    with open(acf_path) as acffp:
        ac = yaml.safe_load(acffp)
    return ac


def check_scan_sess(scansess_in):
    """\
    Check scan session syntax

    Parameters
    ----------
    sessmeta : dict
        The ScanSession metadata.
    """
    sessmeta, scans_obsargs = process_scansess(scansess_in)
    print(yaml.dump(sessmeta, default_flow_style=False), end='')
    sessscans = {'scans': []}
    try:
        for scan_obsargs in scans_obsargs:
            sessscans['scans'].append(scan_obsargs)
    except ValueError as ve:
        print("ValueError:", ve)
    else:
        print(yaml.dump(sessscans, default_flow_style=False), end='')
    lastscan = sessscans['scans'][-1]
    endtime = modeparms.timestr2datetime(lastscan['starttime_guess']) \
              + datetime.timedelta(seconds=lastscan['duration'])
    print('end:', endtime)
    starttime = modeparms.timestr2datetime(sessscans['scans'][0]['starttime_guess'])
    print('duration_total:', endtime - starttime)
    if check_sess_start_passed(sessmeta):
        _LOGGER.warning("Session starttime {} is in the past."
                        .format(sessmeta['start']))


def obs(scansess_in, sac):
    """\
    Observe a Scansession

    Parameters
    ----------
    scansess_in : dict
        The ScanSession specification.
    """
    cli_start = scansess_in['cli_start']
    mockrun = scansess_in['mockrun']
    projectid = scansess_in['projectid']
    file = scansess_in['file']
    issued_at = datetime.datetime.utcnow().isoformat(timespec='seconds')
    # Initialize stationdriver
    try:
        stndrv = StationDriver(sac['LCU'], sac['DRU'], mockrun=mockrun)
    except (ConnectionError, AssertionError) as err:
        _LOGGER.error("Cannot setup {}'s LCU correctly"
                      .format(sac['LCU']['stnid']))
        _LOGGER.error(err)
        sys.exit(1)
    scansess_in['station'] = stndrv.get_stnid()
    scnsess = ScanSession(stndrv)
    stndrv._lcu_interface._fake_slow_conn = 0  # Fake a slower connection
    if stndrv._lcu_interface._fake_slow_conn != 0:
        _LOGGER.warning('Faking slow connection. Delay {}s.'.format(
            stndrv._lcu_interface._fake_slow_conn))
    scnsess.run_scansess(scansess_in)
    cmd = 'obs:' + file
    with open(OBSLOGFILE, 'a') as lgf:
        if mockrun:
            priority_fld = 'M'
        else:
            priority_fld = '0'
        lgf.write("{} {} {} {}".format(issued_at, cli_start, priority_fld,
                                       projectid)
                  + " {} {} {} '{}'\n".format(stndrv.get_stnid(), cmd,
                                              scnsess.session_id,
                                              scansess_in['note']))
    return scnsess


def main_cli():
    """\
    CLI to send hi-level commands to a LOFAR station.
    """
    """Parse a schedule commandline."""
    cmdln_prsr = argparse.ArgumentParser()
    cmdln_prsr.add_argument('-t', '--time', type=str, default=None,
                            help="Start Time: YYYY-mm-ddTHH:MM:SS or 'ASAP'"
                            )
    cmdln_prsr.add_argument('-p', '--project', type=str, default='0',
                            help="Project ID")
    cmdln_prsr.add_argument('-s', '--station', type=str, default=None,
                            help="Station ID")
    cmdln_prsr.add_argument('-n', '--note', type=str, default=None,
                            help="Note about session")
    cmdln_prsr.add_argument('-m', '--mockrun', action='store_true',
                            help="Mockrun")
    cmdln_prsr.add_argument('-c', '--check', action='store_true',
                            help="Check scansession sanity.")
    cmdln_prsr.add_argument('file', help='ScanSession file')
    args = cmdln_prsr.parse_args(sys.argv[1:])

    global _LOGGER
    try:
        sac = get_proj_stn_access_conf(args.project, args.station)
    except RuntimeError as e:
        _LOGGER = logging.getLogger(__name__)
        _LOGGER.error(e)
        sys.exit(1)

    stnid = args.station if args.station else sac['LCU']['stnid']

    # python logging of ilisa exec steps for session added to file
    SESLOGFILE = 'session_{}.log'.format(stnid)
    __FH = logging.FileHandler(SESLOGFILE)
    __FH.setFormatter(logging.Formatter(
        '%(asctime)s.%(msecs)03d | %(levelname)s | %(message)s',
        datefmt=ilisa.operations.DATETIMESTRFMT))
    __FH.setLevel(logging.INFO)
    _LOGGER_ROOT = logging.getLogger(__package__)
    _LOGGER_ROOT.addHandler(__FH)
    _LOGGER = logging.getLogger(__name__)

    _LOGGER.info('ilisa_obs -t{} -p{} -m{} -s{} -c{} {}'.format(
        args.time, args.project, args.mockrun, args.station, args.check,
        args.file)
    )

    with open(args.file, 'r') as f:
        scansess_in = yaml.safe_load(f)
    scansess_in['cli_start'] = args.time
    scansess_in['mockrun'] = args.mockrun
    scansess_in['projectid'] = args.project
    scansess_in['file'] = args.file
    if args.note:
        scansess_in['note'] = args.note
    scansess_in['station'] = stnid
    if args.check:
        check_scan_sess(scansess_in)
        return
    try:
        the_scansess = obs(scansess_in, sac)
    except ValueError as err:
        _LOGGER.error(err)
    _LOGGER.info('Ending scansession {}'.format(the_scansess.session_id))
    shutil.move(SESLOGFILE, os.path.join(the_scansess.get_sesspath()))


if __name__ == "__main__":
    main_cli()
