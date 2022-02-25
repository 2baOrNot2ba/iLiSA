import sys
import os
import time
import datetime
import yaml
import argparse

import ilisa
import ilisa.operations
import ilisa.operations.directions as directions
import ilisa.operations.modeparms as modeparms
import ilisa.operations.programs as programs
from ilisa.operations.stationdriver import StationDriver, waituntil
from ilisa.operations.scan import still_time_fun, LScan

LOGFILE = "ilisa_cmds.log"
LATEST_SESS_DATAFILE = os.path.join(ilisa.operations.USER_CACHE_DIR,
                                    "latest_session_datafiles.txt")


def projid2meta(projectid):
    """Get the project metadata for project with projectid."""
    # Setup projectmeta:
    if projectid is not None:
        projectfile =  "project_"+projectid+".yml"
        projectfile = os.path.join(ilisa.operations.USER_CONF_DIR, projectfile)
        with open(projectfile) as projectfilep:
            projectprofile = yaml.safe_load(projectfilep)
        projectmeta = projectprofile['project']
        accessfiles = projectprofile['accessfiles']
    else:
        projectmeta = {'observer': None, 'name': None}
        accessfiles = None
    return projectmeta, accessfiles

def process_scansess(sesscans_in, stnid, session_id=None):
    """Function for parsing a station session schedule."""

    # Set the session_id to something
    if not session_id:
        session_id = sesscans_in.get('session_id')
    try:
        mockrun = sesscans_in['mockrun']
    except KeyError:
        mockrun = False
    try:
        projectid = sesscans_in['projectid']
    except KeyError:
        projectid = '0'
    try:
        note = sesscans_in['note']
    except KeyError:
        note = None

    # Initialize processed station session schedule
    sessmeta = {'session_id': session_id,
                'projectid': projectid,
                'station': stnid,
                'note': note,
                'mockrun': False,
                'cli_start': sesscans_in['cli_start'],
                'start': sesscans_in.get('start', None)
                #,'scans': []
                }
    if mockrun:
        sessmeta['mockrun'] = True
    utcnow = datetime.datetime.utcnow()
    if sessmeta['cli_start'] and sessmeta['cli_start'] != 'ASAP':
        cli_start = modeparms.timestr2datetime(sessmeta['cli_start'])
        if cli_start < utcnow:
            raise ValueError('CLI starttime already passed {} ago.'.format(
                utcnow - cli_start
            ))
    if sessmeta['start'] and sessmeta['start'] != 'ASAP':
        _start_dattim = modeparms.timestr2datetime(sessmeta['start'])
        if _start_dattim < utcnow:
            raise ValueError('Session starttime already passed')
    # Process the start time for the session
    # # cli_start overrides scan session start (if not None)
    if sessmeta['cli_start']:
        # Command Line Interface requested starttime takes priority
        sessmeta['start'] = sessmeta['cli_start']
    if not sessmeta['start']:
         sessmeta['start'] = 'ASAP'
    starttime_session = modeparms.timestr2datetime(
        modeparms.normalizetimestr(sessmeta['start']))

    def generate_scans():
        # Margin of time between two scans in seconds (fastest safe switching)
        margintime = datetime.timedelta(seconds=13.0)

        # Initialize "previous" values for scan loop
        starttimeprev = starttime_session
        duration_totprev = 0
        for scan in sesscans_in['scans']:
            # Prepare observation arguments:

            # - Starttime computed based on previous starttime
            # - - and after time:
            dur_dprev = datetime.timedelta(seconds=duration_totprev)
            after = scan.get('after')
            if after:
                after_delta = modeparms.hmsstr2deltatime(after)
                if after_delta < (dur_dprev + margintime):
                    raise ValueError(
                        'No time for next scan: after {} - dur {} < marg {}.'
                        .format(dur_dprev, after_delta, margintime))
                starttime = starttimeprev + after_delta
            elif scan == sesscans_in['scans'][0]:
                # First scan should not include margintime between scans
                starttime = starttimeprev + dur_dprev
            else:
                # All other scans should include margintime between scans
                starttime = starttimeprev + dur_dprev + margintime

            # - Duration total
            duration_in = scan.get('duration')
            duration_tot = modeparms.hmsstr2deltatime(str(duration_in)
                                                      ).total_seconds()
            # duration_tot = eval(str(scan['duration']))
            file_dur = scan.get('file_dur')
            # Next scan use current time as previous time and current time
            starttimeprev = starttime
            duration_totprev = duration_tot

            # - Source name
            source = scan.get('source')

            # - Beam
            beam = scan.get('beam', {})
            # -- Freq
            freqspec = beam.get('freqspec')
            # -- Pointing
            pointing = beam.get('pointing')
            # -- direction: alternative to pointing but can't be name
            direction = beam.get('direction')

            # -- Allsky
            allsky = beam.get('allsky', False)

            # Postprocess beam to get direction
            if not direction:
                if pointing:
                    direction = directions.normalizebeamctldir(pointing)
                elif source:
                    direction = directions.std_pointings(source)
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
            try:
                obsprog = scan['obsprog']
            except KeyError:
                obsprog = None

            # Collect observation parameters specified
            obsargs_in = {'beam':
                              {'freqspec': freqspec,
                               'pointing': pointing,
                               'direction': direction,
                               'allsky': allsky},
                          'rec': rec,
                          'acc': acc,
                          'bfs': bfs,
                          'bsx_stat': bsx_stat,
                          'integration': integration,
                          'duration': duration_tot,
                          'file_dur': file_dur,
                          'starttime': starttime,
                          'source': source
                          }
            obsargs_in.update({'obsprog': obsprog})
            #sessmeta['scans'].append(obsargs_in)
            yield obsargs_in
    return sessmeta, generate_scans()


class ScanSession(object):
    """Class that runs a session on a station."""
    def __init__(self, stndrv, session_id=None):
        """Initialize Session."""
        self.stndrv = stndrv
        if session_id is None:
            session_id = self.make_session_id()
        self.session_id = session_id
        if os.path.exists(LATEST_SESS_DATAFILE):
            os.remove(LATEST_SESS_DATAFILE)

    def set_stn_session_id(self, parent_session_id):
        self.stn_sess_id = '{}_{}'.format(parent_session_id, self.stndrv.get_stnid())

    def get_stn_session_id(self):
        return self.stn_sess_id

    def make_session_id(self, ref_dattim=datetime.datetime.utcnow()):
        """Make a session ID based on time of creation.
        session_id has format 'sid<CT>' where <CT> is the datetime
        in the format '%Y%m%dT%H%M%S' of the time of creation.
        """
        session_id = "sid{}".format(ref_dattim.strftime('%Y%m%dT%H%M%S'))
        return session_id

    def get_datastorepath(self):
        return self.stndrv.dru_data_root

    def get_projpath(self):
        projpath = os.path.join(self.get_datastorepath(), 'Projects',
                                'proj{}'.format(self.projectmeta['id']))
        return projpath

    def get_sesspath(self):
        sesspath = os.path.join(self.get_projpath(),
                                'sess_{}'.format(self.get_stn_session_id()))
        return sesspath

    def save_scansess(self, sched):
        sesspath = self.get_sesspath()
        ses_sched_file = os.path.join(sesspath, 'ScanSess.yml')
        with open(ses_sched_file, 'w') as f:
            yaml.dump(sched, f, explicit_start=True)

    def _writescanrecs(self, scanrecs):
        """Write the  scanrec for each recorded ldat."""
        for ldat in scanrecs.keys():
            try:
                scanrecpath = scanrecs[ldat].get_scanrecpath()
            except (AttributeError, KeyError):
                scanrecpath = None
            if scanrecpath:
                scanrecs[ldat].write(scanrecpath)

    def run_scansess(self, sesscans_in, session_id=None):
        """Run a local session given a stn_ses_schedule dict. That is, dispatch to the
        stationdrivers to setup corresponding operations."""
        if session_id:
            self.session_id = session_id
        sessmeta, scans_iter = process_scansess(sesscans_in,
                                                self.stndrv.get_stnid(),
                                                self.session_id)
        # Starttime handling
        session_start = modeparms.timestr2datetime(sessmeta['start'])
        if not session_id:
            self.session_id = self.make_session_id(session_start)
        self.projectmeta, _ = projid2meta(sessmeta['projectid'])
        self.set_stn_session_id(self.session_id)
        # Set where ldata should be put after recording on LCU
        sesspath = self.get_sesspath()
        bfdsesdumpdir = 'proj{}'.format(self.projectmeta['id'])
        # Boot Time handling
        dt2beamctl = datetime.timedelta(
            seconds=self.stndrv._beam_time2startup_hint())

        # Wait until it is time to bootup
        waituntil(session_start, dt2beamctl)
        print("(scansession started @ {})".format(datetime.datetime.utcnow()))

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
                                 'source': scan['source']}
                acc = scan['acc']
                bfs = scan['bfs']
                bsx_stat = scan['bsx_stat']
                integration = scan['integration']
                freqspec = scan['beam']['freqspec']
                freqsetup = modeparms.FreqSetup(freqspec)

                starttime = scan['starttime']
                duration_tot = scan['duration']
                stoptime = starttime + datetime.timedelta(
                    seconds=int(duration_tot))
                print('Will stop @', stoptime)
                stop_cond = still_time_fun(stoptime)
                margin_scan_start = datetime.timedelta(seconds=2)
                startedtime = waituntil(starttime, margin_scan_start)
                # Initialize LScan
                lscan = LScan(self.stndrv, bsx_stat, freqsetup,
                              duration_tot, pointing_spec, integration,
                              starttime, acc=acc, bfs=bfs,
                              destpath=sesspath, destpath_bfs=bfdsesdumpdir,
                              file_dur=scan['file_dur'])
                subscan = iter(lscan)
                # Start the subscan
                next(subscan)
                while stop_cond():
                    try:
                        print('IN TOP SUBSCAN LOOP', datetime.datetime.utcnow())
                        _ = subscan.send(stop_cond())
                    except StopIteration:
                        break
                lscan.close()

                if not bfs and not modeparms._xtract_bsx(bsx_stat) and not acc:
                    _sleepfor = (datetime.timedelta(seconds=duration_tot)
                                 + margin_scan_start)
                    print('Not recording for {}s'.format(_sleepfor))
                    time.sleep(_sleepfor)
                scanresult = lscan.scanresult
            scan['id'] = scanresult.get('scan_id', None)
            scanpath_scdat = scanresult.get('scanpath_scdat', None)
            print("Saved scan here: {}".format(scanpath_scdat))
            # Append latest scanrec paths file
            if 'bsx' in lscan.scanresult['rec']:
                with open(LATEST_SESS_DATAFILE, 'a') as f:
                    f.write(lscan.scanresult['bsx'].scanrecpath)
                    f.write('\n')
            scan_ended_at = datetime.datetime.utcnow()
            duration_actual = scan_ended_at - startedtime
            print("Finished scan @", scan_ended_at,
                  "Request dur=", duration_tot,
                  "Actual dur=", duration_actual)
            scans_done.append(scan)
        sessmeta['scans'] = scans_done
        self.save_scansess(sessmeta)


def get_proj_stn_access_conf(projid, stnid):
    """Get access conf for stnid as per projid"""
    projectmeta, accessfiles = projid2meta(projid)
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
                           .format(stnid, projid))
    userilisadir = ilisa.operations.USER_CONF_DIR
    acf_path = os.path.join(userilisadir, acf_name)
    with open(acf_path) as acffp:
        ac = yaml.safe_load(acffp)
    return ac


def obs(stndrv, args):
    """Observe a scansession from ScanSes file."""
    with open(args.file, 'r') as f:
        scansess_in = yaml.safe_load(f)
    scansess_in['cli_start'] = args.time
    scansess_in['mockrun'] = args.mockrun
    scansess_in['projectid'] = args.project

    scnsess = ScanSession(stndrv)
    if args.check:
        sessmeta, scans_obsargs = process_scansess(scansess_in,
                                                   stndrv.get_stnid(),
                                                   scnsess.session_id)
        print(yaml.dump(sessmeta, default_flow_style=False))
        sessscans={'scans': []}
        try:
            for scan_obsargs in scans_obsargs:
                sessscans['scans'].append(scan_obsargs)
        except ValueError as ve:
            print("ValueError:", ve)
        else:
            print(yaml.dump(sessscans, default_flow_style=False))
        sys.exit()
    scnsess.run_scansess(scansess_in)
    args.cmd = 'obs:' + args.file


def main_cli():
    """\
    CLI to send hi-level commands to a LOFAR station.
    """
    """Parse a schedule commandline."""
    cmdln_prsr = argparse.ArgumentParser()
    cmdln_prsr.add_argument('-t', '--time', type=str, default=None,
                            help="Start Time (format: YYYY-mm-ddTHH:MM:SS)"
                            )
    cmdln_prsr.add_argument('-p', '--project', type=str, default='0',
                            help="Project ID")
    cmdln_prsr.add_argument('-s', '--station', type=str, default=None,
                            help="Station ID")
    cmdln_prsr.add_argument('-m', '--mockrun', action='store_true',
                            help="Mockrun")

    cmdln_prsr.add_argument('-c', '--check', action='store_true',
                            help="Check scansession sanity.")
    cmdln_prsr.add_argument('file', help='ScanSession file')
    args = cmdln_prsr.parse_args(sys.argv[1:])

    issued_at = datetime.datetime.utcnow().isoformat(timespec='seconds')
    try:
        # Run an observation based on ScanSes spec file
        ac = get_proj_stn_access_conf(args.project, args.station)
        # Initialize stationdriver :
        try:
            stndrv = StationDriver(ac['LCU'], ac['DRU'], mockrun=args.mockrun)
        except AssertionError as ass_err:
            print("Error: Cannot setup {}'s LCU correctly"
                  .format(ac['LCU']['stnid']))
            print(ass_err)
            sys.exit(1)
        obs(stndrv, args)
    except RuntimeError as e:
        print('Exiting due to RuntimeError:', e)
        sys.exit(1)
    with open(LOGFILE, 'a') as lgf:
        if args.mockrun:
            priority_fld = 'M'
        else:
            priority_fld = '0'
        lgf.write(f"{issued_at} {args.time} {priority_fld} {args.project}"
                  f" {stndrv.get_stnid()} {args.cmd}\n")


if __name__ == "__main__":
    main_cli()
