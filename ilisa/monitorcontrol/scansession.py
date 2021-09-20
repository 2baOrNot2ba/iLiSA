import os
import time
import datetime

import yaml
import ilisa
import ilisa.monitorcontrol
import ilisa.monitorcontrol.directions as directions
import ilisa.monitorcontrol.modeparms as modeparms
import ilisa.monitorcontrol.programs as programs
from ilisa.monitorcontrol.stationdriver import StationDriver, waituntil,\
    _xtract_bsx, rec_scan_start, rec_scan_stop


def projid2meta(projectid):
    """Get the project metadata for project with projectid."""
    # Setup projectmeta:
    if projectid is not None:
        projectfile =  "project_"+projectid+".yml"
        projectfile = os.path.join(ilisa.monitorcontrol.user_conf_dir, projectfile)
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
    sesscans = {'session_id': session_id,
                     'projectid': projectid,
                     'station': stnid,
                     'note': note,
                     'scans': []}
    if mockrun:
        sesscans['mockrun'] = True

    if sesscans_in['cli_start']:
        # Command Line Interface requested starttime takes priority
        starttime0 = sesscans_in['cli_start']
    else:
        starttime0 = sesscans_in.get('start', 'ASAP')
    if starttime0 == 'ASAP':
        starttime0 = datetime.datetime.utcnow()

    # Initialize "previous" values for scan loop
    starttimeprev = starttime0
    duration_totprev = 0
    for scan in sesscans_in['scans']:
        # Prepare observation arguments:

        # - Starttime computed based on previous starttime
        # - - and after time:
        try:
            after = scan['after']
        except:
            after = None
        if after:
            after_delta = modeparms.hmsstr2deltatime(after)
            starttime = starttimeprev + after_delta
        else:
            dur_delta = datetime.timedelta(seconds=duration_totprev)
            starttime = starttimeprev + dur_delta

        # - Duration total
        duration_in = scan.get('duration')
        duration_tot = modeparms.hmsstr2deltatime(str(duration_in)
                                                  ).total_seconds()
        # duration_tot = eval(str(scan['duration']))
        # Next scan use current time as previous time and current time
        starttimeprev = starttime
        duration_totprev = duration_tot

        # - Beam
        beam = scan.get('beam', {})
        # -- Freq
        freqspec = beam.get('freqspec')
        # -- Pointing
        pointing_in = beam.get('pointing')
        # -- direction: alternative to pointing but can't be name
        direction = beam.get('direction')

        # -- Source name
        source = beam.get('source')
        # -- Allsky
        allsky = beam.get('allsky', False)
        # Postprocess beam to get pointing
        pointing = None
        if pointing_in:
            pointing = directions.normalizebeamctldir(pointing_in)
        elif direction:
            pointing = directions.normalizebeamctldir(direction)
        elif source:
            pointing = directions.std_pointings(source)
        # - Record
        #     defaults
        rec = scan.get('rec', [])
        _rec =  rec
        acc = False
        if 'acc' in _rec:
            acc = True
            _rec.remove('acc')
        bfs = False
        if 'bfs' in rec:
            bfs = True
            _rec.remove('bfs')
        bsx_stat = None
        if len(_rec) > 0:
            bsx_stat = _rec.pop()  # Should only be bsx left
        del _rec
        # - Integration for rec bsx
        integration = scan.get('integration', 1.0)
        if integration > duration_tot:
            raise ValueError("integration {} is longer than duration_scan {}."
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
                           'source' : source,
                           'allsky': allsky},
                      'rec': rec,
                      'acc': acc,
                      'bfs': bfs,
                      'bsx_stat': bsx_stat,
                      'integration': integration,
                      'duration': duration_tot,
                      'starttime': starttime,
                      }
        obsargs_in.update({'obsprog': obsprog})
        sesscans['scans'].append(obsargs_in)
    return sesscans


class ScanSession(object):
    """Class that runs a session on a station."""
    def __init__(self, stndrv, session_id=None):
        """Initialize Session."""
        self.stndrv = stndrv
        if session_id is None:
            session_id = self.make_session_id()
        self.session_id = session_id

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
        return self.stndrv.LOFARdataArchive

    def get_projpath(self):
        projpath = os.path.join(self.get_datastorepath(),
                                'proj{}'.format(self.projectmeta['id']))
        return projpath

    def get_sesspath(self):
        sesspath = os.path.join(self.get_projpath(),
                                'sess_{}'.format(self.get_stn_session_id()))
        return sesspath

    def get_bfdsesdumpdir(self):
        """Determine directory to where beamformed data scan should be dumped."""
        bf_data_dir = os.path.join(self.stndrv.bf_data_dir, 'proj{}'.format(
            self.projectmeta['id']), '')
        return bf_data_dir

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
        stationdrivers to setup corresponding monitorcontrol."""
        if session_id:
            self.session_id = session_id
        sesscans = process_scansess(sesscans_in, self.stndrv.get_stnid(),
                                    self.session_id)
        # Starttime handling
        startscantime = sesscans['scans'][0]['starttime']
        if (startscantime == 'ASAP' or startscantime == 'NOW'
                or not startscantime):
            startscantime = datetime.datetime.utcnow()
        if not session_id:
            self.session_id = self.make_session_id(startscantime)
        self.projectmeta, _ = projid2meta(sesscans['projectid'])
        self.set_stn_session_id(self.session_id)
        # Set where ldata should be put after recording on LCU
        sesspath = self.get_sesspath()
        bfdsesdumpdir = self.get_bfdsesdumpdir()
        self.stndrv.scanpath = sesspath
        self.stndrv.bf_data_dir = bfdsesdumpdir
        # Boot Time handling
        beaminittime = 13
        bootupstart = startscantime - datetime.timedelta(seconds=beaminittime)

        # Wait until it is time to bootup
        print("In scansess: Will start scansession @ {}".format(bootupstart))
        waituntil(bootupstart)
        self.stndrv.goto_observingstate()

        scans_done = []
        for scan in sesscans['scans']:
            freqbndobj = modeparms.FreqSetup(scan['beam']['freqspec'])
            if scan['obsprog'] is not None:
                _mockrun = False
                if not _mockrun:
                    scanresult = programs.record_obsprog(self.stndrv, scan)
                else:
                    scanresult = {}
            else:
                duration_tot = scan['duration']
                # Only pointing used not source name but it's in scan metadata
                pointing = scan['beam']['pointing']
                starttime = scan['starttime']
                #rec = scan['rec']
                acc = scan['acc']
                bfs = scan['bfs']
                bsx_stat = scan['bsx_stat']
                integration = scan['integration']
                freqspec = scan['beam']['freqspec']
                freqsetup = modeparms.FreqSetup(freqspec)
                starttime = waituntil(starttime, datetime.timedelta(seconds=2))
                duration_tot, ldatinfos, ldatinfo_bfs, bfsdatapaths,\
                bfslogpaths =\
                    rec_scan_start(self.stndrv, bsx_stat, freqsetup,
                                   duration_tot, pointing, integration,
                                   starttime, acc=acc, bfs=bfs,
                                   destpath=sesspath)
                if not bfs and not _xtract_bsx(bsx_stat):
                    print('Not recording for {}s'.format(duration_tot + 10))
                    time.sleep(duration_tot + 10)
                rec_scan_stop(self.stndrv, bsx_stat, freqsetup, pointing,
                              starttime, acc, bfs, duration_tot, ldatinfos,
                              ldatinfo_bfs, bfsdatapaths, bfslogpaths)
                scanresult = self.stndrv.scanresult
            scan['id'] = scanresult.pop('scan_id', None)
            scanresult['source'] = scan['beam']['source']
            scanpath_scdat = scanresult.pop('scanpath_scdat', None)
            self._writescanrecs(scanresult)
            print("Saved scan here: {}".format(scanpath_scdat))
            print("Finished scan @ {}".format(datetime.datetime.utcnow()))
            scans_done.append(scan)
        sesscans['scans'] = scans_done
        self.save_scansess(sesscans)


LOGFILE = "ilisa_cmds.log"


def idle(stndrv):
    """Put station into idle state."""
    stndrv.halt_observingstate()


def boot(stndrv):
    """Put station into ready to observe state."""
    stndrv.goto_observingstate()


def handback(stndrv):
    """Handback station to ILT control."""
    idle(stndrv)


def checkobs(stndrv):
    """Check if user can observe on LCU."""
    is_inobsstate = stndrv.is_inobservingstate()
    print("User can observe on station {} now: {}".format(stndrv.get_stnid(),
          is_inobsstate))
    if not is_inobsstate:
        obs_allowed = stndrv.is_observingallowed()
        if obs_allowed:
            reason = "swlevel not 3"
        else:
            reason = "Observing not allowed now"
        print("Reason: {}.".format(reason))


def adm(stndrv, args):
    """Dispatch admin commands."""
    if args.admcmd == 'boot':
        boot(stndrv)
    elif args.admcmd == 'idle':
        idle(stndrv)
    elif args.admcmd == 'handback':
        handback(stndrv)
    elif args.admcmd == 'checkobs':
        checkobs(stndrv)
    args.cmd = args.admcmd


def obs(stndrv, args):
    """Observe a scansession from ScanSes file."""
    with open(args.file, 'r') as f:
        scansess_in = yaml.safe_load(f)
    scansess_in['cli_start'] = args.time
    scansess_in['mockrun'] = args.mockrun
    scansess_in['projectid'] = args.project

    scnsess = ScanSession(stndrv)
    if args.check:
        sesscans = process_scansess(scansess_in, stndrv.get_stnid(),
                                    scnsess.session_id)
        print(yaml.dump(sesscans, default_flow_style=False))
        sys.exit()
    scnsess.run_scansess(scansess_in)
    args.cmd = 'obs:' + args.file


import argparse


def parse_cmdline(argv):
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
    cmdln_sbprsr = cmdln_prsr.add_subparsers(title='LOFAR stand alone commands',
                                             description='Select a command.',
                                             help='LOFAR commands:')

    parser_adm = cmdln_sbprsr.add_parser('adm', help="Admin")
    parser_adm.set_defaults(func=adm)
    parser_adm.add_argument('admcmd', help='Admin command')

    parser_obs = cmdln_sbprsr.add_parser('obs', help="Observe a scansession.")
    parser_obs.set_defaults(func=obs)
    parser_obs.add_argument('-c', '--check', action='store_true',
                            help="Check scansession sanity.")
    parser_obs.add_argument('file', help='ScanSession file')
    args = cmdln_prsr.parse_args(argv)
    if args.time:
        args.time = modeparms.timestr2datetime(args.time)
    return args


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
        raise RuntimeError("Station {} not found for project {}".format(stnid,
                                                                        projid))
    userilisadir = ilisa.monitorcontrol.user_conf_dir
    acf_path = os.path.join(userilisadir, acf_name)
    with open(acf_path) as acffp:
        ac = yaml.safe_load(acffp)
    return ac


def exec_cmdline(args):
    """Run a schedule commandline."""
    ac = get_proj_stn_access_conf(args.project, args.station)
    # Initialize stationdriver :
    try:
        stndrv = StationDriver(ac['LCU'], ac['DRU'], mockrun=args.mockrun)
    except AssertionError as ass_err:
        print("Error: Cannot setup {}'s LCU correctly".format(ac['LCU']['stnid']
                                                              ))
        print(ass_err)
        sys.exit(1)
    args.func(stndrv, args)


import sys


def main():
    """CLI to send hi-level commands to a LOFAR station.
    """
    args = parse_cmdline(sys.argv[1:])
    try:
        exec_cmdline(args)
    except RuntimeError as e:
        print('Exiting due to RuntimeError:', e)
        sys.exit(1)
    with open(LOGFILE, 'a') as lgf:
        if args.mockrun:
            priority_fld = 'M'
        else:
            priority_fld = '0'
        _datim_str = None
        if args.time:
            _datim_str = args.time.strftime(modeparms.DATETIMESTRFMT)
        lgf.write("{} {} {} {} {}\n".format(_datim_str, priority_fld,
                                            args.project, args.station,
                                            args.cmd))


if __name__ == "__main__":
    main()
