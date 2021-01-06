import os
import time
import datetime

import yaml
import ilisa
import ilisa.monitorcontrol
import ilisa.monitorcontrol.directions
import ilisa.monitorcontrol.modeparms as modeparms
import ilisa.monitorcontrol.programs as programs
from ilisa.monitorcontrol.stationdriver import StationDriver, waituntil,\
    _xtract_bsx
from ilisa.monitorcontrol.stationdriver import rec_scan_start, rec_scan_stop


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

    def make_session_id(self):
        """Make a session ID based on time of creation.
        session_id has format 'sid<CT>' where <CT> is the datetime
         in the format '%Y%m%dT%H%M%S' of the time of creation.
        """
        session_id = "sid{}".format(datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%S'))
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

    def process_scansess(self, sesscans_in, session_id=None):
        """Method for parsing a station session schedule."""

        # Set the session_id to something
        if not session_id:
            try:
                session_id = sesscans_in['session_id']
            except KeyError:
                session_id = self.session_id
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
        stnid = self.stndrv.get_stnid()
        sesscans = {'session_id': session_id,
                         'projectid': projectid,
                         'station': stnid,
                         'note': note,
                         'scans': []}
        if mockrun:
            sesscans['mockrun'] = True

        starttime0 = sesscans_in['start']
        if starttime0 == 'NOW':
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
            #   if not aft:
            #       st[n] = st[n-1] + aft + gap
            #   else:
            #       st[n] = st[n-1] + dur[n-1] + gap
            if after:
                after_t = datetime.datetime.strptime(after, "%Hh%Mm%Ss")
                after_delta = datetime.timedelta(hours=after_t.hour,
                                                 minutes=after_t.minute,
                                                 seconds=after_t.second)
                starttime = starttimeprev + after_delta
            else:
                dur_delta = datetime.timedelta(seconds=duration_totprev)
                starttime = starttimeprev + dur_delta
            # - - plus a gap:
            try:
                gap = scan['gap']
            except:
                gap = "00h00m00s"
            gap_t = datetime.datetime.strptime(gap, "%Hh%Mm%Ss")
            gap_delta = datetime.timedelta(hours=gap_t.hour, minutes=gap_t.minute,
                                           seconds=gap_t.second)
            starttime += gap_delta

            # - Duration total
            duration_tot = eval(str(scan['duration_tot']))
            # Next scan use current time as previous time and current time
            starttimeprev = starttime
            duration_totprev = duration_tot

            # - Beam
            # -- Freq
            freqspec = scan['beam']['freqspec']
            # -- Pointing
            try:
                pointing = scan['beam']['pointing']
            except KeyError:
                # No pointing specified so set to None
                pointing = None
            # -- Allsky
            try:
                allsky = scan['beam']['allsky']
            except KeyError:
                allsky = False

            # - Record
            #     defaults
            try:
                scan['rec']
            except KeyError:
                scan['rec'] = []
            rec = scan['rec']
            # - Integration for rec bsx
            try:
                scan['integration']
            except KeyError:
                scan['integration'] = 1.0
            integration = eval(str(scan['integration']))
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
                               'allsky': allsky},
                          'rec': rec,
                          'integration': integration,
                          'duration_tot': duration_tot,
                          'starttime': starttime,
                          }
            obsargs_in.update({'obsprog': obsprog})
            sesscans['scans'].append(obsargs_in)
        return sesscans

    def run_scansess(self, sesscans_in, session_id=None):
        """Run a local session given a stn_ses_schedule dict. That is, dispatch to the
        stationdrivers to setup corresponding monitorcontrol."""
        sesscans = self.process_scansess(sesscans_in, session_id)
        self.projectmeta, _ = projid2meta(sesscans['projectid'])
        self.set_stn_session_id(sesscans['session_id'])
        # Set where ldata should be put after recording on LCU
        sesspath = self.get_sesspath()
        bfdsesdumpdir = self.get_bfdsesdumpdir()
        self.stndrv.scanpath = sesspath
        self.stndrv.bf_data_dir = bfdsesdumpdir
        # Boot Time handling
        nw = datetime.datetime.utcnow()
        startscantime = sesscans['scans'][0]['starttime']
        if startscantime == 'NOW':
            startscantime = nw
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
                duration_tot = scan['duration_tot']
                pointing = scan['beam']['pointing']
                starttime = scan['starttime']
                rec = scan['rec']
                integration = scan['integration']
                acc = False
                if 'acc' in rec:
                    acc = True
                    rec.remove('acc')
                bfs = False
                if 'bfs' in rec:
                    bfs = True
                    rec.remove('bfs')
                rec_type = None
                if len(rec) > 0:
                    rec_type = rec.pop()  # Should only be bsx left
                freqspec = scan['beam']['freqspec']
                freqsetup = modeparms.FreqSetup(freqspec)
                starttime = waituntil(starttime, datetime.timedelta(seconds=2))
                duration_tot, ldatinfos, ldatinfo_bfs, bfsdatapaths,\
                bfslogpaths =\
                    rec_scan_start(self.stndrv, rec_type, freqsetup,
                                   duration_tot, pointing, integration,
                                   starttime, acc=acc, bfs=bfs,
                                   destpath=sesspath)
                if not bfs and not _xtract_bsx(rec_type):
                    print('Recording for {}s'.format(duration_tot + 10))
                    time.sleep(duration_tot + 10)
                rec_scan_stop(self.stndrv, rec_type, freqsetup, pointing,
                              starttime, acc, bfs, duration_tot, ldatinfos,
                              ldatinfo_bfs, bfsdatapaths, bfslogpaths)
                scanresult = self.stndrv.scanresult
            scan['id'] = scanresult.pop('scan_id', None)
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
    is_inobsstate = stndrv.is_in_observingstate()
    print("User can observe on station {} now: {}".format(stndrv.get_stnid(),
          is_inobsstate))
    if not is_inobsstate:
        obs_allowed = stndrv.checkobservingallowed()
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
    scansess_in['start'] = args.time
    scansess_in['mockrun'] = args.mockrun
    scansess_in['projectid'] = args.project

    scnsess = ScanSession(stndrv)
    scnsess.run_scansess(scansess_in)
    args.cmd = 'obs:' + args.file


import argparse


def parse_cmdline(argv):
    """Parse a schedule commandline."""
    cmdln_prsr = argparse.ArgumentParser()
    cmdln_prsr.add_argument('-t', '--time', type=str, default='NOW',
                            help="Start Time (format: YYYY-mm-ddTHH:MM:SS)"
                            )
    cmdln_prsr.add_argument('-p', '--project', type=str, default='0', help="Project ID")
    cmdln_prsr.add_argument('-s', '--station', type=str, default=None, help="Station ID")
    cmdln_prsr.add_argument('-m', '--mockrun', action='store_true', help="Mockrun")
    cmdln_sbprsr = cmdln_prsr.add_subparsers(title='LOFAR stand alone commands',
                                             description='Select a command.',
                                             help='LOFAR commands:')

    parser_adm = cmdln_sbprsr.add_parser('adm', help="Admin")
    parser_adm.set_defaults(func=adm)
    parser_adm.add_argument('admcmd', help='Admin command')

    parser_obs = cmdln_sbprsr.add_parser('obs', help="Observe a scansession.")
    parser_obs.set_defaults(func=obs)
    parser_obs.add_argument('file', help='ScanSession file')
    args = cmdln_prsr.parse_args(argv)
    if args.time == 'NOW':
        # Set time to nearest rounded second from now:
        args.time = datetime.datetime.utcnow()
        args.time = args.time.replace(microsecond=0)
        args.time += datetime.timedelta(seconds=1)
    else:
        try:
            args.time = datetime.datetime.strptime(args.time, '%Y-%m-%dT%H:%M:%S')
        except:
            raise RuntimeError("Wrong datetime format.")
    return args


def exec_cmdline(args):
    """Run a schedule commandline."""
    projectmeta, accessfiles = projid2meta(args.project)
    if args.station is None:
        # Try to get station from accessfiles in project config file
        try:
            args.station = list(accessfiles.keys()).pop()
        except:
            raise RuntimeError("No stations found for project {}".format(args.project))
    try:
        # See if station has an access config file
        acf_name = accessfiles[args.station]
    except:
        raise RuntimeError("Station {} not found for project {}".format(args.station,
                                                                        args.project))
    userilisadir = ilisa.monitorcontrol.user_conf_dir
    acf_path = os.path.join(userilisadir, acf_name)
    with open(acf_path) as acffp:
        ac = yaml.safe_load(acffp)
    # Initialize stationdriver :
    stndrv = StationDriver(ac['LCU'], ac['DRU'], mockrun=args.mockrun)
    args.func(stndrv, args)


import sys

def main():
    """CLI to send hi-level commands to a LOFAR station.
    """
    args = parse_cmdline(sys.argv[1:])
    exec_cmdline(args)
    with open(LOGFILE, 'a') as lgf:
        if args.mockrun:
            priority_fld = 'M'
        else:
            priority_fld = '0'
        lgf.write("{} {} {} {} {}\n".format(args.time, priority_fld, args.project,
                                            args.station, args.cmd))


if __name__ == "__main__":
    main()
