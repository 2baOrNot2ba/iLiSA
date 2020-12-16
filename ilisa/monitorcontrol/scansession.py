import os
import datetime

import yaml
import ilisa
import ilisa.monitorcontrol
import ilisa.monitorcontrol.directions
import ilisa.monitorcontrol.modeparms as modeparms
import ilisa.monitorcontrol.programs as programs
from ilisa.monitorcontrol.stationdriver import waituntil


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
                allsky = scan['beam']['allsky']
                scanresult = self.stndrv.record_scan(
                    freqbndobj, duration_tot, pointing, starttime, rec,
                    integration, allsky=allsky)
            scan['id'] = scanresult.pop('scan_id', None)
            scanpath_scdat = scanresult.pop('scanpath_scdat', None)
            self._writescanrecs(scanresult)
            print("Saved scan here: {}".format(scanpath_scdat))
            print("Finished scan @ {}".format(datetime.datetime.utcnow()))
            scans_done.append(scan)
        sesscans['scans'] = scans_done
        self.save_scansess(sesscans)


