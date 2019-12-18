import os
import datetime
import yaml
import ilisa
import ilisa.observations
import ilisa.observations.dataIO as dataIO
import ilisa.observations.directions
import ilisa.observations.session as session
import ilisa.observations.modeparms as modeparms
import ilisa.observations.stationdriver as stationdriver
import ilisa.observations.programs as programs


def projid2meta(projectid):
    """Get the project metadata for project with projectid."""
    # Setup projectmeta:
    if projectid is not None:
        projectfile =  "project_"+projectid+".yml"
        projectfile = os.path.join(ilisa.observations.user_conf_dir, projectfile)
        with open(projectfile) as projectfilep:
            projectprofile = yaml.load(projectfilep)
        projectmeta = projectprofile['project']
        accessfiles = projectprofile['accessfiles']
    else:
        projectmeta = {'observer': None, 'name': None}
        accessfiles = None
    return projectmeta, accessfiles


class StationSession(object):
    """Class that runs a session on a station."""
    def __init__(self, ac_lcu, ac_dru, session_id=None,
                 mockrun=False, halt_observingstate_when_finished=True):
        """Initialize Session."""

        self.LOFARdataArchive = ac_dru['LOFARdataArchive']

        # Initialize stationdriver :
        self.stndrv = stationdriver.StationDriver(ac_lcu, ac_dru, mockrun=mockrun,
                                                  goto_observingstate_when_starting=False)
        self.stndrv.halt_observingstate_when_finished = halt_observingstate_when_finished
        if session_id is None:
            session_id = session.create_session_id()
        self.session_id = session_id

    def set_stn_session_id(self, parent_session_id):
        self.stn_sess_id = '{}_{}'.format(parent_session_id, self.stndrv.get_stnid())

    def get_stn_session_id(self):
        return self.stn_sess_id

    def get_datastorepath(self):
        return self.LOFARdataArchive

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

    def save_stnsessched(self, sched):
        sesspath = self.get_sesspath()
        ses_sched_file = os.path.join(sesspath, 'StnSess.yml')
        with open(ses_sched_file, 'w') as f:
            yaml.dump(sched, f, explicit_start=True)

    def make_session_id(self):
        """Make a session ID based on time of creation."""
        session_id = "sid{}".format(datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%S'))
        return session_id

    def process_stn_ses_sched(self, stn_ses_sched_in, session_id=None):
        """Method for parsing a station session schedule."""

        # Set the session_id to something
        if not session_id:
            try:
                session_id = stn_ses_sched_in['session_id']
            except KeyError:
                session_id = self.make_session_id()
        try:
            mockrun = stn_ses_sched_in['mockrun']
        except KeyError:
            mockrun = False
        try:
            projectid = stn_ses_sched_in['projectid']
        except KeyError:
            projectid = '0'
        try:
            note = stn_ses_sched_in['note']
        except KeyError:
            note = None

        # Initialize processed station session schedule
        stnid = self.stndrv.get_stnid()
        stn_ses_sched = {'session_id': session_id,
                         'projectid': projectid,
                         'station': stnid,
                         'note': note,
                         'scans': []}
        if mockrun:
            stn_ses_sched['mockrun'] = True

        starttime0 = stn_ses_sched_in['scans'][0]['starttime']
        if starttime0 == 'NOW':
            starttime0 = datetime.datetime.utcnow()

        for scan in stn_ses_sched_in['scans']:
            # Prepare observation arguments:
            # - Starttime
            starttime = scan['starttime']
            #   Check for deltatime:
            if starttime != 'NOW' and type(starttime) is not datetime.datetime:
                if type(starttime) is str:
                    starttime = starttime0 + datetime.timedelta(seconds=
                                                                int(eval(starttime)))
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
            # - Duration total
            duration_tot = eval(str(scan['duration_tot']))
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
            stn_ses_sched['scans'].append(obsargs_in)
        return stn_ses_sched

    def run_lcl_sched(self, stn_ses_sched_in, session_id=None):
        """Run a local session given a stn_ses_schedule dict. That is, dispatch to the
        stationdrivers to setup corresponding observations."""
        stn_ses_sched = self.process_stn_ses_sched(stn_ses_sched_in, session_id)
        self.projectmeta, _ = projid2meta(stn_ses_sched['projectid'])
        self.set_stn_session_id(stn_ses_sched['session_id'])
        sesspath = self.get_sesspath()
        bfdsesdumpdir = self.get_bfdsesdumpdir()
        self.stndrv.scanpath = sesspath
        # Boot Time handling
        nw = datetime.datetime.utcnow()
        startscantime = stn_ses_sched['scans'][0]['starttime']
        if startscantime == 'NOW':
            startscantime = nw
        bootupstart = startscantime - datetime.timedelta(seconds=10)
        # Wait until it is time to bootup
        print("In stnsess: Wait until it is time to bootup")
        self.stndrv._waittoboot(bootupstart)
        scans_done = []
        for scan in stn_ses_sched['scans']:
            freqbndobj = modeparms.FrequencyBand(scan['beam']['freqspec'])
            scanrecs={}
            scanrecs['acc'] = dataIO.ScanRecInfo()
            scanrecs['bfs'] = dataIO.ScanRecInfo()
            scanrecs['bsx'] = dataIO.ScanRecInfo()
            for k in scanrecs.keys():
                scanrecs[k].set_stnid(self.stndrv.get_stnid())
            scanmeta = stationdriver.ScanMeta(sesspath, bfdsesdumpdir, scanrecs)
            if scan['obsprog'] is not None:
                programs.record_obsprog(self.stndrv, scan, scanmeta=scanmeta)
            else:
                duration_tot = scan['duration_tot']
                pointing = scan['beam']['pointing']
                starttime = scan['starttime']
                rec = scan['rec']
                integration = scan['integration']
                allsky = scan['beam']['allsky']
                programs.record_scan(self.stndrv, freqbndobj, duration_tot,
                                     pointing, starttime, rec, integration,
                                     allsky=allsky, scanmeta=scanmeta)
            # Write scanrecinfo files and ldat headers
            scanrecpath = scanmeta.write_scanrecs()
            print("Saved scan here: {}".format(scanrecpath))
            scan['id'] = scanmeta.scan_id
            scans_done.append(scan)
        stn_ses_sched['scans'] = scans_done
        self.save_stnsessched(stn_ses_sched)
