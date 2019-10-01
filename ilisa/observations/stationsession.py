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
        os.makedirs(sesspath)
        ses_sched_file = os.path.join(sesspath, 'STN_SESSION.yml')
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

        # Initialize processed station session schedule
        stnid = self.stndrv.get_stnid()
        stn_ses_sched = {'session_id': session_id,
                         'projectid': projectid,
                         'station': stnid,
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
                pointsrc = scan['beam']['pointing']
            except KeyError:
                # No pointing specified so set to None
                pointsrc = None
            try:
                pointing = ilisa.observations.directions.stdPointings(pointsrc)
            except KeyError:
                try:
                    phi, theta, ref = pointsrc.split(',', 3)
                    # FIXME:  (not always going to be correct)
                    pointing = pointsrc
                except ValueError:
                    raise ValueError(
                        "Error: %s invalid pointing syntax".format(pointsrc))
            # -- Allsky
            try:
                allsky = scan['beam']['allsky']
            except KeyError:
                allsky = False
            # - Record statistics
            try:
                scan['rec_stat']
            except KeyError:
                scan['rec_stat'] = None
            if scan['rec_stat'] is not None:
                rec_stat_type = scan['rec_stat']['type']
                integration = eval(str(scan['rec_stat']['integration']))
            else:
                rec_stat_type = None
                integration = 1
            # - Duration total
            duration_tot = eval(str(scan['duration_tot']))
            if integration > duration_tot:
                raise ValueError("integration {} is longer than duration_scan {}."
                                 .format(integration, duration_tot))
            # - ACC
            try:
                do_acc = scan['acc']
            except KeyError:
                do_acc = False
            # - BFS
            try:
                rec_bfs = scan['rec_bfs']
            except KeyError:
                rec_bfs = False

            # If dedicated observation program chosen, set it up
            # otherwise run main obs program
            try:
                obsprog = scan['obsprog']
            except KeyError:
                obsprog = None

            # Collect observation parameters specified
            obsargs_in = {'beam':
                              {'freqspec': freqspec,
                               'pointsrc': pointsrc,
                               'pointing': pointing,
                               'allsky': allsky},
                          'integration': integration,
                          'duration_tot': duration_tot,
                          'starttime': starttime,
                          'rec_stat_type': rec_stat_type,
                          'rec_bfs': rec_bfs,
                          'do_acc': do_acc
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
        self.save_stnsessched(stn_ses_sched)
        sesspath = self.get_sesspath()
        bfdsesdumpdir = self.get_bfdsesdumpdir()
        # Boot Time handling
        nw = datetime.datetime.utcnow()
        starttime = stn_ses_sched['scans'][0]['starttime']
        if starttime == 'NOW':
            starttime = nw
        bootupstart = starttime - datetime.timedelta(seconds=10)
        # Wait until it is time to bootup
        print("In stnsess: Wait until it is time to bootup")
        st = self.stndrv._waittoboot(bootupstart.strftime("%Y-%m-%dT%H:%M:%S"))
        for scan in stn_ses_sched['scans']:
            freqbndobj = modeparms.FrequencyBand(scan['beam']['freqspec'])
            scanrecs={}
            if scan['do_acc']:
                scanrecs['acc'] = dataIO.ScanRecInfo(self.projectmeta)
            if scan['rec_bfs']:
                scanrecs['bfs'] = dataIO.ScanRecInfo(self.projectmeta)
            if scan['rec_stat_type']:
                scanrecs['bsx'] = dataIO.ScanRecInfo(self.projectmeta)
            for k in scanrecs.keys():
                scanrecs[k].set_stnid(self.stndrv.get_stnid())
            scanmeta = stationdriver.ScanMeta(sesspath, bfdsesdumpdir, scanrecs)
            if scan['obsprog'] is not None:
                self.stndrv.do_obsprog(scan, scanmeta=scanmeta)
            else:
                integration = scan['integration']
                duration_tot = scan['duration_tot']
                pointing = scan['beam']['pointing']
                pointsrc = scan['beam']['pointsrc']
                starttime = scan['starttime']
                rec_stat_type = scan['rec_stat_type']
                rec_bfs = scan['rec_bfs']
                do_acc = scan['do_acc']
                allsky = scan['beam']['allsky']
                self.stndrv.main_scan(freqbndobj, integration, duration_tot, pointing,
                                      pointsrc, starttime, rec_stat_type,
                                      rec_bfs=rec_bfs, do_acc=do_acc, allsky=allsky,
                                      scanmeta=scanmeta)
