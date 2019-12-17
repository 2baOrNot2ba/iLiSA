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
                pointsrc = scan['beam']['pointing']
            except KeyError:
                # No pointing specified so set to None
                pointsrc = None
            try:
                pointing = ilisa.observations.directions.normalizebeamctldir(pointsrc)
            except KeyError:
                raise ValueError("Error: %s invalid pointing syntax".format(pointsrc))
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
                scan['rec'] = None
            bsx_type = None
            integration = 1
            rec_acc = False
            rec_bfs = False
            #     select recording a combination of acc, bfs, or bsx:
            if scan['rec'] is not None and scan['rec'] is not []:
                for chdat in scan['rec']:
                    if chdat == 'acc':
                        rec_acc = True
                    elif chdat == 'bst' or chdat == 'sst' or chdat == 'xst':
                        bsx_type = chdat
                    elif chdat == 'bfs':
                        rec_bfs = True
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
                               'pointsrc': pointsrc,
                               'pointing': pointing,
                               'allsky': allsky},
                          'integration': integration,
                          'duration_tot': duration_tot,
                          'starttime': starttime,
                          'bsx_type': bsx_type,
                          'rec_bfs': rec_bfs,
                          'rec_acc': rec_acc
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
            if scan['rec_acc']:
                scanrecs['acc'] = dataIO.ScanRecInfo()
            if scan['rec_bfs']:
                scanrecs['bfs'] = dataIO.ScanRecInfo()
            if scan['bsx_type']:
                scanrecs['bsx'] = dataIO.ScanRecInfo()
            for k in scanrecs.keys():
                scanrecs[k].set_stnid(self.stndrv.get_stnid())
            scanmeta = stationdriver.ScanMeta(sesspath, bfdsesdumpdir, scanrecs)
            if scan['obsprog'] is not None:
                scan_id, scanpath_sc, scanpath_bf \
                    = programs.record_obsprog(self.stndrv, scan, scanmeta=scanmeta)
            else:
                duration_tot = scan['duration_tot']
                pointing = scan['beam']['pointing']
                pointsrc = scan['beam']['pointsrc']
                starttime = scan['starttime']
                bsx_type = scan['bsx_type']
                integration = scan['integration']
                rec_bfs = scan['rec_bfs']
                rec_acc = scan['rec_acc']
                allsky = scan['beam']['allsky']
                scan_id, scanpath_sc, scanpath_bf \
                    = programs.record_scan(self.stndrv, freqbndobj, duration_tot,
                                                 pointing, pointsrc, starttime,
                                                 bsx_type, integration,
                                                 rec_bfs=rec_bfs, rec_acc=rec_acc,
                                                 allsky=allsky, scanmeta=scanmeta)
            print("Saved scans here: ", scanpath_sc, scanpath_bf)
            scan['id'] = scan_id
            scans_done.append(scan)
        stn_ses_sched['scans'] = scans_done
        self.save_stnsessched(stn_ses_sched)
