"""This package provides Session class that handles sessions of observations."""
import os
import time
import datetime
import ilisa
import ilisa.observations.stationdriver as stationdriver
import yaml

def create_session_id():
    """Create a session ID. session_id has format 'sid<CT>' where <CT> is the datetime
     in the format '%Y%m%dT%H%M%S' of the time of creation."""
    session_id = "sid{}".format(datetime.datetime.utcnow().strftime('%Y%m%dT%H%M%S'))
    return session_id

class Session(object):

    obslogfile = "obs.log"

    def __init__(self, stnroster='stations_roster.conf', session_id=None, mockrun=False,
                 halt_observingstate_when_finished=True):
        """Initialize Session."""

        # Find accessconfig files:
        accessconffiles = []
        if stnroster is None:
            accessconffiles.append(os.path.join(ilisa.user_conf_dir, 'access_config.yml'))
        else:
            stnroster = os.path.join(ilisa.user_conf_dir, stnroster)
            with open(stnroster, 'r') as stnrosterf:
                for l in stnrosterf:
                    li = l.strip()
                    if li.startswith('#') or not li:
                        continue
                    accessconffile = os.path.join(ilisa.user_conf_dir, li)
                    accessconffiles.append(accessconffile)

        # Initialize stationdrivers :
        self.stationdrivers = {}
        for accessconffile in accessconffiles:
            with open(accessconffile) as cfigfilep:
                accessconf = yaml.load(cfigfilep)
            stndrv = stationdriver.StationDriver(accessconf, mockrun=mockrun,
                                                 goto_observingstate_when_starting=False)
            stndrv.halt_observingstate_when_finished = halt_observingstate_when_finished
            # stationdrivers is a dict with stnid keys and corresp. stationdriver object
            self.stationdrivers[stndrv.get_stnid()] = stndrv

        self.set_session_id(session_id)

    def set_session_id(self, session_id=None):
        """Set the session ID. If session_id is None then a session ID based on time of
        creation will be used."""
        if session_id is None:
            session_id = create_session_id()
        self.session_id = session_id

    def goto_observingstate(self):
        for stndrv in self.stationdrivers.values():
            try:
                stndrv.goto_observingstate()
            except RuntimeError as e:
                raise RuntimeError('Could not goto observing state on station: {}'
                                   .format(e))

    def set_halt_observingstate(self):
        for stndrv in self.stationdrivers.values():
            stndrv.halt_observingstate_when_finished = True

    def log(self, com_ses_sched):
        with open(self.obslogfile,'a') as f:
            yaml.dump(com_ses_sched, f)

    def _waittostart(self, when):
        if when != 'NOW':
            now = datetime.datetime.utcnow()
            try:
                st = datetime.datetime.strptime(when, "%Y-%m-%dT%H:%M:%S")
            except ValueError:
                raise ValueError('Start time format not understood')
            maxminsetuptime = datetime.timedelta(
                seconds=0)
            d = (st - maxminsetuptime) - now
            timeuntilboot = d.total_seconds()
            if timeuntilboot < 0.:
                timeuntilboot = 0
            print("Will start after " + str(
                timeuntilboot) + " seconds...")
            time.sleep(timeuntilboot)

    def implement_scanschedule(self, sessionsched):
        """Implement the scan schedule dict. That is, dispatch to the
        stationdrivers to setup corresponding observations."""
        stn_ses_sched = self.split_stn_sched(sessionsched)
        self.log(sessionsched)
        for stn in stn_ses_sched.keys():
            self.stationdrivers[stn].run_lcl_session(stn_ses_sched[stn], self.session_id)

    def split_stn_sched(self, ses_sched_in):
        """Is a method for splitting a common session schedule by station."""
        def _parse_stations(stns_str):
            # From the stations string, generate list of stations:
            if stns_str == '*' or stns_str == 'ALL':
                stns = self.stationdrivers.keys()
            else:
                stns = stns_str.split(',')
            return stns

        nonscanfields = []
        for field in ses_sched_in.keys():
            if field != 'scans':
                nonscanfields.append(field)

        # Initialize list of station session schedules
        stn_ses_sched = {}
        stns = set()
        for scan in ses_sched_in['scans']:
            scan_stns = _parse_stations(scan['stations'])
            stns.update(scan_stns)
        scansnostations = []
        for scan in ses_sched_in['scans']:
            del scan['stations']
            scansnostations.append(scan)
        for stn in stns:
            stn_ses_sched[stn] = {field: ses_sched_in[field] for field in nonscanfields}
            stn_ses_sched[stn]['scans'] = scansnostations
        return stn_ses_sched
