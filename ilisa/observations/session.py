"""This package provides Session class that handles sessions of observations."""
import os
import time
import datetime
import ilisa.observations.stationdriver as stationdriver
import ilisa.observations.dataIO as dataIO
import ilisa.observations.programs as programs
import ilisa.observations.modeparms as modeparms
from functools import wraps
import yaml


class Session(object):

    obslogfile = "obs.log"
    userilisadir = os.path.expanduser('~/.iLiSA/')

    def __init__(self, stnroster='stations_roster.conf', projectid=None,
                 halt_observingstate_when_finished=True):
        """Initialize Session."""

        # Setup projectmeta:
        if projectid is not None:
            projectfile = projectid + "_project.yml"
            projectfile = os.path.join(self.userilisadir, projectfile)
            with open(projectfile) as projectfilep:
                projectprofile = yaml.load(projectfilep)
            projectmeta = projectprofile['PROJECTPROFILE']
        else:
            projectmeta = {'observer': None, 'projectname': None}
        self.stnsesinfo = dataIO.StationSessionInfo(projectmeta)

        # Find accessconfig files:
        accessconffiles = []
        if stnroster is None:
            accessconffiles.append(os.path.join(self.userilisadir, 'access_config.yml'))
        else:
            stnroster = os.path.join(self.userilisadir, stnroster)
            with open(stnroster, 'r') as stnrosterf:
                for l in stnrosterf:
                    li = l.strip()
                    if li.startswith('#') or not li:
                        continue
                    accessconffile = os.path.join(self.userilisadir, li)
                    accessconffiles.append(accessconffile)

        # Initialize stationdrivers :
        self.stationdrivers = {}
        for accessconffile in accessconffiles:
            with open(accessconffile) as cfigfilep:
                accessconf = yaml.load(cfigfilep)
            stndrv = stationdriver.StationDriver(accessconf, self.stnsesinfo,
                                                 goto_observingstate_when_starting=False)
            stndrv.halt_observingstate_when_finished = halt_observingstate_when_finished
            # stationdrivers is a dict with stnid keys and corresp. stationdriver object
            self.stationdrivers[stndrv.get_stnid()] = stndrv

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

    def readlogfile(self):
        sessionlogs = []
        with open(self.obslogfile, 'r') as f:
            for sessionentry in yaml.load_all(f):
                sessionlogs.append(sessionentry)
        projectsessions = {}
        for sessionentry in sessionlogs:
            projectsessions[sessionentry['projectname']] = {}
        for sessionentry in sessionlogs:
            projectsessions[sessionentry['projectname']][sessionentry['sessionnr']] = \
                {'StartTime':sessionentry['StartTime'],
                 'Duration': sessionentry['Duration'],
                 'Command':  sessionentry['Command'],
                 'DataPaths':sessionentry['DataPaths']}


    done_logsessionbegin = False

    def logsessionbegin(self):
        """Log that the observing session is beginning."""
        if not self.done_logsessionbegin:
            with open(self.obslogfile, 'a') as ologfile:
                ologfile.write("\n---\n")
                ologfile.write("projectname: {}\n"
                               .format(self.stnsesinfo.projectmeta['projectname']))
                # ologfile.write("    SessionNr: {}\n".format(calltime.isoformat()))
        self.done_logsessionbegin = True

    def log_obs(obsf):
        @wraps(obsf)
        def logit(self, *args, **kwargs):
            self.logsessionbegin()
            calltime = datetime.datetime.utcnow()
            retval_obsf = obsf(self, *args, **kwargs)
            rettime = datetime.datetime.utcnow()
            elemind  = "    "
            listind = "  - "
            with open(self.obslogfile, 'a') as f:
                f.write(listind+"StartTime: {}\n".format(calltime.isoformat()))
                f.write(elemind+"Duration: {}\n".format(rettime - calltime))
                argsstr = ', '.join(map(repr,args))
                kwargsstr = ', '.join(['{}={!r}'.format(k, v) for k, v in kwargs.items()])
                cmdargstr = ', '.join([argsstr, kwargsstr])
                f.write(elemind+"Command: {}({})\n".format(obsf.__name__, cmdargstr))
                f.write(elemind+"DataPaths: {}\n".format(retval_obsf))
            return retval_obsf
        return logit

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

    @log_obs
    def do_bfs(self, band, duration, pointsrc, when='NOW', shutdown=True):
        """Record bfs data. Beamformed stream is capture with udpcapture on backend."""
        datafolder_urls = []
        for stndrv in self.stationdrivers.values():
            try:
                datafolder_url = stndrv.do_bfs(band, duration, pointsrc, when=when,
                                               shutdown=shutdown)
            except RuntimeError as rte:
                print("Error in do_bfs(): {}".format(rte))
                datapath = None
            datafolder_urls.append(datafolder_url)
        return datafolder_urls

    @log_obs
    def do_acc(self, band, duration, pointsrc, when='NOW'):
        """Record acc data for one of the LOFAR bands over a duration on all stations.
        """
        self._waittostart(when)
        for stndrv in self.stationdrivers.values():
            acc_url, sst_url = stndrv.do_acc(band, duration, pointsrc)
            print("Saved ACC data in folder: {}".format(acc_url))
            print("Saved SST data in folder: {}".format(sst_url))
        return [acc_url, sst_url]

    @log_obs
    def do_tbb(self, duration, band):
        """Record Transient Buffer Board (TBB) data from one of the LOFAR bands for
        duration seconds on all stations.
        """
        for stndrv in self.stationdrivers.values():
            stndrv.do_tbb(duration, band)
        return "."

    def _executeblock(self, stn, scans):
        for scan in scans:
            prg = programs.BasicObsPrograms(self.stationdrivers[stn])

            # Prepare observation arguments:
            # - Starttime
            starttime = scan['starttime']
            # - Beam
            # -- Freq
            freqbndobj = modeparms.FrequencyBand(scan['beam']['freqspec'])
            # -- Pointing
            pointsrc = scan['beam']['pointing']
            try:
                pointing = modeparms.stdPointings(pointsrc)
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

            # Collect observation parameters specified
            obsargs_in = {'starttime': starttime,
                          'freqbndobj': freqbndobj,
                          'pointsrc': pointsrc,
                          'pointing': pointing,
                          'allsky': allsky,
                          'rec_stat_type': rec_stat_type,
                          'integration': integration,
                          'duration_tot': duration_tot,
                          'rec_bfs': rec_bfs
                          }

            # If dedicated observation program chosen, set it up
            # otherwise run main obs program
            try:
                scan['obsprog']
            except KeyError:
                scan['obsprog'] = None
            if scan['obsprog'] is not None:
                obsfun, obsargs_sig = prg.getprogram(scan['obsprog'])
                # Map only args required by
                obsargs = {k: obsargs_in[k] for k in obsargs_sig}
                self.stationdrivers[stn].do_obsprog(starttime, obsfun, obsargs)
            else:
                self.stationdrivers[stn].main_scan(freqbndobj, integration, duration_tot,
                                                   pointing, pointsrc, starttime,
                                                   rec_stat_type, rec_bfs=rec_bfs,
                                                   do_acc=do_acc, allsky=allsky)

    def implement_scanschedule(self, sessionsched):
        """Implement the scan schedule dict. That is, dispatch to the
        stationdrivers to setup corresponding observations."""
        #print sessionsched
        eb = sessionsched
        # From the scan['stations'] string, generate list of stations:
        if eb['stations'] == '*' or eb['stations'] == 'ALL':
            stns = self.stationdrivers.keys()
        else:
            stns = eb['stations'].split(',')
        for stn in stns:
            try:
                self._executeblock(stn, eb['scans'])
            except RuntimeError:
                break
