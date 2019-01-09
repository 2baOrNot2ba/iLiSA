"""This package provides Session class that handles observations."""
import os
import math
import time
import datetime
import ilisa.observations.stationdriver as stationdriver
import ilisa.observations.modeparms as modeparms
from functools import wraps
import yaml


class Session(object):

    obslogfile = "obs.log"

    def __init__(self, accessconffile=None, projectid=None,
                 halt_observingstate_when_finished=True):
        if accessconffile is None:
            accessconffile = os.path.expanduser('~/.iLiSA/access_config.yml')
        with open(accessconffile) as cfigfilep:
            accessconf = yaml.load(cfigfilep)
        userilisadir = os.path.expanduser('~/.iLiSA/')
        projectmeta = {'observer': None, 'projectname': None}
        if projectid is not None:
            projectfile = projectid + "_project.yml"
            projectfile = os.path.join(userilisadir, projectfile)
            with open(projectfile) as projectfilep:
                projectprofile = yaml.load(projectfilep)
            projectmeta = projectprofile['PROJECTPROFILE']
        self.observer = projectmeta['observer']
        self.project = projectmeta['projectname']
        self.stationdrivers = []
        stndrv = stationdriver.StationDriver(accessconf, projectmeta,
                                             goto_observingstate_when_starting=False)
        stndrv.halt_observingstate_when_finished = halt_observingstate_when_finished
        self.stationdrivers.append(stndrv)

    def goto_observingstate(self):
        for stndrv in self.stationdrivers:
            try:
                stndrv.goto_observingstate()
            except RuntimeError as e:
                raise RuntimeError('Could not goto observing state on station: {}'
                                   .format(e))

    def set_halt_observingstate(self):
        for stndrv in self.stationdrivers:
            stndrv.halt_observingstate_when_finished = True

    def readlogfile(self):
        sessionlogs = []
        with open(self.obslogfile, 'r') as f:
            for sessionentry in yaml.load_all(f):
                sessionlogs.append(sessionentry)
        projectsessions = {}
        for sessionentry in sessionlogs:
            projectsessions[sessionentry['ProjectName']] = {}
        for sessionentry in sessionlogs:
            projectsessions[sessionentry['ProjectName']][sessionentry['SessionNr']] = \
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
                ologfile.write("ProjectName: {}\n".format(self.project))
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
        datapaths = []
        for stndrv in self.stationdrivers:
            try:
                datapath = stndrv.do_bfs(band, duration, pointsrc, when=when,
                                         shutdown=shutdown)
            except RuntimeError as rte:
                print("Error in do_bfs(): {}".format(rte))
                datapath = None
            dataurl = "{}:{}".format(stndrv.get_stnid(), datapath)
            datapaths.append(dataurl)
        return datapaths

    @log_obs
    def do_acc(self, band, duration, pointsrc, when='NOW'):
        """Record acc data for one of the LOFAR bands over a duration on all stations.
        """
        self._waittostart(when)
        duration = int(eval(duration))
        for stndrv in self.stationdrivers:
            accdestdir = stndrv.do_acc(band, duration, pointsrc)
            print("Saved ACC data in folder: {}".format(accdestdir))
        return accdestdir

    @log_obs
    def do_bsxST(self, statistic, freqbnd, integration, duration, pointsrc,
                 when='NOW', allsky=False):
        """Records bst,sst,xst data in one of the LOFAR bands and creates a header file
        with observational settings on all stations.
        """

        duration = int(math.ceil(eval(duration)))
        frqbndobj = modeparms.FrequencyBand(freqbnd)
        self._waittostart(when)
        datapaths = []
        for stndrv in self.stationdrivers:
            if allsky and 'HBA' in frqbndobj.antsets[0]:
                stndrv.do_SEPTON(statistic, frqbndobj, integration, duration)
            else:
                datapath = None
                try:
                    datapath = stndrv.bsxST(statistic, frqbndobj, integration, duration,
                                            pointsrc)
                except RuntimeError as rte:
                    print("Error in do_bsxST(): {}".format(rte))
                dataurl = "{}:{}".format(stndrv.get_stnid(), datapath)
                datapaths.append(dataurl)
        return datapaths

    @log_obs
    def do_tbb(self, duration, band):
        """Record Transient Buffer Board (TBB) data from one of the LOFAR bands for
        duration seconds on all stations.
        """
        duration = float(eval(duration))
        for stndrv in self.stationdrivers:
            stndrv.do_tbb(duration, band)
        return "."
