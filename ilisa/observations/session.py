"""This package provides Session class that handles observations."""
import os
import sys
import math
import time
import datetime
import ilisa.observations.stationdriver as stationdriver
import ilisa.observations.modeparms as modeparms
import ilisa.observations.beamformedstreams.bfbackend as bfbackend
import ilisa.observations.dataIO as dataIO
from functools import wraps
import yaml


class Session(object):

    obslogfile = "obs.log"

    def __init__(self, accessconffile=None, projectprofile=None,
                 halt_observingstate_when_finished=True):
        if accessconffile is None:
            accessconffile = os.path.expanduser('~/.iLiSA/access_config.yml')
        with open(accessconffile) as cfigfilep:
            accessconf = yaml.load(cfigfilep)
        if projectprofile is None:
            userilisadir = os.path.expanduser('~/.iLiSA/')
            userilisadirfiles = os.listdir(userilisadir)
            for userilisafile in userilisadirfiles:
                if userilisafile.endswith('projprof.yml'):
                    projectprofile = os.path.join(userilisadir, userilisafile)
        with open(projectprofile) as projectprofilep:
            projectprofile = yaml.load(projectprofilep)
        projectmeta = projectprofile['PROJECTPROFILE']
        self.observer = projectmeta['observer']
        self.project = projectmeta['projectname']
        self.logsessionbegin()
        self.stationdrivers = []
        stndrv = stationdriver.StationDriver(accessconf, projectmeta, goto_observingstate_when_starting=False)
        stndrv.halt_observingstate_when_finished = halt_observingstate_when_finished
        self.stationdrivers.append(stndrv)

    def goto_observingstate(self):
        for stndrv in self.stationdrivers:
            try:
                stndrv.goto_observingstate()
            except RuntimeError as e:
                raise RuntimeError('Could not goto observing state on station: {}'
                                   .format(e))

    def logsessionbegin(self):
        """Log that the observing session is beginning."""
        with open(self.obslogfile, 'a') as ologfile:
            ologfile.write("\nProject: {}\n".format(self.project))

    def log_obs(obsf):
        @wraps(obsf)
        def logit(self, *args, **kwargs):
            calltime = datetime.datetime.utcnow()
            retval_obsf = obsf(self, *args, **kwargs)
            rettime = datetime.datetime.utcnow()
            with open(self.obslogfile, 'a') as f:
                f.write("{}; {}; ".format(calltime.isoformat(), rettime - calltime))
                f.write("{}{}; ".format(obsf.__name__, args, kwargs))
                f.write("{}\n".format(retval_obsf))
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

        for stndrv in self.stationdrivers:
            datapath = stndrv.do_bfs(band, duration, pointsrc, when=when,
                                     shutdown=shutdown)
        return datapath

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
