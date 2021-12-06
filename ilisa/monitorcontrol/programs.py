import sys
import os
import time
import datetime
import inspect

import ilisa.monitorcontrol.directions
import ilisa.monitorcontrol.modeparms as modeparms
import ilisa.pipelines.bfbackend as bfbackend
import ilisa.monitorcontrol.data_io as dataIO
from ilisa.monitorcontrol.stationdriver import waituntil


class ObsPrograms(object):

    def __init__(self, stationdriver):
        self.stationdriver = stationdriver

    def getprogram(self, programname):
        """Return function pointer and non-default args of a program by name.
        """
        programpointer = getattr(self, programname)
        fullargspec = inspect.getfullargspec(programpointer)
        # Remove args which have defaults:
        defargstart = None
        defarg = fullargspec.defaults
        if defarg is not None:
            defargstart = -len(defarg)
        programargs = fullargspec.args[1:defargstart]
        return programpointer, programargs

    def _streambeams_mltfreq(self, freqbndobj, pointing,
                             recDuration=float('inf'), attenuation=0,
                             DUMMYWARMUP=False):
        """Form beams with station."""
        bits = freqbndobj.bits
        if DUMMYWARMUP:
            print("Warning warnup not currently implemented")
        beamctl_cmds = []
        for bandbeamidx in range(len(freqbndobj.rcumodes)):
            _antset = freqbndobj.antsets[bandbeamidx]
            rcumode = freqbndobj.rcumodes[bandbeamidx]
            beamlets = freqbndobj.beamlets[bandbeamidx]
            subbands =  freqbndobj.subbands_spw[bandbeamidx]
            rcusel = freqbndobj.rcusel[bandbeamidx]
            beamctl_main = self.stationdriver._run_beamctl(beamlets, subbands,
                                                           rcumode, pointing,
                                                           rcusel)
            beamctl_cmds.append(beamctl_main)
        rcu_setup_cmd = self.stationdriver._rcusetup(bits, attenuation)
        return rcu_setup_cmd, beamctl_cmds

    def do_bfs_OW(self, freqbndobj, duration, pointing, bfdsesdumpdir,
                  starttime):
        """Record BeamFormed Streams (BFS) with particular beamlet allocation.
        """

        band = freqbndobj.rcubands[0]
        ###
        bits = 8  # 8
        attenuation = None
        # Subbands allocation
        if band == '10_90' or band == '30_90':
            # LBA
            lanes = (0, 1)  # (0,1)
            beamlets = '0:243'  # '0:243'
            subbands = '164:407'  # '164:407'
        elif band == '110_190':
            # HBAlo
            lanes = (0, 1, 2, 3)  # Normally (0,1,2,3) for all 4 lanes.
            beamlets = '0:487'
            subbands = '12:499'
        elif band == '210_250':
            # HBAhi
            lanes = (0, 1)
            beamlets = '0:243'
            subbands = '12:255'
        else:
            raise ValueError(
                "Wrong band: "
                "should be 10_90 (LBA), 110_190 (HBAlo) or 210_250 (HBAhi).")

        # Wait until it is time to start
        starttime_req = starttime
        warmuptime = 14
        pause = 0  # Sufficient?
        beaminittime = 13
        #self.stationdriver._waittoboot(rectime, pause)
        margin = datetime.timedelta(seconds=(warmuptime + pause + beaminittime))
        starttime = waituntil(starttime_req, margin)
        rectime = starttime

        # Necessary since fork creates multiple instances of myobs and each one
        # will call it's __del__ on completion and __del__ shutdown...
        shutdown = self.stationdriver.halt_observingstate_when_finished
        self.stationdriver.halt_observingstate_when_finished = False
        self.stationdriver.exit_check = False

        dir_bmctl = ilisa.monitorcontrol.directions.normalizebeamctldir(pointing)

        # BEGIN Dummy or hot beam start: (takes about 14sec)
        print("Running warmup beam... @ {}".format(datetime.datetime.utcnow()))
        # Setting bits is necessary:
        self.stationdriver._rcusetup(bits, attenuation)  
        self.stationdriver._run_beamctl(beamlets, subbands, band, dir_bmctl)
        self.stationdriver.stop_beam()
        # END Dummy or hot start

        print("Pause {}s after boot.".format(pause))
        time.sleep(pause)

        # Real beam start:
        print("Beam started @ UT {}".format(datetime.datetime.utcnow()))
        rcu_setup_cmds = self.stationdriver._rcusetup(bits, attenuation)
        beamctl_cmds = self.stationdriver._run_beamctl(beamlets, subbands, band,
                                                      dir_bmctl)
        rspctl_cmds = []
        beamstart = datetime.datetime.utcnow()
        timeleft = rectime - beamstart
        if timeleft.total_seconds() < 0.:
            rectime = beamstart
        print("(Beam started) Time left before recording: {}".format(
            timeleft.total_seconds()))
        bfsnametime = starttime.strftime("%Y%m%d_%H%M%S")
        obsinfo = dataIO.LDatInfo('bfs', rcu_setup_cmds, beamctl_cmds, rspctl_cmds,
                                  self.stationdriver.get_stnid())
        obsinfo.filenametime = bfsnametime

        REC = True
        if REC == True:
            port0 = self.stationdriver.bf_port0
            stnid = self.stationdriver.get_stnid()
            compress = True
            datafiles, _logfiles = bfbackend.rec_bf_streams(
                rectime, duration, lanes, band, bfdsesdumpdir, port0,
                stnid, compress)
            bfsdatapaths = []
            for lane in lanes:
                datafileguess = datafiles.pop()
                if not datafileguess:
                    _outdumpdir, _outarg, datafileguess, _dumplogname = \
                        bfbackend.bfsfilepaths(lane, rectime,
                                               modeparms.band2rcumode(band),
                                               bfdsesdumpdir, port0, stnid)
                bfsdatapaths.append(datafileguess)
        else:
            print("Not recording")
            time.sleep(duration)
        sys.stdout.flush()
        self.stationdriver.stop_beam()
        self.stationdriver.halt_observingstate_when_finished = shutdown
        return [obsinfo]


def record_obsprog(stationdriver, scan):
    """At starttime execute the observation program specified by the obsfun
    method pointer and run with arguments specified by obsargs dict.
    """
    scan_flat = dict(scan)
    del scan_flat['beam']
    for k in scan['beam'].keys():
        scan_flat[k] = scan['beam'][k]
    freqbndobj = modeparms.FreqSetup(scan['beam']['freqspec'])
    scan_flat['freqbndobj'] = freqbndobj
    prg = ObsPrograms(stationdriver)
    obsfun, obsargs_sig = prg.getprogram(scan['obsprog'])
    scan_flat['bfdsesdumpdir'] = stationdriver.bf_data_dir
    # Map only args required by
    obsargs = {k: scan_flat[k] for k in obsargs_sig}
    # Setup Calibration tables on LCU:
    CALTABLESRC = 'default'   # FIXME put this in args
    ## (Only BST uses calibration tables)
    # Choose between 'default' or 'local'
    stationdriver.set_caltable(CALTABLESRC)

    # Prepare for obs program.
    try:
        stationdriver.goto_observingstate()
    except RuntimeError as e:
        raise RuntimeError(e)

    # Run the observation program:
    obsinfolist = obsfun(**obsargs)
    # Stop program beam
    stationdriver.stop_beam()
    scan_id = None
    scanpath_scdat = None
    scanresult = {}
    if obsinfolist is not None:
        datatype = 'sop:' + scan['obsprog']
        scanrec = dataIO.ScanRecInfo()
        scanrec.set_stnid(stationdriver.get_stnid())
        scanrec.set_scanrecparms(datatype, freqbndobj.arg,
                                 scan['duration'], scan['beam']['pointing'])
        beamstarted = datetime.datetime.strptime(obsinfolist[0].filenametime,
                                                 "%Y%m%d_%H%M%S")
        scan_id = stationdriver.get_scanid(beamstarted)
        scanpath_scdat = os.path.join(stationdriver.scanpath, scan_id)
        # Add caltables used
        caltabinfos = stationdriver.get_caltableinfos(freqbndobj.rcumodes)
        # Add obsinfos to scanrecs
        for obsinfo in obsinfolist:
            # obsinfo.caltabinfos = caltabinfos
            scanrec.add_obs(obsinfo)
        scanrec.set_caltabinfos(caltabinfos)
        # Move data to archive
        stationdriver.movefromlcu(stationdriver.get_lcuDumpDir() + "/*.dat",
                                  scanpath_scdat)
        scanrec.path = scanpath_scdat
        scanresult[datatype] = scanrec
    scanresult['scan_id'] = scan_id
    scanresult['scanpath_scdat'] = scanpath_scdat
    return scanresult
