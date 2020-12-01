import sys
import os
import time
import datetime
import shutil
import warnings
import inspect

import ilisa.observations.directions
import ilisa.observations.modeparms as modeparms
import ilisa.pipelines.bfbackend as bfbackend
import ilisa.observations.dataIO as dataIO
from ilisa.observations.stationdriver import waituntil


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
            subbands =  freqbndobj.sb_range[bandbeamidx]
            rcusel = freqbndobj.rcusel[bandbeamidx]
            beamctl_main = self.stationdriver._run_beamctl(beamlets, subbands,
                                                           rcumode, pointing,
                                                           rcusel)
            beamctl_cmds.append(beamctl_main)
        rcu_setup_cmd = self.stationdriver._rcusetup(bits, attenuation)
        return rcu_setup_cmd, beamctl_cmds

    def do_bfs_OW(self, freqbndobj, duration_tot, pointing, bfdsesdumpdir,
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

        dir_bmctl = ilisa.observations.directions.normalizebeamctldir(pointing)

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
        rcu_setup_cmd = self.stationdriver._rcusetup(bits, attenuation)
        beamctl_cmd = self.stationdriver._run_beamctl(beamlets, subbands, band,
                                                      dir_bmctl)
        beamstart = datetime.datetime.utcnow()
        timeleft = rectime - beamstart
        if timeleft.total_seconds() < 0.:
            rectime = beamstart
        print("(Beam started) Time left before recording: {}".format(
            timeleft.total_seconds()))
        bfsnametime = starttime.strftime("%Y%m%d_%H%M%S")
        obsinfo = dataIO.LDatInfo('bfs', bfsnametime, beamctl_cmd,
                                  rcu_setup_cmd, caltabinfos="")

        REC = True
        if REC == True:
            port0 = self.stationdriver.bf_port0
            stnid = self.stationdriver.get_stnid()
            datafiles, _logfiles = bfbackend.rec_bf_streams(
                rectime, duration_tot, lanes, band, bfdsesdumpdir, port0,
                stnid)
            bfsdatapaths = []
            for lane in lanes:
                datafileguess = datafiles.pop()
                if not datafileguess:
                    _outdumpdir, _outarg, datafileguess, _dumplogname = \
                        bfbackend.bfsfilepaths(lane, rectime, band,
                                               bfdsesdumpdir, port0, stnid)
                bfsdatapaths.append(datafileguess)
        else:
            print("Not recording")
            time.sleep(duration_tot)
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
    freqbndobj = modeparms.FrequencyBand(scan['beam']['freqspec'])
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
    scanresult = {}
    if obsinfolist is not None:
        datatype = 'sop:' + scan['obsprog']
        scanrec = dataIO.ScanRecInfo()
        scanrec.set_stnid(stationdriver.get_stnid())
        scanrec.set_scanrecparms(datatype, freqbndobj.arg,
                                 scan['duration_tot'], scan['beam']['pointing'],
                                 allsky=False)
        beamstarted = datetime.datetime.strptime(obsinfolist[0].filenametime,
                                                 "%Y%m%d_%H%M%S")
        scan_id = stationdriver.get_scanid(beamstarted)
        scanpath_scdat = os.path.join(stationdriver.scanpath, scan_id)
        # Add caltables used
        caltabinfos = []
        for rcumode in freqbndobj.rcumodes:
            caltabinfo = stationdriver.get_caltableinfo(rcumode)
            caltabinfos.append(caltabinfo)
        # Add obsinfos to scanrecs
        for obsinfo in obsinfolist:
            obsinfo.caltabinfos = caltabinfos
            scanrec.add_obs(obsinfo)
        # Move data to archive
        stationdriver.movefromlcu(stationdriver.get_lcuDumpDir() + "/*.dat",
                                  scanpath_scdat)
        scanrec.path = scanpath_scdat
        scanresult[datatype] = scanrec
    scanresult['scan_id'] = scan_id
    scanresult['scanpath_scdat'] = scanpath_scdat
    return scanresult


def record_scan(stationdriver, freqbndobj, duration_tot, pointing,
                starttime='NOW', rec=(None,), integration=1.0, allsky=False,
                duration_frq=None):
    """Run a generic scan.
    
    Note
    ----


    Parameters
    ----------
        stationdriver: StationDriver()
            A StationDriver object.
        freqbndobj: FrequencyBand
            Frequency specification of scan.
        integration: float
            Integration time in seconds.
        duration_tot: float
            Total duration of scan in seconds.
        pointing: str
            Pointing source direction of scan. Can be a source name or a
            beamctl direction str.
        starttime: str
            The time at which scan should start.
        rec: tuple
            The types of LOFAR data to record.
            This is a tuple of up to three of the following data types:
                'acc', 'bfs', 'bst', 'sst', 'xst'
            of which 3 latter types are mutual exclusive.

            'acc': Record autocovariance-cubes (ACC) files. The actual total
                duration will at most be duration_tot. (Usually it will be
                shorter so it fits within the cadence of whole ACC aquisition,
                which is 512+7=519 seconds). ACC files are the covariance of
                all array elements with each as a function of subband.

            'bfs': Record BeamFormed Stream data.

            'bst': Record Beamlet STatistics data.

            'sst': Record Subband STatistics data.

            'xst': Record Xrosslet STatistics data.

            None: No recording of data.
        duration_frq: float
            Duration in seconds of a sampled frequency within a frequency sweep
            scan.
        allsky: bool, optional
            HBA allsky mode.
    """
    rec_acc = False; rec_bfs = False; bsx_type = None
    if rec:
        for recreq in rec:
            if recreq == 'acc':
                rec_acc = True
            elif recreq == 'bfs':
                rec_bfs = True
            elif recreq == 'bst' or recreq == 'sst' or recreq == 'xst':
                bsx_type = recreq
    arraytype = freqbndobj.antsets[-1][:3]
    # Mode logic
    todo_tof = False  # This is the case when arraytype=='LBA'
    if allsky and arraytype == 'HBA':
        if pointing is not None:
            raise ValueError('hba-allsky and beam cannot run simultaneously.')
        if bsx_type == 'bst':
            raise ValueError('bst and hba-allsky cannot be combined.')
        if rec_acc:
            raise ValueError('acc and hba-allsky cannot be combined.')
        if rec_bfs:
            raise ValueError('bfs and hba-allsky cannot be combined.')
        todo_tof = True
    if not allsky and pointing is None:
        pointing = 'Z'

    stnid = stationdriver.get_stnid()

    # Initialize scanresult
    scanresult = {'rec': []}

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

    duration_tot_req = duration_tot
    band = freqbndobj.rcubands[0]
    rcumode = freqbndobj.rcumodes[0]

    if todo_tof:
        septonconf = stationdriver.setup_tof()
    else:
        septonconf = None

    if rec_acc:
        scanresult['rec'].append('acc')
        scanresult['acc'] = dataIO.ScanRecInfo()
        scanresult['acc'].set_stnid(stnid)
        # Also duration of ACC sweep since each sb is 1 second.
        dur1acc = modeparms.TotNrOfsb  # Duration of one ACC
        interv2accs = 7  # time between end of one ACC and start of next one
        acc_cadence = dur1acc+interv2accs  # =519s time between start of 2 ACCs
        (nraccs, timrest) = divmod(duration_tot_req, acc_cadence)
        if timrest > dur1acc:
            nraccs += 1
        duration_tot = nraccs*acc_cadence-interv2accs
        if duration_tot != duration_tot_req:
            print("""Note: will use total duration {}s to fit with ACC
                  cadence.""".format(duration_tot))
        stationdriver.acc_mode(enable=True, mock_dur=duration_tot)

    # Wait until it is time to start
    starttime_req = starttime
    starttime = waituntil(starttime_req)

    # Necessary since fork creates multiple instances of myobs and each one
    # will call it's __del__ on completion and __del__ shutdown...
    shutdown = stationdriver.halt_observingstate_when_finished
    stationdriver.halt_observingstate_when_finished = False
    stationdriver.exit_check = False

    caltabinfos = [""]
    # Get metadata about caltables to be used
    if not allsky:
        for rcumode in freqbndobj.rcumodes:
            caltabinfo = stationdriver.get_caltableinfo(rcumode)
            caltabinfos.append(caltabinfo)

    if pointing is not None:
        # Real beam start:
        print("Now running real beam... @ {}".format(
                                                datetime.datetime.utcnow()))
        dir_bmctl = ilisa.observations.directions.normalizebeamctldir(pointing)
        rcu_setup_cmd, beamctl_cmds = stationdriver.streambeams(freqbndobj,
                                                                dir_bmctl)
        beamstarted = datetime.datetime.utcnow()
        rectime = beamstarted
        scan_id = stationdriver.get_scanid(beamstarted)
    else:
        rcu_setup_cmd, beamctl_cmds = "", ""
        rectime = starttime
        scan_id = stationdriver.get_scanid()

    if rec_bfs:
        scanresult['rec'].append('bfs')
        scanresult['bfs'] = dataIO.ScanRecInfo()
        scanresult['bfs'].set_stnid(stnid)
        scanpath_bfdat = os.path.join(stationdriver.bf_data_dir, scan_id)
        lanes = tuple(freqbndobj.getlanes().keys())
        datafiles, logfiles = stationdriver.dru_interface.rec_bf_proxy(rectime,
            duration_tot, lanes, band, scanpath_bfdat, stationdriver.bf_port0,
                                                                       stnid)
        bfsnametime = rectime.strftime("%Y%m%d_%H%M%S")
        this_obsinfo = dataIO.LDatInfo('bfs', bfsnametime, beamctl_cmds,
                                       rcu_setup_cmd, caltabinfos)
        scanresult['bfs'].add_obs(this_obsinfo)
        bfsdatapaths = []
        bfslogpaths = []
        for lane in lanes:
            datafileguess = datafiles.pop()
            dumplogname = logfiles.pop()
            if not datafileguess:
                _outdumpdir, _outarg, datafileguess, dumplogname = \
                    bfbackend.bfsfilepaths(lane, rectime, band, scanpath_bfdat,
                                           stationdriver.bf_port0,
                                           stationdriver.get_stnid())
            bfsdatapaths.append(datafileguess)
            bfslogpaths.append(dumplogname)
    else:
        scanpath_bfdat = None
        print("Not recording bfs")
    sys.stdout.flush()

    if bsx_type is not None:
        scanresult['rec'].append('bsx')
        scanresult['bsx'] = dataIO.ScanRecInfo()
        scanresult['bsx'].set_stnid(stnid)
        # Record statistic for duration_tot seconds
        if bsx_type == 'bst' or bsx_type == 'sst':
            rspctl_cmd = stationdriver.rec_bsx(bsx_type, integration,
                                               duration_tot)
            obsdatetime_stamp = stationdriver.get_data_timestamp(-1)
            curr_obsinfo = dataIO.LDatInfo(bsx_type, obsdatetime_stamp,
                                           beamctl_cmds, rspctl_cmd,
                                           caltabinfos)
            scanresult['bsx'].add_obs(curr_obsinfo)
        elif bsx_type == 'xst':
            nrsubbands = freqbndobj.nrsubbands()
            if duration_frq is None:
                if nrsubbands > 1:
                    duration_frq = integration
                else:
                    duration_frq = duration_tot
            # TODO Consider that specified duration is not the same as actual
            #  duration. Each step in frequency sweep take about 6s for 1s int.
            (rep, _rst) = divmod(duration_tot, duration_frq * nrsubbands)
            rep = int(rep)
            if rep == 0:
                warnings.warn("""Total duration too short for 1 full repetition.
Will increase total duration to get 1 full repetition.""")
                duration_tot = duration_frq * nrsubbands
                rep = 1
            # Repeat rep times (freq sweep)
            for _itr in range(rep):
                # Start freq sweep
                for sb_rcumode in freqbndobj.sb_range:
                    if ':' in sb_rcumode:
                        sblo, sbhi = sb_rcumode.split(':')
                        subbands = range(int(sblo), int(sbhi) + 1)
                    else:
                        subbands = [int(sb) for sb in sb_rcumode.split(',')]
                    for subband in subbands:
                        # Record data
                        rspctl_cmd = stationdriver.rec_bsx(bsx_type,
                                                           integration,
                                                           duration_frq,
                                                           subband)
                        obsdatetime_stamp = stationdriver.get_data_timestamp(-1)
                        curr_obsinfo = dataIO.LDatInfo('xst', obsdatetime_stamp,
                                                       beamctl_cmds, rspctl_cmd,
                                                       caltabinfos=caltabinfo,
                                                       septonconf=septonconf)
                        scanresult['bsx'].add_obs(curr_obsinfo)
        else:
            raise Exception('LOFAR statistic datatype "{}" unknown.\
                            (Known are bst, sst, xst)'.format(bsx_type))

    if not rec_acc and not rec_bfs and bsx_type is None:
        print("Will run beam with no active recording for {} seconds.".format(
                                                                duration_tot))
        # Since we're not recording anything, just do nothing for the
        # duration_tot.
        time.sleep(duration_tot)

    # Finished recording
    stationdriver.stop_beam()

    if todo_tof:
        stationdriver.stop_tof()

    # Work out where station-correlated data should be stored:
    scanpath_scdat = os.path.join(stationdriver.scanpath, scan_id)
    # and create the directory: (may not have ldat if no rec but will have info
    # files)
    os.makedirs(scanpath_scdat)

    if rec_acc:
        # Switch back to normal state i.e. turn-off ACC dumping:
        stationdriver.acc_mode(enable=False)

        # Create obsinfo each ACC file
        _, acc_files =  stationdriver.getdatalist()
        for acc_file in acc_files:
            obsid, _ = acc_file.split('_acc_')
            this_obsinfo =  dataIO.LDatInfo('acc', obsid, beamctl_cmds,
                                            rcu_setup_cmd)
            scanresult['acc'].add_obs(this_obsinfo)

        # Set scanrecinfo
        acc_integration = 1.0
        scanresult['acc'].set_scanrecparms('acc', band, duration_tot, pointing,
                                                  acc_integration, allsky)
        scanresult['acc'].set_scanpath(scanpath_scdat)
        scanrecpath = scanresult['acc'].get_scanrecpath()

        # Transfer data from LCU to DAU
        if os.path.exists(scanrecpath):
            print("Appropriate directory exists already (will put data here)")
        else:
            print("Creating directory "+scanrecpath+" for ACC "
                  + str(duration_tot)+" s rcumode="+str(rcumode)+" calibration")
            os.makedirs(scanrecpath)

        # Move ACC dumps to storage
        accsrcfiles = stationdriver.get_ACCsrcDir() + "/*.dat"
        stationdriver.movefromlcu(accsrcfiles, scanrecpath)

    if bsx_type is not None:
        # Set scanrecinfo
        scanresult['bsx'].set_scanrecparms(bsx_type, freqbndobj.arg,
                                           duration_tot, pointing,
                                           integration, allsky=allsky)
        # Move data to archive
        scanresult['bsx'].set_scanpath(scanpath_scdat)
        scanrecpath = scanresult['bsx'].get_scanrecpath()
        stationdriver.movefromlcu(stationdriver.get_lcuDumpDir() + "/*.dat",
                                  scanrecpath)

    if rec_bfs:
        scanresult['bfs'].set_stnid(stationdriver.get_stnid())
        scanresult['bfs'].set_scanrecparms('bfs', band, duration_tot, pointing,
                                           allsky=allsky)
        # Make a project folder for BFS data
        scanresult['bfs'].set_scanpath(scanpath_scdat)
        scanrecpath = scanresult['bfs'].get_scanrecpath()
        # Create BFS destination folder on DPU:
        os.makedirs(scanrecpath)
        if stationdriver.dru_interface.hostname == 'localhost':
            # Make soft links to actual BFS files and move logs to scanrec folder
            for lane in lanes:
                if bfsdatapaths[lane] is not None:
                    os.symlink(bfsdatapaths[lane], os.path.join(scanrecpath,
                                            os.path.basename(bfsdatapaths[lane])))
                if bfslogpaths[lane] is not None:
                    shutil.move(bfslogpaths[lane], scanrecpath)

    stationdriver.cleanup()
    # Necessary due to possible forking
    stationdriver.halt_observingstate_when_finished = shutdown
    scanresult['scan_id'] = scan_id
    scanresult['scanpath_scdat'] = scanpath_scdat
    return scanresult
