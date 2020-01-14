import sys
import os
import time
import datetime
import shutil
import warnings
import inspect

import ilisa.observations.directions
import ilisa.observations.modeparms as modeparms
import ilisa.observations.beamformedstreams.bfbackend as bfbackend
import ilisa.observations.dataIO as dataIO


class ObsPrograms(object):

    def __init__(self, stationdriver):
        self.stationdriver = stationdriver
        self.lcu_interface = stationdriver.lcu_interface

    def getprogram(self, programname):
        programpointer = getattr(self, programname)
        defargstart = None
        # # FIXME Only return args without defaults
        # defarg = inspect.getargspec(programpointer).defaults
        # if defarg is not None:
        #     defargstart = -len(defarg)
        programargs = inspect.getargspec(programpointer).args[1:defargstart]
        return programpointer, programargs

    def _streambeams_mltfreq(self, freqbndobj, pointing, recDuration=float('inf'),
                             attenuation=0, DUMMYWARMUP=False):
        """Form beams with station."""
        bits = freqbndobj.bits
        if DUMMYWARMUP:
            print("Warning warnup not currently implemented")
        beamctl_cmds = []
        for bandbeamidx in range(len(freqbndobj.rcumodes)):
            antset = freqbndobj.antsets[bandbeamidx]
            rcumode = freqbndobj.rcumodes[bandbeamidx]
            beamletIDs = freqbndobj.beamlets[bandbeamidx]
            subbands =  freqbndobj.sb_range[bandbeamidx]
            rcusel = freqbndobj.rcusel[bandbeamidx]
            beamctl_main = self.lcu_interface.run_beamctl(beamletIDs, subbands,
                                                              rcumode, pointing, rcusel)
            beamctl_cmds.append(beamctl_main)
        rcu_setup_cmd = self.lcu_interface.rcusetup(bits, attenuation)
        return rcu_setup_cmd, beamctl_cmds

    def do_bfs_OW(self, freqbndobj, duration_tot, pointing, bfdsesdumpdir,
                  starttime='NOW'):
        """Record BeamFormed Streams (BFS) with particular beamlet allocation."""

        band = freqbndobj.rcubands[0]
        ###
        bits = 8  # 8
        attenuation = None
        # Subbands allocation
        if band == '10_90' or band == '30_90':
            # LBA
            lanes = (0, 1)  # (0,1)
            beamletIDs = '0:243'  # '0:243'
            subbandNrs = '164:407'  # '164:407'
        elif band == '110_190':
            # HBAlo
            lanes = (0, 1, 2, 3)  # Normally (0,1,2,3) for all 4 lanes.
            beamletIDs = '0:487'
            subbandNrs = '12:499'
        elif band == '210_250':
            # HBAhi
            lanes = (0, 1)
            beamletIDs = '0:243'
            subbandNrs = '12:255'
        else:
            raise ValueError(
                "Wrong band: should be 10_90 (LBA), 110_190 (HBAlo) or 210_250 (HBAhi).")

        # Wait until it is time to start
        pause = 5  # Sufficient?
        if starttime != "NOW":
            rectime = starttime
        else:
            rectime = datetime.datetime.utcnow()
        self.stationdriver._waittoboot(rectime, pause)

        # Necessary since fork creates multiple instances of myobs and each one
        # will call it's __del__ on completion and __del__ shutdown...
        shutdown = self.stationdriver.halt_observingstate_when_finished
        self.stationdriver.halt_observingstate_when_finished = False
        self.stationdriver.exit_check = False

        dir_bmctl = ilisa.observations.directions.normalizebeamctldir(pointing)

        # BEGIN Dummy or hot beam start: (takes about 10sec)
        print("Running warmup beam... @ {}".format(datetime.datetime.utcnow()))
        self.lcu_interface.rcusetup(bits, attenuation)  # setting bits is necessary
        self.lcu_interface.run_beamctl(beamletIDs, subbandNrs, band, dir_bmctl)
        self.lcu_interface.stop_beam()
        # END Dummy or hot start

        print("Pause {}s after boot.".format(pause))
        time.sleep(pause)

        # Real beam start:
        print("Now running real beam... @ {}".format(datetime.datetime.utcnow()))
        rcu_setup_cmd = self.lcu_interface.rcusetup(bits, attenuation)
        beamctl_cmd = self.lcu_interface.run_beamctl(beamletIDs, subbandNrs, band,
                                                     dir_bmctl)
        beamstart = datetime.datetime.utcnow()
        timeleft = rectime - beamstart
        if timeleft.total_seconds() < 0.:
            rectime = beamstart
        print("(Beam started) Time left before recording: {}".format(
            timeleft.total_seconds()))
        bfsnametime = starttime.strftime("%Y%m%d_%H%M%S")
        obsinfo = dataIO.LDatInfo('bfs', bfsnametime, beamctl_cmd, rcu_setup_cmd,
                                  caltabinfos="")

        REC = True
        if REC == True:
            port0 = self.stationdriver.bf_port0
            stnid = self.lcu_interface.stnid
            bfbackend.rec_bf_streams(rectime, duration_tot, lanes, band, bfdsesdumpdir,
                                     port0, stnid)
            bfsdatapaths = []
            for lane in lanes:
                outdumpdir, outarg, datafileguess, dumplogname = \
                    bfbackend.bfsfilepaths(lane, rectime, band, bfdsesdumpdir, port0,
                                           stnid)
                bfsdatapaths.append(datafileguess)
        else:
            print("Not recording")
            time.sleep(duration_tot)
        sys.stdout.flush()
        self.lcu_interface.stop_beam()
        self.stationdriver.halt_observingstate_when_finished = shutdown
        return [obsinfo]


def record_obsprog(stationdriver, scan):
    """At starttime execute the observation program specified by the obsfun method
    pointer and run with arguments specified by obsargs dict.
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
    stationdriver.lcu_interface.selectCalTable(CALTABLESRC)

    # Prepare for obs program.
    try:
        stationdriver.goto_observingstate()
    except RuntimeError as e:
        raise RuntimeError(e)

    # Run the observation program:
    obsinfolist = obsfun(**obsargs)
    # Stop program beam
    stationdriver.lcu_interface.stop_beam()
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
            caltabinfo = stationdriver.lcu_interface.getCalTableInfo(rcumode)
            caltabinfos.append(caltabinfo)
        # Add obsinfos to scanrecs
        for obsinfo in obsinfolist:
            obsinfo.caltabinfos = caltabinfos
            scanrec.add_obs(obsinfo)
        # Move data to archive
        stationdriver.movefromlcu(stationdriver.lcu_interface.lcuDumpDir + "/*.dat",
                                  scanpath_scdat)
        scanrec.path = scanpath_scdat
        scanresult[datatype] = scanrec
    scanresult['scan_id'] = scan_id
    scanresult['scanpath_scdat'] = scanpath_scdat
    return scanresult


def record_scan(stationdriver, freqbndobj, duration_tot, pointing,
                starttime='NOW', rec=[], integration=1.0, allsky=False,
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
            Pointing source direction of scan. Can be a source name or a beamctl
            direction str.
        starttime: str
            The time at which scan should start.
        rec: [str]
            The types of LOFAR data to record.
            Is a list of up to three of the following data types:
                'acc', 'bfs', 'bst', 'sst', 'xst'
            of which 3 latter types are mutual exclusive.

            'acc': Record autocovariance-cubes (ACC) files. The actual total duration will
            at most be duration_tot. (Usually it will be shorter so it fits within the
            cadence of whole ACC aquisition, which is 512+7=519 seconds). ACC files are
            the covariance of all array elements with each as a function of subband.

            'bfs': Record BeamFormed Stream data.

            'bst': Record Beamlet STatistics data.

            'sst': Record Subband STatistics data.

            'xst': Record Xrosslet STatistics data.
        duration_frq: float
            Duration in seconds of a sampled frequency within a frequency sweep scan.
        allsky: bool, optional
            HBA allsky mode.
    """
    rec_acc = False; rec_bfs = False; bsx_type = None
    for recreq in rec:
        if recreq == 'acc':
            rec_acc = True
        elif recreq == 'bfs':
            rec_bfs = True
        else:
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

    scanresult = {}
    scanresult['acc'] = dataIO.ScanRecInfo()
    scanresult['bfs'] = dataIO.ScanRecInfo()
    scanresult['bsx'] = dataIO.ScanRecInfo()
    for k in scanresult.keys():
        scanresult[k].set_stnid(stationdriver.get_stnid())

    # Setup Calibration tables on LCU:
    CALTABLESRC = 'default'   # FIXME put this in args
    ## (Only BST uses calibration tables)
    # Choose between 'default' or 'local'
    stationdriver.lcu_interface.selectCalTable(CALTABLESRC)

    # Prepare for obs program.
    try:
        stationdriver.goto_observingstate()
    except RuntimeError as e:
        raise RuntimeError(e)

    duration_tot_req = duration_tot
    band = freqbndobj.rcubands[0]
    rcumode = freqbndobj.rcumodes[0]

    if todo_tof:
        septonconf = stationdriver._setup_tof()
    else:
        septonconf = None

    if rec_acc:
        # Also duration of ACC sweep since each sb is 1 second.
        dur1acc = modeparms.TotNrOfsb  # Duration of one ACC
        interv2accs = 7  # time between end of one ACC and start of next one
        acc_cadence = dur1acc+interv2accs  # =519s time between start of two ACCs
        (nraccs, timrest) = divmod(duration_tot_req, acc_cadence)
        if timrest > dur1acc:
            nraccs += 1
        duration_tot = nraccs*acc_cadence-interv2accs
        if duration_tot != duration_tot_req:
            print("""Note: will use total duration {}s to fit with ACC
                  cadence.""".format(duration_tot))

        # Make sure swlevel=<2
        stationdriver.lcu_interface.set_swlevel(2)

        # Set CalServ.conf to dump ACCs:
        stationdriver.lcu_interface.acc_mode(enable=True)

        # Boot to swlevel 3 so the calserver service starts
        stationdriver.lcu_interface.set_swlevel(3)

    # Wait until it is time to start
    now = datetime.datetime.utcnow()
    if starttime == "NOW":
        rectime = now
    else:
        rectime = starttime
    timeleft = rectime - now
    secondsleft = int(timeleft.total_seconds())
    if secondsleft < 0:
        secondsleft = 0
    print("Waiting {}s to start scan ".format(secondsleft))
    time.sleep(secondsleft)

    # Necessary since fork creates multiple instances of myobs and each one
    # will call it's __del__ on completion and __del__ shutdown...
    shutdown = stationdriver.halt_observingstate_when_finished
    stationdriver.halt_observingstate_when_finished = False
    stationdriver.exit_check = False

    caltabinfos = [""]
    # Get metadata about caltables to be used
    if not allsky:
        for rcumode in freqbndobj.rcumodes:
            caltabinfo = stationdriver.lcu_interface.getCalTableInfo(rcumode)
            caltabinfos.append(caltabinfo)

    if pointing is not None:
        # Real beam start:
        print("Now running real beam... @ {}".format(datetime.datetime.utcnow()))
        dir_bmctl = ilisa.observations.directions.normalizebeamctldir(pointing)
        rcu_setup_cmd, beamctl_cmds = stationdriver.streambeams(freqbndobj, dir_bmctl)
        beamstarted = datetime.datetime.utcnow()
        rectime = beamstarted
        scan_id = stationdriver.get_scanid(beamstarted)
    else:
        rcu_setup_cmd, beamctl_cmds = "", ""
        scan_id = stationdriver.get_scanid()

    if rec_bfs:
        scanpath_bfdat = os.path.join(stationdriver.bf_data_dir, scan_id)
        stnid = stationdriver.lcu_interface.stnid
        lanes = tuple(freqbndobj.getlanes().keys())
        bfbackend.rec_bf_streams(rectime,
                                 duration_tot, lanes, band,
                                 scanpath_bfdat, stationdriver.bf_port0, stnid)
        bfsnametime = rectime.strftime("%Y%m%d_%H%M%S")
        this_obsinfo = dataIO.LDatInfo('bfs', bfsnametime, beamctl_cmds, rcu_setup_cmd,
                                       caltabinfos)
        scanresult['bfs'].add_obs(this_obsinfo)
    else:
        scanpath_bfdat = None
        print("Not recording bfs")
    sys.stdout.flush()

    if bsx_type is not None:
        # Record statistic for duration_tot seconds
        if bsx_type == 'bst':
            rspctl_cmd = stationdriver.lcu_interface.rec_bst(integration, duration_tot)
            # beamlet statistics also generate empty *01[XY].dat so remove:
            stationdriver.lcu_interface.rm(stationdriver.lcu_interface.lcuDumpDir + "/*01[XY].dat")
            obsdatetime_stamp = stationdriver.get_data_timestamp(-1)
            curr_obsinfo = dataIO.LDatInfo('bst', obsdatetime_stamp, beamctl_cmds,
                                           rspctl_cmd, caltabinfos)
            scanresult['bsx'].add_obs(curr_obsinfo)
        elif bsx_type == 'sst':
            caltabinfo = ""
            rspctl_cmd = stationdriver.lcu_interface.rec_sst(integration, duration_tot)
            obsdatetime_stamp = stationdriver.get_data_timestamp(-1)
            curr_obsinfo = dataIO.LDatInfo('sst', obsdatetime_stamp, beamctl_cmds,
                                           rspctl_cmd, caltabinfo)
            scanresult['bsx'].add_obs(curr_obsinfo)
        elif bsx_type == 'xst':
            caltabinfo = ""  # No need for caltab info for xst data
            nrsubbands = freqbndobj.nrsubbands()
            if duration_frq is None:
                if nrsubbands > 1:
                    duration_frq = integration
                else:
                    duration_frq = duration_tot
            # TODO Consider that specified duration is not the same as actual
            #  duration. Each step in frequency sweep take about 6s for 1s int.
            (rep, rst) = divmod(duration_tot, duration_frq * nrsubbands)
            rep = int(rep)
            if rep == 0:
                warnings.warn("""Total duration too short for 1 full repetition.
Will increase total duration to get 1 full repetition.""")
                duration_tot = duration_frq * nrsubbands
                rep = 1
            # Repeat rep times (freq sweep)
            for itr in range(rep):
                # Start freq sweep
                for sb_rcumode in freqbndobj.sb_range:
                    if ':' in sb_rcumode:
                        sblo, sbhi = sb_rcumode.split(':')
                        subbands = range(int(sblo), int(sbhi) + 1)
                    else:
                        subbands = [int(sb) for sb in sb_rcumode.split(',')]
                    for subband in subbands:
                        # Record data
                        rspctl_cmd = stationdriver.lcu_interface.rec_xst(subband,
                                                                         integration,
                                                                         duration_frq)
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
        # Since we're not recording anything, just do nothing for the duration_tot.
        time.sleep(duration_tot)

    # Finished recording
    stationdriver.lcu_interface.stop_beam()

    if todo_tof:
        stationdriver.lcu_interface.set_swlevel(3)

    # Work out where station-correlated data should be stored:
    scanpath_scdat = os.path.join(stationdriver.scanpath, scan_id)
    # and create the directory: (may not have ldat if no rec but will have info files)
    os.makedirs(scanpath_scdat)

    if rec_acc:
        # Switch back to normal state i.e. turn-off ACC dumping:
        stationdriver.lcu_interface.set_swlevel(2)
        stationdriver.lcu_interface.acc_mode(enable=False)
        stationdriver.lcu_interface.set_swlevel(3)

        # Transfer data from LCU to DAU
        obsdatetime_stamp = stationdriver.get_data_timestamp(ACC=True)
        accsrcfiles = stationdriver.lcu_interface.ACCsrcDir + "/*.dat"
        scanrecpath = \
            os.path.join(scanpath_scdat,
                         '{}_{}_rcu{}_dur{}'.format(stationdriver.lcu_interface.stnid,
                                                    obsdatetime_stamp, rcumode,
                                                    duration_tot))
        if int(rcumode) > 4:  # rcumodes more than 4 need pointing
            scanrecpath += "_" + pointing
        scanrecpath += "_acc"
        if os.path.exists(scanrecpath):
            print("Appropriate directory exists already (will put data here)")
        else:
            print("Creating directory "+scanrecpath+" for ACC "+str(duration_tot)
                  + " s rcumode="+str(rcumode)+" calibration")
            os.makedirs(scanrecpath)

        # Move ACC dumps to storage
        stationdriver.movefromlcu(accsrcfiles, scanrecpath)
        accdestfiles = os.listdir(scanrecpath)

        # Create obsinfo each ACC file
        for destfile in accdestfiles:
            obsid, _ = destfile.split('_acc_')
            this_obsinfo =  dataIO.LDatInfo('acc', obsid, beamctl_cmds, rcu_setup_cmd)
            scanresult['acc'].add_obs(this_obsinfo)

        # Set scanrecinfo
        acc_integration = 1.0
        scanresult['acc'].set_stnid(stationdriver.get_stnid())
        scanresult['acc'].set_scanrecparms('acc', band, duration_tot, pointing,
                                                  acc_integration, allsky)
        scanresult['acc'].path = scanrecpath

    if bsx_type is not None:
        # Move data to archive
        scanrecfolder = scanresult['bsx'].scanrecfolder()
        scanrecpath = os.path.join(scanpath_scdat, scanrecfolder)
        stationdriver.movefromlcu(stationdriver.lcu_interface.lcuDumpDir + "/*.dat",
                                  scanrecpath)

        # Set scanrecinfo
        scanresult['bsx'].set_stnid(stationdriver.get_stnid())
        scanresult['bsx'].set_scanrecparms(bsx_type, freqbndobj.arg,
                                                  duration_tot, pointing,
                                                  integration, allsky=allsky)
        scanresult['bsx'].path = scanrecpath

    if rec_bfs:
        # Make a project folder for BFS data
        scanrecfolder = scanresult['bfs'].scanrecfolder()
        scanrecpath = os.path.join(scanpath_scdat, scanrecfolder)
        print("Creating BFS destination folder on DPU:\n{}".format(scanrecpath))
        os.makedirs(scanrecpath)
        bfsdatapaths = []
        bfslogpaths = []
        for lane in lanes:
            outdumpdir, outarg, datafileguess, dumplogname = \
                bfbackend.bfsfilepaths(lane, rectime,
                                       band, scanpath_bfdat,
                                       stationdriver.bf_port0, stationdriver.get_stnid())
            bfsdatapaths.append(datafileguess)
            bfslogpaths.append(dumplogname)
        # Make soft links to actual BFS files and move logs to scanrec folder
        for lane in lanes:
            os.symlink(bfsdatapaths[lane], os.path.join(scanrecpath,
                                                        os.path.basename(
                                                            bfsdatapaths[lane])))
            shutil.move(bfslogpaths[lane], scanrecpath)
        scanresult['bfs'].set_stnid(stationdriver.get_stnid())
        scanresult['bfs'].set_scanrecparms('bfs', band, duration_tot, pointing,
                                                  allsky=allsky)
        scanresult['bfs'].path = scanrecpath

    stationdriver.lcu_interface.cleanup()
    # Necessary due to possible forking
    stationdriver.halt_observingstate_when_finished = shutdown
    scanresult['scan_id'] = scan_id
    scanresult['scanpath_scdat'] = scanpath_scdat
    return scanresult
