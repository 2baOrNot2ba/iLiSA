import sys
import os
import time
import datetime
import shutil
import warnings
import copy
import inspect

import ilisa.observations.directions
import ilisa.observations.modeparms as modeparms
import ilisa.observations.beamformedstreams.bfbackend as bfbackend
import ilisa.observations.dataIO as dataIO


class BasicObsPrograms(object):

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

    def do_bfs_OW(self, freqbndobj, duration_tot, pointsrc, bfdsesdumpdir,
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
        pointing = ilisa.observations.directions.normalizebeamctldir(pointsrc)

        # Wait until it is time to start
        pause = 5  # Sufficient?
        if starttime != "NOW":
            starttimestr = starttime.strftime("%Y-%m-%dT%H:%M:%S")
        else:
            starttimestr = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
        st = self.stationdriver._waittoboot(starttimestr, pause)

        # From swlevel 0 it takes about 1:30min? to reach swlevel 3
        print("Booting @ {}".format(datetime.datetime.utcnow()))

        # Necessary since fork creates multiple instances of myobs and each one
        # will call it's __del__ on completion and __del__ shutdown...
        shutdown = self.stationdriver.halt_observingstate_when_finished
        self.stationdriver.halt_observingstate_when_finished = False
        self.stationdriver.exit_check = False

        # BEGIN Dummy or hot beam start: (takes about 10sec)
        # TODO: This seems necessary, otherwise beamctl will not start up next time,
        #       although it should not have to necessary.)
        print("Running warmup beam... @ {}".format(datetime.datetime.utcnow()))
        self.lcu_interface.run_beamctl(beamletIDs, subbandNrs, band, pointing)
        self.lcu_interface.rcusetup(bits,
                                    attenuation)  # setting bits also seems necessary
        self.lcu_interface.stop_beam()
        # END Dummy or hot start

        print("Pause {}s after boot.".format(pause))
        time.sleep(pause)

        # Real beam start:
        print("Now running real beam... @ {}".format(datetime.datetime.utcnow()))
        beamctl_cmds = self.lcu_interface.run_beamctl(beamletIDs, subbandNrs, band,
                                                     pointing)
        rcu_setup_cmd = self.lcu_interface.rcusetup(bits, attenuation)
        beamstart = datetime.datetime.utcnow()
        timeleft = st - beamstart
        if timeleft.total_seconds() < 0.:
            starttimestr = beamstart.strftime("%Y-%m-%dT%H:%M:%S")
        print("(Beam started) Time left before recording: {}".format(
            timeleft.total_seconds()))

        REC = True
        if REC == True:
            port0 = self.stationdriver.bf_port0
            stnid = self.lcu_interface.stnid
            bfbackend.rec_bf_streams(starttimestr, duration_tot, lanes, band,
                                     bfdsesdumpdir, port0, stnid)
        else:
            print("Not recording")
            time.sleep(duration_tot)
        sys.stdout.flush()
        self.lcu_interface.stop_beam()
        self.stationdriver.halt_observingstate_when_finished = shutdown
        return None


def record_scan(stationdriver, freqbndobj, integration, duration_tot, pointing, pointsrc,
                starttime='NOW', rec_stat_type=None, rec_bfs=False, duration_frq=None,
                do_acc=False, allsky=False, scanmeta=None):
    """Run a generic scan.

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
        pointing:
            Pointing direction of scan.
        pointsrc: str
            The field name.
        starttime: str
            The time at which scan should start.
        rec_stat_type: str
            The type of LOFAR statistic to be recorded. Can be either: 'bst', 'sst',
            or 'xst'.
        rec_bfs: bool
            Record the beam-formed stream.
        duration_frq: float
            Duration in seconds of a sampled frequency within a frequency sweep scan.
        do_acc: bool
            Perform calibration observation mode on station. Also known as ACC
            mode. The actual total duration will at most be duration_tot_req.
            (Usually it will be shorter so it fits within the cadence of whole ACC
            aquisition, which is 512+7=519 seconds).

            Some technical details: swlevel needs to cycle down to 2 (or less) and then to 3.
            If swlevel is kept at 3 (i.e. exit_obsstate=False), then ACC will continue to be
            produced, until swlevel goes below 2.

            ACC files are autocovariance-cubes: the covariance of all array
            elements with each as a function of subband. These files are generated
            by the MAC service called CalServer. It run at swlevel 3 and is
            configured in the file lofar/etc/CalServer.conf. Note subband
            integration is always 1s, so ACC file is dumped after 512 seconds.
        allsky: bool, optional
            HBA allsky mode.
        scanmeta: ScanMeta, optional
            The metadata for the scan. Contains: the directory where the scan should
            be stored, the directory where the bfs data should be dumped, and
            a dict of ScanRec for potential scan recordings. The ScanRec dict will
            get updated with the actual ScanRec settings used.
    """
    req_allsky = allsky; del allsky
    arraytype = freqbndobj.antsets[-1][:3]
    # Mode logic
    todo_tof = False  # This is the case when arraytype=='LBA'
    if req_allsky and arraytype == 'HBA':
        if pointing is not None:
            raise ValueError('hba-allsky and beam cannot run simultaneously.')
        if rec_stat_type == 'bst':
            raise ValueError('bst and hba-allsky cannot be combined.')
        if do_acc:
            raise ValueError('acc and hba-allsky cannot be combined.')
        if rec_bfs:
            raise ValueError('bfs and hba-allsky cannot be combined.')
        todo_tof = True
    if not req_allsky and pointing is None:
        pointsrc = 'Z'

    pointing = ilisa.observations.directions.normalizebeamctldir(pointsrc)

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

    if do_acc:
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
        #sst_integration = int(acc_cadence)

        # Make sure swlevel=<2
        stationdriver.lcu_interface.set_swlevel(2)

        # Set CalServ.conf to dump ACCs:
        stationdriver.lcu_interface.acc_mode(enable=True)

        # Boot to swlevel 3 so the calserver service starts
        stationdriver.set_swlevel(3)

    # Wait until it is time to start
    pause = 5  # Sufficient?
    if starttime != "NOW":
        starttimestr = starttime.strftime("%Y-%m-%dT%H:%M:%S")
    else:
        starttimestr = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
    st = stationdriver._waittoboot(starttimestr, pause)

    # Necessary since fork creates multiple instances of myobs and each one
    # will call it's __del__ on completion and __del__ shutdown...
    shutdown = stationdriver.halt_observingstate_when_finished
    stationdriver.halt_observingstate_when_finished = False
    stationdriver.exit_check = False

    caltabinfos = []
    if pointing is not None:
        # Get metadata about caltables to be used
        for rcumode in freqbndobj.rcumodes:
            caltabinfo = stationdriver.lcu_interface.getCalTableInfo(rcumode)
            caltabinfos.append(caltabinfo)
        # warmup used to be here
        print("Pause {}s after boot.".format(pause))
        time.sleep(pause)

        # Real beam start:
        print("Now running real beam... @ {}".format(datetime.datetime.utcnow()))
        rcu_setup_cmd, beamctl_cmds = stationdriver.streambeams(freqbndobj, pointing)
        beamstarted = datetime.datetime.utcnow()
        timeleft = st - beamstarted
        if timeleft.total_seconds() < 0.:
            starttimestr = beamstarted.strftime("%Y-%m-%dT%H:%M:%S")
        print("(Beam started) Time left before recording: {}".format(
            timeleft.total_seconds()))
    else:
        print("No pointing...")
        rcu_setup_cmd = ""
        beamctl_cmds = ""
        beamstarted = None

    if rec_bfs:
        bfdsesdumpdir = scanmeta.bfdsesdumpdir
        stnid = stationdriver.lcu_interface.stnid
        lanes = tuple(freqbndobj.getlanes().keys())
        bfbackend.rec_bf_streams(starttimestr, duration_tot, lanes, band,
                                 bfdsesdumpdir, stationdriver.bf_port0, stnid)
        bfsdatapaths = []
        bfslogpaths = []
        for lane in lanes:
            outdumpdir, outarg, datafileguess, dumplogname =\
                bfbackend.bfsfilepaths(lane, starttimestr, band, bfdsesdumpdir,
                                       stationdriver.bf_port0, stnid)
            bfsdatapaths.append(datafileguess)
            bfslogpaths.append(dumplogname)
    else:
        print("Not recording bfs")
    sys.stdout.flush()

    if rec_stat_type is not None:
        # Record statistic for duration_tot seconds
        obsinfolist = []
        if rec_stat_type == 'bst':
            rspctl_cmd = stationdriver.lcu_interface.rec_bst(integration, duration_tot)
            # beamlet statistics also generate empty *01[XY].dat so remove:
            stationdriver.lcu_interface.rm(stationdriver.lcu_interface.lcuDumpDir + "/*01[XY].dat")
            obsdatetime_stamp = stationdriver.get_data_timestamp(-1)
            curr_obsinfo = dataIO.ObsInfo()
            curr_obsinfo.setobsinfo_fromparams('bst', obsdatetime_stamp, beamctl_cmds,
                                               rspctl_cmd, caltabinfos)
            obsinfolist.append(curr_obsinfo)
        elif rec_stat_type == 'sst':
            caltabinfo = ""
            rspctl_cmd = stationdriver.lcu_interface.rec_sst(integration, duration_tot)
            obsdatetime_stamp = stationdriver.get_data_timestamp(-1)
            curr_obsinfo = dataIO.ObsInfo()
            curr_obsinfo.setobsinfo_fromparams('sst', obsdatetime_stamp, beamctl_cmds,
                                               rspctl_cmd, caltabinfo)
            obsinfolist.append(curr_obsinfo)
        elif rec_stat_type == 'xst':
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
                        rspctl_cmd = stationdriver.lcu_interface.rec_xst(subband, integration,
                                                                duration_frq)
                        obsdatetime_stamp = stationdriver.get_data_timestamp(-1)
                        curr_obsinfo = dataIO.ObsInfo()
                        curr_obsinfo.setobsinfo_fromparams('xst', obsdatetime_stamp,
                                                           beamctl_cmds, rspctl_cmd,
                                                           caltabinfos=caltabinfo,
                                                           septonconf=septonconf)
                        obsinfolist.append(curr_obsinfo)
        else:
            raise Exception('LOFAR statistic datatype "{}" unknown.\
                            (Known are bst, sst, xst)'.format(rec_stat_type))
    else:
        obsinfolist = None

    if not do_acc and not rec_bfs and rec_stat_type is None:
        # Since we're not recording anything, just do nothing for the duration_tot.
        time.sleep(duration_tot)

    # Finished recording
    stationdriver.lcu_interface.stop_beam()

    if todo_tof:
        stationdriver.lcu_interface.set_swlevel(3)

    # Work out where data should be stored:
    scanpath = stationdriver.get_scanpath(scanmeta.sesspath, beamstarted)

    if do_acc:
        # Switch back to normal state i.e. turn-off ACC dumping:
        stationdriver.lcu_interface.set_swlevel(2)
        stationdriver.lcu_interface.acc_mode(enable=False)
        stationdriver.lcu_interface.set_swlevel(3)

        # Transfer data from LCU to DAU
        obsdatetime_stamp = stationdriver.get_data_timestamp(ACC=True)
        accsrcfiles = stationdriver.lcu_interface.ACCsrcDir + "/*.dat"
        acc_destfolder = \
            os.path.join(scanpath,
                         '{}_{}_rcu{}_dur{}'.format(stationdriver.lcu_interface.stnid,
                                                    obsdatetime_stamp, rcumode,
                                                    duration_tot))
        if int(rcumode) > 4:  # rcumodes more than 4 need pointing
            acc_destfolder += "_"+pointsrc
        acc_destfolder += "_acc"
        if os.path.exists(acc_destfolder):
            print("Appropriate directory exists already (will put data here)")
        else:
            print("Creating directory "+acc_destfolder+" for ACC "+str(duration_tot)
                  + " s rcumode="+str(rcumode)+" calibration")
            os.makedirs(acc_destfolder)

        # Move ACC dumps to storage
        stationdriver.movefromlcu(accsrcfiles, acc_destfolder)
        accdestfiles = os.listdir(acc_destfolder)

        # - Create project header
        acc_integration = 1.0
        scanmeta.scanrecs['acc'].set_stnid(stationdriver.get_stnid())
        scanmeta.scanrecs['acc'].set_obsfolderinfo('acc', obsdatetime_stamp, band,
                                          acc_integration, duration_tot, pointing)
        scanmeta.scanrecs['acc'].write_scan_rec(acc_destfolder)

        # - Create header for each ACC file
        for destfile in accdestfiles:
            filedatestr, filetimestr, _ = destfile.split('_', 2)
            filedtstr = '_'.join([filedatestr, filetimestr])
            scanmeta.scanrecs['acc'].new_obsinfo()
            scanmeta.scanrecs['acc'].obsinfos[-1].setobsinfo_fromparams('acc',
                filedtstr, beamctl_cmds, '')
            scanmeta.scanrecs['acc'].obsinfos[-1].create_LOFARst_header(
                acc_destfolder)
        scanmeta.scanrecs['acc'].datapath = acc_destfolder

    if rec_stat_type is not None and obsinfolist is not None:
        obsinfo = obsinfolist[0]

        # Move data to archive
        bsxSTobsEpoch, datapath = obsinfo.getobsdatapath(scanpath)
        stationdriver.movefromlcu(stationdriver.lcu_interface.lcuDumpDir + "/*.dat", datapath)
        for curr_obsinfo in obsinfolist:
            curr_obsinfo.create_LOFARst_header(datapath)
        # Prepare metadata for session on this station
        scanmeta.scanrecs['bsx'].set_stnid(stationdriver.get_stnid())
        scanmeta.scanrecs['bsx'].set_obsfolderinfo(obsinfo.LOFARdatTYPE, bsxSTobsEpoch,
                                     freqbndobj.arg, obsinfo.integration,
                                     obsinfo.duration_scan, obsinfo.pointing)
        scanmeta.scanrecs['bsx'].datapath = datapath
        scanmeta.scanrecs['bsx'].write_scan_rec(datapath)

    if rec_bfs:
        # Make a project folder for BFS data
        headertime = datetime.datetime.strptime(starttimestr, "%Y-%m-%dT%H:%M:%S"
                                                ).strftime("%Y%m%d_%H%M%S")
        scanmeta.scanrecs['bfs'].new_obsinfo()
        scanmeta.scanrecs['bfs'].obsinfos[-1].setobsinfo_fromparams('bfs', headertime,
                                            beamctl_cmds, rcu_setup_cmd, caltabinfos)
        bsxSTobsEpoch, datapath = \
            scanmeta.scanrecs['bfs'].obsinfos[-1].getobsdatapath(scanpath)
        print("Creating BFS destination folder on DPU:\n{}".format(datapath))
        os.makedirs(datapath)
        bfsdatapaths = []
        bfslogpaths = []
        for lane in lanes:
            outdumpdir, outarg, datafileguess, dumplogname = \
                bfbackend.bfsfilepaths(lane, starttimestr, band, bfdsesdumpdir,
                                       stationdriver.bf_port0, stationdriver.get_stnid())
            bfsdatapaths.append(datafileguess)
            bfslogpaths.append(dumplogname)
        # Make soft links to actual BFS files
        for lane in lanes:
            os.symlink(bfsdatapaths[lane], os.path.join(datapath,
                                                        os.path.basename(
                                                            bfsdatapaths[lane])))
            shutil.move(bfslogpaths[lane], datapath)
        # TODO: Move dump log files
        scanmeta.scanrecs['bfs'].set_stnid(stationdriver.get_stnid())
        scanmeta.scanrecs['bfs'].obsinfos[-1].create_LOFARst_header(datapath)
        integration = None
        scanmeta.scanrecs['bfs'].set_obsfolderinfo('bfs', headertime, band,
                                                   integration, duration_tot,
                                                   pointing)
        scanmeta.scanrecs['bfs'].datapath = datapath
        scanmeta.scanrecs['bfs'].write_scan_rec(datapath)

    stationdriver.lcu_interface.cleanup()
    # Necessary due to possible forking
    stationdriver.halt_observingstate_when_finished = shutdown
    return scanpath
