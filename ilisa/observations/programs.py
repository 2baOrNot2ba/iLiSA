import sys
import os
import time
import datetime
import copy
import inspect
import ilisa.observations.modeparms as modeparms
import ilisa.observations.dataIO as dataIO
import ilisa.observations.beamformedstreams.bfbackend as bfbackend

class BasicObsPrograms(object):

    def __init__(self, stationdriver):
        self.stationdriver = stationdriver
        self.lcu_interface = stationdriver.lcu_interface

    def getprogram(self, programname):
        programpointer = getattr(self, programname)
        programargs = inspect.getargspec(programpointer).args[1:]
        return programpointer, programargs

    def _fix_obsargs(self, kwargs_in):
        """Check and transform observational arguments to useful arguments."""
        kwargs_out = {}
        try:
            kwargs_out['pointing'] = modeparms.stdPointings(kwargs_in['pointsrc'])

        except KeyError:
            try:
                phi, theta, ref = kwargs_in['pointsrc'].split(',', 3)
                # FIXME Not always going to be correct
            except ValueError:
                raise ValueError("Error: %s invalid pointing syntax"\
                                 .format(kwargs_in['pointsrc']))
            else:
                kwargs_out['pointing'] = kwargs_in['pointsrc']
        if kwargs_in['integration'] > kwargs_in['duration_tot']:
            raise (ValueError, "integration {} is longer than duration_scan {}."
                               .format(kwargs_in['integration'],
                                       kwargs_in['duration_tot']))
        else:
            kwargs_out['integration'] = kwargs_in['integration']
            kwargs_out['duration_tot'] = kwargs_in['duration_tot']
        kwargs_out['freqbndobj'] = modeparms.FrequencyBand(kwargs_in['freqbndarg'])

        return kwargs_out

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

    def do_bstnew(self, freqbndobj, integration, duration_tot, pointing):
        """Do a Beamlet STatistic (bst) recording on station. frequency in Hz.
        """
        duration_scan = duration_tot
        obsinfolist = []

        # Start beamforming
        (rcu_setup_cmd, beamctl_cmds) = self._streambeams_mltfreq(freqbndobj, pointing)

        # Record data
        rspctl_cmd = self.lcu_interface.rec_bst(integration, duration_scan)

        # beamlet statistics also generate empty *01[XY].dat so remove:
        self.lcu_interface.rm(self.lcu_interface.lcuDumpDir+"/*01[XY].dat")

        # Collect some obsinfo (calinfo will be added later):
        obsdatetime_stamp = self.stationdriver.get_data_timestamp()
        curr_obsinfo = dataIO.ObsInfo()
        curr_obsinfo.setobsinfo_fromparams('bst', obsdatetime_stamp, beamctl_cmds,
                                           rspctl_cmd)
        obsinfolist.append(curr_obsinfo)
        return obsinfolist

    def do_sstnew(self, freqbndobj, integration, duration_tot, pointing):
        """Run an sst static."""

        duration_scan =  duration_tot
        obsinfolist = []

        # Use this to do a system temperature measurement.
        SYS_TEMP_MEAS = False

        # Start beamforming
        rcu_setup_CMD, beamctl_cmds = \
            self._streambeams_mltfreq(freqbndobj, pointing)
        if SYS_TEMP_MEAS:
            lbbalnaoff_CMD = self.lcu_interface.turnoffLBA_LNAs()
            #beamctl_CMD += "\n"+lbbalnaoff_CMD

        # Get some metadata about operational settings:
        caltabinfo = ""    # No need for caltab info

        # Record data
        rspctl_cmd = self.lcu_interface.rec_sst(integration, duration_scan)
        obsdatetime_stamp = self.stationdriver.get_data_timestamp(-1)
        curr_obsinfo = dataIO.ObsInfo()
        curr_obsinfo.setobsinfo_fromparams('sst', obsdatetime_stamp,
                                           beamctl_cmds, rspctl_cmd,
                                           caltabinfo)
        obsinfolist.append(curr_obsinfo)
        return obsinfolist

    def do_xstnew(self, freqbndobj, integration, duration_tot, pointing,
                  duration_scan=None):
        """New: Run an xst statistic towards the given pointing. This corresponds to
        a crosscorrelation of all elements at the given frequency for
        integration seconds over a duration_scan of seconds."""

        caltabinfo = ""  # No need for caltab info for xst data
        obsinfolist = []
        nrsubbands = freqbndobj.nrsubbands()
        if duration_scan is None:
            if nrsubbands > 1:
                duration_scan = integration
            else:
                duration_scan = duration_tot
        # FIXME When duration_tot is too small for 1 rep this will fail badly.
        (rep, rst) = divmod(duration_tot, duration_scan*nrsubbands)
        rep = int(rep)
        # Start beamforming
        rcu_setup_CMD, beamctl_cmds = self._streambeams_mltfreq(freqbndobj, pointing)
        # Repeat rep times (the freq sweep)
        for itr in range(rep):
            # Start freq sweep
            for sb_rcumode in freqbndobj.sb_range:
                if ':' in sb_rcumode:
                    sblo, sbhi = sb_rcumode.split(':')
                    subbands = range(int(sblo),int(sbhi)+1)
                else:
                    subbands = [int(sb) for sb in sb_rcumode.split(',')]
                for subband in subbands:
                    # Record data
                    rspctl_cmd = self.lcu_interface.rec_xst(subband, integration,
                                                            duration_scan)
                    obsdatetime_stamp = self.stationdriver.get_data_timestamp(-1)
                    curr_obsinfo = dataIO.ObsInfo()
                    curr_obsinfo.setobsinfo_fromparams('xst', obsdatetime_stamp,
                                                       beamctl_cmds, rspctl_cmd,
                                                       caltabinfo)
                    obsinfolist.append(curr_obsinfo)
        return obsinfolist

    def do_accnew(self, freqbndobj, duration_tot, pointing, pointsrc):
        """Perform calibration observation mode on station. Also known as ACC
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

        Parameters
        ----------
            freqbndobj: FrequencyBand
            duration_tot: int
            pointsrc: str
        """
        exit_obsstate = False
        duration_tot_req = duration_tot
        band = freqbndobj.rcubands[0]
        rcumode = freqbndobj.rcumodes[0]

        # Get timings
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
        sst_integration = int(acc_cadence)

        # Make sure swlevel=<2
        self.lcu_interface.set_swlevel(2)

        # Set CalServ.conf to dump ACCs:
        self.lcu_interface.acc_mode(enable=True)

        # Boot to swlevel 3 so the calserver service starts
        self.lcu_interface.set_swlevel(3)

        # Beamlet & Subband allocation does not matter here
        # since niether ACC or SST cares
        beamctl_CMD = self.lcu_interface.run_beamctl('0', '255', rcumode, pointing)

        # Run for $duration_tot seconds
        rspctl_CMD = self.lcu_interface.rec_sst(sst_integration, duration_tot)
        self.lcu_interface.stop_beam()

        # Switch back to normal state i.e. turn-off ACC dumping:
        self.lcu_interface.acc_mode(enable=False)

        if exit_obsstate:
            self.lcu_interface.set_swlevel(0)

        # Transfer data from LCU to DAU
        obsdatetime_stamp = self.stationdriver.get_data_timestamp(ACC=True)
        accsrcfiles = self.lcu_interface.ACCsrcDir + "/*.dat"
        acc_destfolder = \
            os.path.join(self.stationdriver.LOFARdataArchive, 'acc',
                         '{}_{}_rcu{}_dur{}'.format(self.lcu_interface.stnid,
                                                    obsdatetime_stamp, rcumode,
                                                    duration_tot))
        if int(rcumode) > 3:
            acc_destfolder += "_"+pointsrc
        acc_destfolder += "_acc"
        if os.path.isdir(acc_destfolder):
            print("Appropriate directory exists already (will put data here)")
        else:
            print("Creating directory "+acc_destfolder+" for ACC "+str(duration_tot)
                  + " s rcumode="+str(rcumode)+" calibration")
            os.mkdir(acc_destfolder)

        # Move ACC dumps to storage
        self.stationdriver.movefromlcu(accsrcfiles, acc_destfolder)
        accdestfiles = os.listdir(acc_destfolder)

        # - Create project header
        sesinfo_acc = copy.deepcopy(self.stationdriver.stnsesinfo)
        acc_integration = 1.0
        sesinfo_acc.set_obsfolderinfo('acc', obsdatetime_stamp, band, acc_integration,
                                     duration_tot, pointing)
        sesinfo_acc.write_session_header(acc_destfolder)

        # - Create header for each ACC file
        for destfile in accdestfiles:
            filedatestr, filetimestr, _ = destfile.split('_', 2)
            filedtstr = '_'.join([filedatestr, filetimestr])
            sesinfo_acc.new_obsinfo()
            sesinfo_acc.obsinfos[-1].setobsinfo_fromparams('acc', filedtstr, beamctl_CMD,
                                                           "")
            sesinfo_acc.obsinfos[-1].create_LOFARst_header(acc_destfolder)

        # Move concurrent data to storage
        sesinfo_sst = copy.deepcopy(self.stationdriver.stnsesinfo)
        sesinfo_sst.new_obsinfo()
        obsdatetime_stamp = self.stationdriver.get_data_timestamp()
        sesinfo_sst.obsinfos[-1].setobsinfo_fromparams('sst', obsdatetime_stamp,
                                                       beamctl_CMD, rspctl_CMD)
        bsxSTobsEpoch, sst_destfolder = \
            sesinfo_sst.obsinfos[-1].getobsdatapath(self.stationdriver.LOFARdataArchive)
        self.stationdriver.movefromlcu(self.lcu_interface.lcuDumpDir + "/*",
                                       sst_destfolder, recursive=True)
        sesinfo_sst.obsinfos[-1].create_LOFARst_header(sst_destfolder)
        sesinfo_sst.set_obsfolderinfo('sst', obsdatetime_stamp, band, sst_integration,
                                      duration_tot, pointing)
        sesinfo_sst.write_session_header(sst_destfolder)
        self.lcu_interface.cleanup()

        acc_url = "{}:{}".format(self.stationdriver.get_stnid(), acc_destfolder)
        sst_url = "{}:{}".format(self.stationdriver.get_stnid(), sst_destfolder)
        return None

    def do_bfsnew(self, freqbndobj, duration_tot, pointsrc, starttime='NOW'):
        """Record BeamFormed Streams (BFS)."""

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
        pointing = modeparms.normalizebeamctldir(pointsrc)
        caltabinfo = self.lcu_interface.getCalTableInfo(modeparms.band2rcumode(band))

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
        beamctl_CMD = self.lcu_interface.run_beamctl(beamletIDs, subbandNrs, band,
                                                     pointing)
        rcu_setup_CMD = self.lcu_interface.rcusetup(bits, attenuation)
        nw = datetime.datetime.utcnow()
        timeleft = st - nw
        if timeleft.total_seconds() < 0.:
            starttimestr = nw.strftime("%Y-%m-%dT%H:%M:%S")
        print("(Beam started) Time left before recording: {}".format(
            timeleft.total_seconds()))

        REC = True
        if REC == True:
            bf_data_dir = self.stationdriver.bf_data_dir
            port0 = self.stationdriver.bf_port0
            stnid = self.lcu_interface.stnid
            bfbackend.rec_bf_streams(starttimestr, duration_tot, lanes, band, bf_data_dir,
                                     port0, stnid)
        else:
            print("Not recording")
            time.sleep(duration_tot)
        sys.stdout.flush()
        self.lcu_interface.stop_beam()
        headertime = datetime.datetime.strptime(starttimestr, "%Y-%m-%dT%H:%M:%S"
                                                ).strftime("%Y%m%d_%H%M%S")
        stnsesinfo = copy.deepcopy(self.stationdriver.stnsesinfo)
        stnsesinfo.new_obsinfo()
        stnsesinfo.obsinfos[-1].setobsinfo_fromparams('bfs', headertime, beamctl_CMD,
                                                      rcu_setup_CMD, caltabinfo)
        bsxSTobsEpoch, datapath = stnsesinfo.obsinfos[-1].getobsdatapath(
            self.stationdriver.LOFARdataArchive)
        print("Creating BFS destination folder on DPU:\n{}".format(datapath))
        os.mkdir(datapath)
        stnsesinfo.obsinfos[-1].create_LOFARst_header(datapath)
        integration = None
        stnsesinfo.set_obsfolderinfo('bfs', headertime, band, integration, duration_tot,
                                     pointing)
        stnsesinfo.write_session_header(datapath)
        self.stationdriver.halt_observingstate_when_finished = shutdown  # Necessary due to forking
        data_url = "{}:{}".format(self.stationdriver.get_stnid(), datapath)
        return None
