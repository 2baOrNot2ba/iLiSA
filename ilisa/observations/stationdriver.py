#!/usr/bin/python
"""Package that provides functions for setting up and running observations
via one LCUinterface instance and one DRUinterface instance.
This package knows about the data archive and
should not run anything directly on LCU."""
import shutil
import sys

import time
import datetime
import subprocess
import os
import warnings
from pathlib import Path
import multiprocessing

import ilisa.observations.directions
from ilisa.observations.lcuinterface import LCUInterface
from ilisa.observations.druinterface import DRUinterface
import ilisa.observations.modeparms as modeparms
import ilisa.observations.dataIO as dataIO
from ilisa.pipelines import bfbackend


class StationDriver(object):
    """StationDriver is a client type class that allows one to observe with LCU
    and record data and metadata from these observations on a DRU."""

    def checkobservingallowed(self):
        """
        Check whether observations are allowed. This occurs is someone else is
        using the station. Note that if mockrun then this method will return
        True.
        """
        serviceuser = self._lcu_interface.who_servicebroker()

        if serviceuser is None or serviceuser == self._lcu_interface.user:
            stationswitchmode = self._lcu_interface.getstationswitchmode()
            if stationswitchmode == 'local':
                return True
            else:
                print("Warning: Station is not in stand-alone mode.")
                return False
        else:
            print("""Warning: Someone else ({}) is using LCU
                              (You are running as {})"""\
                               .format(serviceuser, self._lcu_interface.user))
            return False

    def __init__(self, accessconf_lcu, accessconf_dru, mockrun=False):
        """\
        Initialize a StationDriver object, which has access to a station via
        a LCUInterface object configured with setting given by accessconf dict.
        """

        self.mockrun = mockrun
        self._lcu_interface = LCUInterface(accessconf_lcu)
        bf_ports = self.get_laneports()
        self.bf_port0 = int(bf_ports[0])
        self.dru_interface = DRUinterface(accessconf_dru, bf_ports)
        if self.mockrun:
            self._lcu_interface.DryRun = self.mockrun

        self.LOFARdataArchive = accessconf_dru['LOFARdataArchive']
        self.bf_data_dir =      accessconf_dru['BeamFormDataDir']
        self.bf_logfile =       accessconf_dru['BeamFormLogFile']
        self.tbbraw2h5cmd =     accessconf_dru['TBBraw2h5Cmd']
        self.tbbh5dumpdir =     accessconf_dru['TBBh5dumpDir']
        # Path to folder that will contain scans:
        self.scanpath = os.path.join(self.LOFARdataArchive, 'Scans')
        self.exit_check = True
        self.halt_observingstate_when_finished = False
        self.cleanup()
        # Initialize beamctl_cmds & rcu_setup_cmds
        self.rcu_setup_cmds = []
        self.beamctl_cmds = []

    def is_in_observingstate(self):
        """Check if station is in main observing state for user.
        Returns True if it is else False.
        """
        if not self.checkobservingallowed():
            return False
        swlevel = self._lcu_interface.get_swlevel()
        if swlevel == 3:
            return True
        else:
            return False

    def goto_observingstate(self, warmup=False):
        """Put station into the main observing state.

        Parameters
        ----------

        warmup: bool, optional
                Do a beam warmup after reaching swlevel 3.
        """
        if not self.checkobservingallowed():
            raise RuntimeError('Observations not allowed')
        swlevel_changed = self._lcu_interface.set_swlevel(3)
        if swlevel_changed and warmup:
            # Dummy or hot beam start: (takes about 10sec)
            # This seems necessary: first beamctl after going to swlevel 3
            # seems to crash.
            print("Running warmup beam... @ {}"\
                .format(datetime.datetime.utcnow()))
            self.streambeams(modeparms.FrequencyBand('10_90'),
                             '0.,1.5707963,AZELGEO')
            self._lcu_interface.stop_beam()
            print("Finished warmup beam... @ {}"\
                .format(datetime.datetime.utcnow()))

    def halt_observingstate(self):
        """Halt observing state on station."""
        if self.checkobservingallowed():
            self._lcu_interface.set_swlevel(0)
            # Cleanup any data left on LCU.
            self.cleanup()

    def __del__(self):
        """May shutdown observation mode on station."""
        if self.halt_observingstate_when_finished:
            self.halt_observingstate()
        elif self.exit_check:
            if self.checkobservingallowed():
                swlevel = self._lcu_interface.get_swlevel()
                if swlevel != 0:
                    print("Warning: You are leaving station in swlevel {} != 0"
                          .format(swlevel))

    def movefromlcu(self, source, dest, recursive=False):
        """Move file(s) off LCU to DPU."""
        if not os.path.exists(dest):
            os.makedirs(dest)
        move_cmdline = ["scp"]
        if recursive:
            move_cmdline.append("-r")
        src_arg = self._lcu_interface.lcuURL + ":" + source
        move_cmdline.append(src_arg)
        dst_arg = dest
        move_cmdline.append(dst_arg)
        cmdprompt = "spawn on DPU>"
        if self._lcu_interface.verbose:
            print("{} {}".format(cmdprompt, " ".join(move_cmdline)))
        proc = subprocess.Popen(move_cmdline)
        proc.wait()  # Since cleanup() after might zap data before completion
        # # Override DryRun temporarily to really remove source files
        # dryrun = self._lcu_interface.DryRun
        # self._lcu_interface.DryRun = False
        # self._lcu_interface.rm(source)
        # self._lcu_interface.DryRun = dryrun

    def getdatalist(self):
        dd_dir, acc_dir = self._lcu_interface.getdatalist()
        return dd_dir, acc_dir

    def get_data_timestamp(self, order=0, ACC=False):
        """
        Get timestamp of datafiles on LCU. order is the temporal order of
        the data file, order=0 is oldest, order=-1 is newest.
        """
        dd_dir, acc_dir = self._lcu_interface.getdatalist()
        # Assumes only one file in datadump dir with
        # format YYYYmmdd_HHMMSS_[bsx]st.dat
        if not ACC:
             the_dir = dd_dir
        else:
            the_dir = acc_dir
        obsdate, obstime, _obssuff = the_dir[order].split('_', 2)
        obsdatetime_stamp = obsdate+'_'+obstime
        return obsdatetime_stamp

    def get_stnid(self):
        """Return the station id that this StationDriver is managing."""
        return self._lcu_interface.stnid

    def get_laneports(self):
        """Return the UDP ports the 4 data lanes are sent to on the DRU."""
        if self.mockrun:
            return 1,2,3,4
        rspdriver_conf = self._lcu_interface.get_RSPDriver_conf()
        port0 = rspdriver_conf['RSPDriver']['LANE_00_DSTPORT']
        port1 = rspdriver_conf['RSPDriver']['LANE_01_DSTPORT']
        port2 = rspdriver_conf['RSPDriver']['LANE_02_DSTPORT']
        port3 = rspdriver_conf['RSPDriver']['LANE_03_DSTPORT']
        return port0, port1, port2, port3

    def get_ACCsrcDir(self):
        """Get ACC src directory from LCU."""
        return self._lcu_interface.ACCsrcDir

    def get_lcuDumpDir(self):
        """Get LCU dump directory from LCU."""
        return  self._lcu_interface.lcuDumpDir

    def cleanup(self, forced=False):
        """Clean up on LCU."""
        dryrun = self._lcu_interface.DryRun
        if forced:
            # Override DryRun temporarily to really remove source files
            self._lcu_interface.DryRun = False
        self._lcu_interface.cleanup()
        if forced:
            self._lcu_interface.DryRun = dryrun

    def acc_mode(self, enable=True, mock_dur=None):
        """Enable or Disable ACC mode."""
        if enable:
            # Make sure swlevel=<2
            self._lcu_interface.set_swlevel(2)

            # Set CalServ.conf to dump ACCs:
            self._lcu_interface.acc_mode(enable=True)

            # Boot to swlevel 3 so the calserver service starts
            self._lcu_interface.set_swlevel(3)

            # Possibly make mock acc statistics:
            if self.mockrun and mock_dur is not None:
                self._lcu_interface.mockstatistics('acc', 1.0, mock_dur)
        else:
            self._lcu_interface.set_swlevel(2)
            self._lcu_interface.acc_mode(enable=False)
            self._lcu_interface.set_swlevel(3)

    def set_caltable(self, which):
        """Select a calibration table on LCU to use."""
        self._lcu_interface.selectCalTable(which)

    def get_caltableinfo(self, rcumode):
        """Get Calibration Table Info."""
        caltabinfo = self._lcu_interface.getCalTableInfo(rcumode)
        return caltabinfo

    def __load_user_rcu_disable_list(self, rcumode):
        """Load list of rcus that the user wants to disable.
        This is analogous to the rcu_disable_list() method of LCUinterface,
        but is defined by user on DRU.
        """
        path2disableddir = ilisa.observations.user_conf_dir
        filename = os.path.join(path2disableddir, self._lcu_interface.stnid,
                                "DISABLED",
                                "disabled-mode{}.txt".format(rcumode))
        fp = Path(filename)
        if fp.is_file():
            filecontents = fp.read_text()
        else:
            filecontents = ""
        disabledrcus = filecontents.split(',')
        return disabledrcus

    def allowed_rcus(self, rcumode):
        """Get allowed rcus as a flag argument."""
        nrofrcus = modeparms.nrofrcus
        #disabledrcu_lcu = self._lcu_interface.rcu_disable_list(rcumode)
        disabledrcu_usr = self.__load_user_rcu_disable_list(rcumode)
        disabledrcu_tot = disabledrcu_usr
        if disabledrcu_tot[0] != '':
            disabledrcus = [int(rcustr) for rcustr in disabledrcu_tot]
            allrcus = range(nrofrcus)
            enabledrcus = [rcu for rcu in allrcus if rcu not in disabledrcus]
            enabledrcuflagstr = str(enabledrcus[0]) + ":"
            for rcuidx in range(1, len(enabledrcus) - 1):
                if enabledrcus[rcuidx] - enabledrcus[rcuidx - 1] != 1:
                    enabledrcuflagstr += str(enabledrcus[rcuidx - 1]) + "," \
                                         + str(enabledrcus[rcuidx]) + ":"
            enabledrcuflagstr += str(enabledrcus[-1])
        else:
            enabledrcuflagstr = "0:{}".format(nrofrcus - 1)
        return enabledrcuflagstr

    def _rcusetup(self, bits, attenuation):
        """Setup RCUs on LCU."""
        rcu_setup_cmds = self._lcu_interface.rcusetup(bits, attenuation)
        self.rcu_setup_cmds = rcu_setup_cmds
        return rcu_setup_cmds

    def _run_beamctl(self, beamlets, subbands, band, anadigdir, rcusel='all'):
        """Run beamctl command on LCU."""
        beamctl_cmd = self._lcu_interface.run_beamctl(beamlets, subbands,
                                                      band, anadigdir, rcusel)
        self.beamctl_cmds.append(beamctl_cmd)
        return beamctl_cmd

    def streambeams(self, freqbndobj, pointing, recDuration=float('inf'),
                    attenuation=0, DUMMYWARMUP=False):
        """Form beams with station."""
        if DUMMYWARMUP:
            print("Warning warnup not currently implemented")
        bits = freqbndobj.bits
        rcu_setup_cmd = self._rcusetup(bits, attenuation)
        beamctl_cmds = []
        for bandbeamidx in range(len(freqbndobj.rcumodes)):
            _antset = freqbndobj.antsets[bandbeamidx]
            rcumode = freqbndobj.rcumodes[bandbeamidx]
            beamlets = freqbndobj.beamlets[bandbeamidx]
            subbands =  freqbndobj.sb_range[bandbeamidx]
            # Select RCUs
            rcus_allowed_set = set(modeparms.seqarg2list(
                self.allowed_rcus(rcumode)))
            rcus_desired_set = set(modeparms.seqarg2list(
                freqbndobj.rcusel[bandbeamidx]))
            rcu_list = list(rcus_desired_set.intersection(rcus_allowed_set))
            rcusel = ','.join(map(str, list(rcu_list)))
            # Run beamctl
            beamctl_main = self._run_beamctl(beamlets, subbands, rcumode,
                                             pointing, rcusel)
            beamctl_cmds.append(beamctl_main)
        return rcu_setup_cmd, beamctl_cmds

    def stop_beam(self):
        """Turn off beam."""
        self._lcu_interface.stop_beam()
        self.rcu_setup_cmds = []
        self.beamctl_cmds = []

    def rec_bsx(self, bsxtype, integration, duration, subband=0):
        """Record BSX data."""
        rspctl_cmd = self._lcu_interface.run_rspctl_statistics(bsxtype,
                                                               integration,
                                                               duration,
                                                               subband)
        return rspctl_cmd

    def start_scanrec(self, bsxtype, integration, duration, subband=0):
        """\
        Start scan on LCU
        """
        rcu_setup_cmds = self.rcu_setup_cmds
        beamctl_cmds = self.beamctl_cmds
        rspctl_cmd = self._lcu_interface.run_rspctl_statistics(bsxtype,
                                                               integration,
                                                               duration,
                                                               subband)
        obsdatetime_stamp = self.get_data_timestamp(-1)
        ldatinfo = dataIO.LDatInfo(bsxtype, obsdatetime_stamp, beamctl_cmds,
                                    rspctl_cmd, caltabinfos="", septonconf="")
        return ldatinfo
    
    def stop_scanrec(self, ldatinfo, freqbndobj):
        """\
        Stop scan on LCU
        """
        # Set scanrecinfo
        scanresult = {'rec': ['bsx']}
        scanresult['bsx'] = dataIO.ScanRecInfo()
        scanresult['bsx'].set_stnid(self.get_stnid())
        scanresult['bsx'].add_obs(ldatinfo)
        scanresult['bsx'].set_scanrecparms(ldatinfo.ldat_type, freqbndobj.arg,
                                           ldatinfo.duration_scan,
                                           ldatinfo.pointing,
                                           ldatinfo.integration)
        # Move data to archive
        scanresult['bsx'].set_scanpath(self.scanpath)
        scanrecpath = scanresult['bsx'].get_scanrecpath()
        self.movefromlcu(self.get_lcuDumpDir() + "/*.dat", scanrecpath)
        return scanresult

    def _waittoboot(self, starttime, pause=0):
        """Before booting, wait until time given by starttime which includes
        a pause . """
        nw = datetime.datetime.utcnow()
        #st = datetime.datetime.strptime(starttime, "%Y-%m-%dT%H:%M:%S")

        maxminsetuptime = datetime.timedelta(seconds=105 + pause)
        # maxminsetuptime is the Longest minimal time
        # before observation start to set up (This includes a dummy beam warmup)
        d = (starttime - maxminsetuptime) - nw
        timeuntilboot = d.total_seconds()
        if timeuntilboot < 0.:
            timeuntilboot = 0
        print("Will boot to observe state after " + str(timeuntilboot)
              + " seconds...")
        time.sleep(timeuntilboot)
        # From swlevel 0 it takes about 1:30min? to reach swlevel 3
        print("Booting @ {}".format(datetime.datetime.utcnow()))

    def setup_tof(self, elemsOn=modeparms.elOn_gILT):
        """Setup (HBA) tiling off mode."""
        # NOTE: LCU must be in swlevel=2 to run SEPTON!
        self._lcu_interface.set_swlevel(2)
        # self.stationcontroller.turnoffElinTile_byTile(elemsOn) # Alternative
        self._lcu_interface.turnoffElinTile_byEl(elemsOn)
        septonconf = modeparms.elementMap2str(elemsOn)
        return septonconf

    def stop_tof(self):
        """Stop tiling off mode."""
        self._lcu_interface.set_swlevel(2)

    def _setupTBBs(self):
        """Setup transient buffer boards and start recording."""
        print("In setupTBBs: Freeing TBBs")
        self._lcu_interface.run_tbbctl(free=True)
        print("In setupTBBs: Setting TBB transient mode on rspctl")
        self._lcu_interface.run_rspctl(tbbmode='transient')
        time.sleep(1)
        print("In setupTBBs: Allocating TBBs")
        self._lcu_interface.run_tbbctl(alloc=True)
        time.sleep(1)
        print("In setupTBBs: Setting TBB transient mode on tbbctl")
        self._lcu_interface.run_tbbctl(mode='transient')
        time.sleep(1)
        print("In setupTBBs: Start TBB recording")
        self._lcu_interface.run_tbbctl(record=True)
        print("In setupTBBs: Finished setting up TBBs & started recording")

    def _freezeTBBdata(self):
        """Stop TBB recording."""
        print("In freezeTBBdata: Stopping TBB recording")
        self._lcu_interface.run_tbbctl(stop=True)
        print("In freezeTBBdata: Stopping any dummy beam")
        self._lcu_interface.stop_beam()

    def _startTBBdataStream(self, duration_scan):
        """Stream duration_scan seconds of TBB data out of the LCU to
        datataking node."""
        # Set delay between subsequent frames. One delay unit is 5us.
        udpdelay = 500  # (Previously 100)
        Nqfreq = 100e6
        nrpages = str(int(duration_scan*2*Nqfreq/1024))
        # One page is 1024 samples.
        # Normal sampling frequency is 200MHz
        self._lcu_interface.run_tbbctl(select='0:15,16:31,32:47',
                                       storage='lofarA1')
        self._lcu_interface.run_tbbctl(select='48:63,64:79,80:95',
                                       storage='lofarA2')
        self._lcu_interface.run_tbbctl(select='96:111,112:127,128:143',
                                       storage='lofarA3')
        self._lcu_interface.run_tbbctl(select='144:159,160:175,176:191',
                                       storage='lofarA4')

        self._lcu_interface.run_tbbctl(cepdelay=str(udpdelay))

        self._lcu_interface.run_tbbctl(select='0:15', readall=nrpages,
                                       backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='48:63', readall=nrpages,
                                       backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='96:111', readall=nrpages,
                                       backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='144:159', readall=nrpages,
                                       backgroundJOB='locally')

        self._lcu_interface.run_tbbctl(select='16:31', readall=nrpages,
                                       backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='64:79', readall=nrpages,
                                       backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='112:127', readall=nrpages,
                                       backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='160:175', readall=nrpages,
                                       backgroundJOB='locally')

        self._lcu_interface.run_tbbctl(select='32:47', readall=nrpages, 
                                       backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='80:95', readall=nrpages,
                                       backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='128:143', readall=nrpages,
                                       backgroundJOB='locally')
        # Last one is not put in background so the parent process blocks
        # until finished.
        self._lcu_interface.run_tbbctl(select='176:191', readall=nrpages)

    def do_tbb(self, duration_scan, band, start_after=0, observer="",
               project=""):
        """Record duration_scan seconds of TBB data from rcumode.
        """
        observationID = "Null"

        # Prepare for obs program.
        try:
            self.goto_observingstate()
        except RuntimeError as e:
            raise RuntimeError(e)

        # Start a beam
        pointing = ilisa.observations.directions.std_pointings('Z')
        freqband = modeparms.FrequencyBand(band)
        # FrequencyBand obtained from band spec sets 8 bit mode,
        # so create a new FrequencyBand object with only center frequency
        freqlo, freqhi = freqband.edgefreqs()
        freq0 = (freqlo+freqhi)/2.0
        actualfb = modeparms.FrequencyBand(freq0)
        antset = actualfb.antsets[0]
        self.streambeams(actualfb, pointing)

        print("Setting up TBBs")
        self._setupTBBs()

        print("Will start TBB recording in {}s".format(start_after))
        time.sleep(start_after)

        print("Will freeze TBB recording in {}s".format(duration_scan))
        time.sleep(duration_scan)  # Arbitrary time to trigger
        print("Sending trigger to TBBs")
        self._freezeTBBdata()

        # Start data capture process locally
        dalcap = \
            multiprocessing.Process(target=capture_data_DAL1,
                                    args=(self.tbbraw2h5cmd, self.tbbh5dumpdir,
                                          observer, antset, project,
                                          observationID, False))
        dalcap.start()

        print("Streaming {}s of TBB data out of LCU".format(duration_scan))
        self._startTBBdataStream(float(duration_scan))
        dalcap.join()

    def get_scanid(self, beamstarted=None):
        """
        Determine scan ID.
        
        The ID is the MJD of the data file stamp time or when the beam started
        (datetime).
        """
        try:
            scan_dt = datetime.datetime.strptime(self.get_data_timestamp(-1),
                                                 "%Y%m%d_%H%M%S")
        except Exception:
            if beamstarted is None:
                beamstarted = datetime.datetime.utcnow()
            scan_dt = beamstarted
        scan_mjd_id = modeparms.dt2mjd(scan_dt)
        scan_id = "scan_{}".format(scan_mjd_id)
        return scan_id

    def record_scan(self, freqbndobj, duration_tot, pointing,
                    starttime='NOW', rec=(None,), integration=1.0, allsky=False,
                    duration_frq=None):
        """Run a generic scan.

        Parameters
        ----------
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

        Returns
        -------
        scanresult : dict
            Results of scan.
        """
        rec_acc = False
        rec_bfs = False
        bsx_type = None
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

        stnid = self.get_stnid()

        # Initialize scanresult
        scanresult = {'rec': []}

        # Setup Calibration tables on LCU:
        CALTABLESRC = 'default'  # FIXME put this in args
        ## (Only BST uses calibration tables)
        # Choose between 'default' or 'local'
        self.set_caltable(CALTABLESRC)

        # Prepare for obs program.
        try:
            self.goto_observingstate()
        except RuntimeError as e:
            raise RuntimeError(e)

        duration_tot_req = duration_tot
        band = freqbndobj.rcubands[0]
        rcumode = freqbndobj.rcumodes[0]

        if todo_tof:
            septonconf = self.setup_tof()
        else:
            septonconf = None

        if rec_acc:
            scanresult['rec'].append('acc')
            scanresult['acc'] = dataIO.ScanRecInfo()
            scanresult['acc'].set_stnid(stnid)
            # Also duration of ACC sweep since each sb is 1 second.
            dur1acc = modeparms.TotNrOfsb  # Duration of one ACC
            interv2accs = 7  # time between end of one ACC and start of next one
            acc_cadence = dur1acc + interv2accs  # =519s time between start of 2 ACCs
            (nraccs, timrest) = divmod(duration_tot_req, acc_cadence)
            if timrest > dur1acc:
                nraccs += 1
            duration_tot = nraccs * acc_cadence - interv2accs
            if duration_tot != duration_tot_req:
                print("""Note: will use total duration {}s to fit with ACC
                      cadence.""".format(duration_tot))
            self.acc_mode(enable=True, mock_dur=duration_tot)

        # Wait until it is time to start
        starttime_req = starttime
        starttime = waituntil(starttime_req)

        # Necessary since fork creates multiple instances of myobs and each one
        # will call it's __del__ on completion and __del__ shutdown...
        shutdown = self.halt_observingstate_when_finished
        self.halt_observingstate_when_finished = False
        self.exit_check = False

        caltabinfos = [""]
        # Get metadata about caltables to be used
        if not allsky:
            for rcumode in freqbndobj.rcumodes:
                caltabinfo = self.get_caltableinfo(rcumode)
                caltabinfos.append(caltabinfo)

        if pointing is not None:
            # Real beam start:
            print("Now running real beam... @ {}".format(
                datetime.datetime.utcnow()))
            dir_bmctl = ilisa.observations.directions.normalizebeamctldir(pointing)
            rcu_setup_cmd, beamctl_cmds = self.streambeams(freqbndobj,
                                                           dir_bmctl)
            beamstarted = datetime.datetime.utcnow()
            rectime = beamstarted
            scan_id = self.get_scanid(beamstarted)
        else:
            rcu_setup_cmd, beamctl_cmds = "", ""
            rectime = starttime
            scan_id = self.get_scanid()

        if rec_bfs:
            scanresult['rec'].append('bfs')
            scanresult['bfs'] = dataIO.ScanRecInfo()
            scanresult['bfs'].set_stnid(stnid)
            scanpath_bfdat = os.path.join(self.bf_data_dir, scan_id)
            lanes = tuple(freqbndobj.getlanes().keys())
            datafiles, logfiles = self.dru_interface.rec_bf_proxy(rectime,
                                                                  duration_tot,
                                                                  lanes, band,
                                                                  scanpath_bfdat,
                                                                  self.bf_port0,
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
                                               self.bf_port0,
                                               self.get_stnid())
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
                rspctl_cmd = self.rec_bsx(bsx_type, integration,
                                          duration_tot)
                obsdatetime_stamp = self.get_data_timestamp(-1)
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
                            rspctl_cmd = self.rec_bsx(bsx_type, integration,
                                                      duration_frq, subband)
                            obsdatetime_stamp = self.get_data_timestamp(-1)
                            curr_obsinfo =\
                                dataIO.LDatInfo('xst', obsdatetime_stamp,
                                                beamctl_cmds, rspctl_cmd,
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
        self.stop_beam()

        if todo_tof:
            self.stop_tof()

        # Work out where station-correlated data should be stored:
        scanpath_scdat = os.path.join(self.scanpath, scan_id)
        # and create the directory: (may not have ldat if no rec but will have info
        # files)
        os.makedirs(scanpath_scdat)

        if rec_acc:
            # Switch back to normal state i.e. turn-off ACC dumping:
            self.acc_mode(enable=False)

            # Create obsinfo each ACC file
            _, acc_files = self.getdatalist()
            for acc_file in acc_files:
                obsid, _ = acc_file.split('_acc_')
                this_obsinfo = dataIO.LDatInfo('acc', obsid, beamctl_cmds,
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
                print("Creating directory " + scanrecpath + " for ACC "
                      + str(duration_tot) + " s rcumode=" + str(rcumode) + " calibration")
                os.makedirs(scanrecpath)

            # Move ACC dumps to storage
            accsrcfiles = self.get_ACCsrcDir() + "/*.dat"
            self.movefromlcu(accsrcfiles, scanrecpath)

        if bsx_type is not None:
            # Set scanrecinfo
            scanresult['bsx'].set_scanrecparms(bsx_type, freqbndobj.arg,
                                               duration_tot, pointing,
                                               integration, allsky=allsky)
            # Move data to archive
            scanresult['bsx'].set_scanpath(scanpath_scdat)
            scanrecpath = scanresult['bsx'].get_scanrecpath()
            self.movefromlcu(self.get_lcuDumpDir() + "/*.dat",
                             scanrecpath)

        if rec_bfs:
            scanresult['bfs'].set_stnid(self.get_stnid())
            scanresult['bfs'].set_scanrecparms('bfs', band, duration_tot, pointing,
                                               allsky=allsky)
            # Make a project folder for BFS data
            scanresult['bfs'].set_scanpath(scanpath_scdat)
            scanrecpath = scanresult['bfs'].get_scanrecpath()
            # Create BFS destination folder on DPU:
            os.makedirs(scanrecpath)
            if self.dru_interface.hostname == 'localhost':
                # Make soft links to actual BFS files and move logs to scanrec folder
                for lane in lanes:
                    if bfsdatapaths[lane] is not None:
                        os.symlink(bfsdatapaths[lane], os.path.join(scanrecpath,
                                                                    os.path.basename(
                                                                        bfsdatapaths[
                                                                            lane])))
                    if bfslogpaths[lane] is not None:
                        shutil.move(bfslogpaths[lane], scanrecpath)

        self.cleanup()
        # Necessary due to possible forking
        self.halt_observingstate_when_finished = shutdown
        scanresult['scan_id'] = scan_id
        scanresult['scanpath_scdat'] = scanpath_scdat
        return scanresult


def capture_data_DAL1(tbbraw2h5cmd, TBBh5dumpDir, observer, antennaSet,
                      project, observationID, background=False):
    """Start process on DPU to capture streamed TBB data on LCU."""
    if background:
        ground = 'bg'
        grdcmd = '&'
    else:
        ground = 'fg'
        grdcmd = ''
    print("Starting TBBraw2h5 ({})".format(ground))
    if observer == '':
        observer = 'Null'
    if project == '':
        project = 'Null'
    cmdline = ("cd "+TBBh5dumpDir+"; "
                     + " export LD_LIBRARY_PATH=/mnt/old/usr/lib/ ; "
                     + tbbraw2h5cmd
                     + " "
                     + " --observer="+observer
                     + " --antennaSet="+antennaSet
                     + " --project="+project
                     + " --observationID="+observationID
                     + " -P 31664 -P 31665 -P 31666 "
                     + " -P 31667 -P 31668 -P 31669 "
                     + " -P 31670 -P 31671 -P 31672 "
                     + " -P 31673 -P 31674 -P 31675 "
                     + grdcmd)
    print("DPU> {}".format(cmdline))
    subprocess.call(cmdline, shell=True)


def waituntil(starttime_req, margin=datetime.timedelta(seconds=0)):
    """
    Wait until datetime. If datetime is 'now' then this is interpreted as
    current time.
    """
    now = datetime.datetime.utcnow()
    if starttime_req == "NOW":
        starttime = now
    else:
        starttime = starttime_req
    timeleft = (starttime - margin) - now
    secondsleft = int(timeleft.total_seconds())
    if secondsleft < 0:
        secondsleft = 0
    print("Waiting {}s before starting.".format(secondsleft))
    time.sleep(secondsleft)
    return starttime
