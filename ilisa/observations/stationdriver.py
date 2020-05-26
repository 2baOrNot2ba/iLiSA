#!/usr/bin/python
"""Package that provides functions for setting up and running observations
via station controller objects. This package knows about the data archive and
should not run anything directly on LCU."""


import time
import datetime
import subprocess
import os
import multiprocessing

import ilisa.observations.directions
import ilisa.observations.lcuinterface as stationcontrol
import ilisa.observations.modeparms as modeparms


class StationDriver(object):
    """StationDriver is a client type class that allows one to observe with LCU
    and record data and metadata from these observations on a DRU."""

    def checkobservingallowed(self):
        """Check whether observations are allowed. This occurs is someone else
        is using the station. Note that if mockrun then this method will return True."""
        serviceuser = self._lcu_interface.who_servicebroker()

        if serviceuser is None or serviceuser == self._lcu_interface.user:
            stationswitchmode = self._lcu_interface.getstationswitchmode()
            if stationswitchmode == 'local':
                return True
            else:
                print("Warning: Station is not in stand-alone mode.")
                return False
        else:
            print("Warning: Someone else ({}) is using LCU".format(serviceuser))
            print("         (You are running as {})".format(self._lcu_interface.user))
            return False

    def __init__(self, accessconf_lcu, accessconf_dru, mockrun=False):
        """Initialize a StationDriver object, which has access to a station via
        a LCUInterface object configured with setting given by accessconf dict.
        """

        self.mockrun = mockrun
        self._lcu_interface = stationcontrol.LCUInterface(accessconf_lcu)
        if self.mockrun:
            self._lcu_interface.DryRun = self.mockrun

        self.LOFARdataArchive = accessconf_dru['LOFARdataArchive']
        self.bf_data_dir =      accessconf_dru['BeamFormDataDir']
        self.bf_logfile =       accessconf_dru['BeamFormLogFile']
        self.tbbraw2h5cmd =     accessconf_dru['TBBraw2h5Cmd']
        self.tbbh5dumpdir =     accessconf_dru['TBBh5dumpDir']
        bf_ports = self.get_laneports()
        self.bf_port0 = int(bf_ports[0])
        # Path to folder that will contain scans:
        self.scanpath = os.path.join(self.LOFARdataArchive, 'Scans')
        self.exit_check = True
        self.halt_observingstate_when_finished = False
        self.cleanup()

    def isin_observingstate(self):
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
            print("Running warmup beam... @ {}".format(datetime.datetime.utcnow()))
            self.streambeams(modeparms.FrequencyBand('10_90'), '0.,1.5707963,AZELGEO')
            self._lcu_interface.stop_beam()
            print("Finished warmup beam... @ {}".format(datetime.datetime.utcnow()))

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
                if int(swlevel) != 0:
                    print("Warning: You are leaving station in swlevel {} != 0"
                          .format(swlevel))

    def movefromlcu(self, source, dest, recursive=False):
        """Move file(s) off LCU to DPU."""
        if not os.path.exists(dest):
            os.makedirs(dest)
        movecmd = "scp"
        if recursive:
            movecmd += " -r"
        fullcmd = movecmd +" " + self._lcu_interface.lcuURL + ":" + source + " " + dest
        cmdprompt = "on DPU>"
        if self._lcu_interface.verbose:
            print("{} {}".format(cmdprompt, fullcmd))
        subprocess.call(fullcmd, shell=True)
        # # Override DryRun temporarily to really remove source files
        # dryrun = self._lcu_interface.DryRun
        # self._lcu_interface.DryRun = False
        # self._lcu_interface.rm(source)
        # self._lcu_interface.DryRun = dryrun

    def getdatalist(self):
        dd_dir, acc_dir = self._lcu_interface.getdatalist()
        return dd_dir, acc_dir

    def get_data_timestamp(self, order=0, ACC=False):
        """Get timestamp of datafiles on LCU. order is the temporal order of the data
        file, order=0 is oldest, order=-1 is newest."""
        dd_dir, acc_dir = self._lcu_interface.getdatalist()
        # Assumes only one file in datadump dir with
        # format YYYYmmdd_HHMMSS_[bsx]st.dat
        if not ACC:
             the_dir = dd_dir
        else:
            the_dir = acc_dir
        obsdate, obstime, obssuff = the_dir[order].split('_', 2)
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
        if forced:
            # Override DryRun temporarily to really remove source files
            dryrun = self._lcu_interface.DryRun
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

    def _rcusetup(self, bits, attenuation):
        """Setup RCUs on LCU."""
        rcu_setup_cmds = self._lcu_interface.rcusetup(bits, attenuation)
        return rcu_setup_cmds

    def _run_beamctl(self, beamlets, subbands, band, anadigdir):
        """Run beamctl command on LCU."""
        beamctl_cmd = self._lcu_interface.run_beamctl(beamlets, subbands, band, anadigdir)
        return beamctl_cmd

    def streambeams(self, freqbndobj, pointing, recDuration=float('inf'),
                    attenuation=0, DUMMYWARMUP=False):
        """Form beams with station."""
        bits = freqbndobj.bits
        if DUMMYWARMUP:
            print("Warning warnup not currently implemented")
        rcu_setup_cmd = self._lcu_interface.rcusetup(bits, attenuation)
        beamctl_cmds = []
        for bandbeamidx in range(len(freqbndobj.rcumodes)):
            antset = freqbndobj.antsets[bandbeamidx]
            rcumode = freqbndobj.rcumodes[bandbeamidx]
            beamletIDs = freqbndobj.beamlets[bandbeamidx]
            subbands =  freqbndobj.sb_range[bandbeamidx]
            rcusel = freqbndobj.rcusel[bandbeamidx]
            beamctl_main = self._lcu_interface.run_beamctl(beamletIDs, subbands,
                                                           rcumode, pointing, rcusel)
            beamctl_cmds.append(beamctl_main)
        return rcu_setup_cmd, beamctl_cmds

    def stop_beam(self):
        """Turn off beam."""
        self._lcu_interface.stop_beam()

    def rec_bst(self, integration, duration):
        """Record BST data."""
        rspctl_cmd = self._lcu_interface.rec_bst(integration, duration)
        # beamlet statistics also generate empty *01[XY].dat so remove:
        #self._lcu_interface.rm(self._lcu_interface.lcuDumpDir + "/*01[XY].dat")
        return rspctl_cmd

    def rec_sst(self, integration, duration):
        """Record SST data."""
        rspctl_cmd = self._lcu_interface.rec_sst(integration, duration)
        return rspctl_cmd

    def rec_xst(self, subband, integration, duration):
        """Record XST data."""
        rspctl_cmd = self._lcu_interface.rec_xst(subband, integration, duration)
        return rspctl_cmd

    def _waittoboot(self, starttime, pause=0):
        """Before booting, wait until time given by starttime which includes
        a pause . """
        nw = datetime.datetime.utcnow()
        #st = datetime.datetime.strptime(starttime, "%Y-%m-%dT%H:%M:%S")

        maxminsetuptime = datetime.timedelta(seconds=105 + pause)  # Longest minimal time
        # before observation
        # start to set up (This includes a dummy beam warmup)
        d = (starttime - maxminsetuptime) - nw
        timeuntilboot = d.total_seconds()
        if timeuntilboot < 0.:
            timeuntilboot = 0
        print("Will boot to observe state after " + str(timeuntilboot) + " seconds...")
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
        nrpages = str(int(duration_scan*2*Nqfreq/1024))  # One page is 1024 samples.
                                                    # Normal sampling frequency
                                                    # is 200MHz.
        self._lcu_interface.run_tbbctl(select='0:15,16:31,32:47', storage='lofarA1')
        self._lcu_interface.run_tbbctl(select='48:63,64:79,80:95', storage='lofarA2')
        self._lcu_interface.run_tbbctl(select='96:111,112:127,128:143', storage='lofarA3')
        self._lcu_interface.run_tbbctl(select='144:159,160:175,176:191', storage='lofarA4')

        self._lcu_interface.run_tbbctl(cepdelay=str(udpdelay))

        self._lcu_interface.run_tbbctl(select='0:15', readall=nrpages, backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='48:63', readall=nrpages, backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='96:111', readall=nrpages, backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='144:159', readall=nrpages, backgroundJOB='locally')

        self._lcu_interface.run_tbbctl(select='16:31', readall=nrpages, backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='64:79', readall=nrpages, backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='112:127', readall=nrpages, backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='160:175', readall=nrpages, backgroundJOB='locally')

        self._lcu_interface.run_tbbctl(select='32:47', readall=nrpages, backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='80:95', readall=nrpages, backgroundJOB='locally')
        self._lcu_interface.run_tbbctl(select='128:143', readall=nrpages, backgroundJOB='locally')
        # Last one is not put in background so the parent process blocks
        # until finished.
        self._lcu_interface.run_tbbctl(select='176:191', readall=nrpages) # backgroundJOB='locally')

    def do_tbb(self, duration_scan, band, start_after=0, observer="", project=""):
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
        """Determine scan ID.
        The ID is the MJD of the data file stamp time or when the beam started (datetime).
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

# TBBh5dumpDir="/home/tobia/lofar/data/tbb/h5/"
# "/mnt/lane0/TBB/" #Should end with "/" character.
# h5fileprefix="L"

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
    """Wait until datetime. If datetime is 'now' then this is interpreted as current
     time."""
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