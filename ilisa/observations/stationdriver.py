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
import ilisa.observations.beamformedstreams.bfbackend as bfbackend

class ScanMeta(object):
    """Class containing metadata of a scan."""
    def __init__(self, sesspath, bfdsesdumpdir, scanrecs, scan_id=None):
        self.sesspath = sesspath
        self.bfdsesdumpdir = bfdsesdumpdir
        self.scanrecs = scanrecs
        self.scan_id = scan_id

    def write_scanrecs(self):
        """Write the scanrec for each recorded ldat."""
        for ldat in self.scanrecs.keys():
            try:
                scanrecpath = self.scanrecs[ldat].path
            except (AttributeError, KeyError):
                scanrecpath = None
            if scanrecpath:
                self.scanrecs[ldat].write(scanrecpath)
        return scanrecpath


class StationDriver(object):
    """observing.Session class is a client type class typically running on the
    compute node."""

    def cleanup(self):
        pass

    def checkobservingallowed(self):
        """Check whether observations are allowed. This occurs is someone else
        is using the station. Note that if mockrun then this method will return True."""
        serviceuser = self.lcu_interface.who_servicebroker()

        if serviceuser is None or serviceuser == self.lcu_interface.user:
            stationswitchmode = self.lcu_interface.getstationswitchmode()
            if stationswitchmode == 'local':
                return True
            else:
                print("Warning: Station is not in stand-alone mode.")
                return False
        else:
            print("Warning: Someone else ({}) is using LCU".format(serviceuser))
            print("         (You are running as {})".format(self.lcu_interface.user))
            return False

    def __init__(self, accessconf_lcu, accessconf_dru, mockrun=False,
                 goto_observingstate_when_starting=False):
        """Initialize a StationDriver object, which has access to a station via
        a LCUInterface object configured with setting given by accessconf dict.
        When goto_observingstate_when_starting is True, boot the station to swlevel 3,
        but only if observations are allowed on the station.
        """

        self.mockrun = mockrun
        if self.mockrun:
            accessconf_lcu['DryRun'] = True  # mockrun overrides DryRun
        self.lcu_interface = stationcontrol.LCUInterface(accessconf_lcu)

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
        # Check whether the station is being used by someone else:
        if goto_observingstate_when_starting:
            try:
                self.goto_observingstate()
            except RuntimeError:
                raise RuntimeError('Observations not allowed on this station')

    def goto_observingstate(self, warmup=False):
        """Put station into the main observing state.

        Parameters
        ----------

        warmup: bool, optional
                Do a beam warmup after reaching swlevel 3.
        """
        if not self.checkobservingallowed():
            raise RuntimeError('Observations not allowed')
        swlevel_changed = self.lcu_interface.set_swlevel(3)
        if swlevel_changed and warmup:
            # Dummy or hot beam start: (takes about 10sec)
            # This seems necessary: first beamctl after going to swlevel 3
            # seems to crash.
            print("Running warmup beam... @ {}".format(datetime.datetime.utcnow()))
            self.streambeams(modeparms.FrequencyBand('10_90'), '0.,1.5707963,AZELGEO')
            self.lcu_interface.stop_beam()
            print("Finished warmup beam... @ {}".format(datetime.datetime.utcnow()))

    def halt_observingstate(self):
        """Halt observing state on station."""
        if self.checkobservingallowed():
            self.lcu_interface.set_swlevel(0)
            self.cleanup()
            # Cleanup any data left on LCU.
            self.lcu_interface.cleanup()

    def __del__(self):
        """May shutdown observation mode on station."""
        if self.halt_observingstate_when_finished:
            self.halt_observingstate()
        elif self.exit_check:
            if self.checkobservingallowed():
                swlevel = self.lcu_interface.get_swlevel()
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
        fullcmd = movecmd +" " + self.lcu_interface.lcuURL + ":" + source + " " + dest
        cmdprompt = "on DPU>"
        if self.lcu_interface.verbose:
            print("{} {}".format(cmdprompt, fullcmd))
        subprocess.call(fullcmd, shell=True)
        # Override DryRun temporarily to really remove source files
        dryrun = self.lcu_interface.DryRun
        self.lcu_interface.DryRun = False
        self.lcu_interface.rm(source)
        self.lcu_interface.DryRun = dryrun

    def get_data_timestamp(self, order=0, ACC=False):
        """Get timestamp of datafiles on LCU. order is the temporal order of the data
        file, order=0 is oldest, order=-1 is newest."""
        dd_dir, acc_dir = self.lcu_interface.getdatalist()
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
        return self.lcu_interface.stnid

    def get_laneports(self):
        """Return the UDP ports the 4 data lanes are sent to on the DRU."""
        if self.mockrun:
            return 1,2,3,4
        rspdriver_conf = self.lcu_interface.get_RSPDriver_conf()
        port0 = rspdriver_conf['RSPDriver']['LANE_00_DSTPORT']
        port1 = rspdriver_conf['RSPDriver']['LANE_01_DSTPORT']
        port2 = rspdriver_conf['RSPDriver']['LANE_02_DSTPORT']
        port3 = rspdriver_conf['RSPDriver']['LANE_03_DSTPORT']
        return port0, port1, port2, port3

    def streambeams(self, freqbndobj, pointing, recDuration=float('inf'),
                    attenuation=0, DUMMYWARMUP=False):
        """Form beams with station."""
        bits = freqbndobj.bits
        if DUMMYWARMUP:
            print("Warning warnup not currently implemented")
        rcu_setup_cmd = self.lcu_interface.rcusetup(bits, attenuation)
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
        return rcu_setup_cmd, beamctl_cmds

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

    def _setup_tof(self, elemsOn=modeparms.elOn_gILT):
        """Setup (HBA) tiling off mode."""
        # NOTE: LCU must be in swlevel=2 to run SEPTON!
        self.lcu_interface.set_swlevel(2)
        # self.stationcontroller.turnoffElinTile_byTile(elemsOn) # Alternative
        self.lcu_interface.turnoffElinTile_byEl(elemsOn)
        septonconf = modeparms.elementMap2str(elemsOn)
        return septonconf

    def _setupTBBs(self):
        """Setup transient buffer boards and start recording."""
        print("In setupTBBs: Freeing TBBs")
        self.lcu_interface.run_tbbctl(free=True)
        print("In setupTBBs: Setting TBB transient mode on rspctl")
        self.lcu_interface.run_rspctl(tbbmode='transient')
        time.sleep(1)
        print("In setupTBBs: Allocating TBBs")
        self.lcu_interface.run_tbbctl(alloc=True)
        time.sleep(1)
        print("In setupTBBs: Setting TBB transient mode on tbbctl")
        self.lcu_interface.run_tbbctl(mode='transient')
        time.sleep(1)
        print("In setupTBBs: Start TBB recording")
        self.lcu_interface.run_tbbctl(record=True)
        print("In setupTBBs: Finished setting up TBBs & started recording")

    def _freezeTBBdata(self):
        """Stop TBB recording."""
        print("In freezeTBBdata: Stopping TBB recording")
        self.lcu_interface.run_tbbctl(stop=True)
        print("In freezeTBBdata: Stopping any dummy beam")
        self.lcu_interface.stop_beam()

    def _startTBBdataStream(self, duration_scan):
        """Stream duration_scan seconds of TBB data out of the LCU to
        datataking node."""
        # Set delay between subsequent frames. One delay unit is 5us.
        udpdelay = 500  # (Previously 100)
        Nqfreq = 100e6
        nrpages = str(int(duration_scan*2*Nqfreq/1024))  # One page is 1024 samples.
                                                    # Normal sampling frequency
                                                    # is 200MHz.
        self.lcu_interface.run_tbbctl(select='0:15,16:31,32:47', storage='lofarA1')
        self.lcu_interface.run_tbbctl(select='48:63,64:79,80:95', storage='lofarA2')
        self.lcu_interface.run_tbbctl(select='96:111,112:127,128:143', storage='lofarA3')
        self.lcu_interface.run_tbbctl(select='144:159,160:175,176:191', storage='lofarA4')

        self.lcu_interface.run_tbbctl(cepdelay=str(udpdelay))

        self.lcu_interface.run_tbbctl(select='0:15', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='48:63', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='96:111', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='144:159', readall=nrpages, backgroundJOB='locally')

        self.lcu_interface.run_tbbctl(select='16:31', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='64:79', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='112:127', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='160:175', readall=nrpages, backgroundJOB='locally')

        self.lcu_interface.run_tbbctl(select='32:47', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='80:95', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='128:143', readall=nrpages, backgroundJOB='locally')
        # Last one is not put in background so the parent process blocks
        # until finished.
        self.lcu_interface.run_tbbctl(select='176:191', readall=nrpages) # backgroundJOB='locally')

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
