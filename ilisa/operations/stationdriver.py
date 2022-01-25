#!/usr/bin/python
"""Package that provides functions for setting up and running operations
via one LCUinterface instance and one DRUinterface instance.
This package knows about the data archive and
should not run anything directly on LCU."""
import shutil
import threading

import time
import datetime
import subprocess
import os
import warnings
from pathlib import Path
import multiprocessing

import ilisa.operations
import ilisa.operations.directions as directions
from ilisa.operations.lcuinterface import LCUInterface
from ilisa.operations.druinterface import DRUinterface
import ilisa.operations.modeparms as modeparms
import ilisa.operations.data_io as data_io
from ilisa.pipelines import bfbackend


class StationDriver(object):
    """StationDriver is a client type class that allows one to observe with LCU
    and record data and metadata from these operations on a DRU."""

    def is_observingallowed(self):
        """
        Check whether operations are allowed.

        Observing is not allowed if someone else is using the station,
        or if station is not in local mode, or if a beamctl is running.
        Note that if mockrun then this method will return True.

        Returns
        -------
        bool
        """
        serviceuser = self._lcu_interface.who_servicebroker()

        if serviceuser is None or serviceuser == self._lcu_interface.user:
            stationswitchmode = self._lcu_interface.getstationswitchmode()
            if stationswitchmode == 'local':
                if not self._lcu_interface.is_beam_on():
                    return True
                return False
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
        try:
            self._lcu_interface = LCUInterface(accessconf_lcu)
        except AssertionError as a_err:
            self._lcu_interface = None
            raise a_err
        bf_ports = self.get_laneports()
        self.bf_port0 = bf_ports[0]
        self._dru_interface = DRUinterface(accessconf_dru, bf_ports)
        if self.mockrun:
            self._lcu_interface.DryRun = self.mockrun

        self.LOFARdataArchive = accessconf_dru['LOFARdataArchive']
        self.bf_data_dir =      accessconf_dru['BeamFormDataDir']
        self.tbbraw2h5cmd =     accessconf_dru['TBBraw2h5Cmd']
        self.tbbh5dumpdir =     accessconf_dru['TBBh5dumpDir']
        # Path where bsx and acc data from LCU can be staged on driver node
        # (aka Central Control Unit, ccu) before further processing
        self.ccu_cache_bsx = os.path.join(self.LOFARdataArchive, 'Incoming',
                                          'BSX_data')
        self.ccu_cache_acc = os.path.join(self.LOFARdataArchive, 'Incoming',
                                          'ACC_data')
        # Path of symbolic link to latest scan
        self.link2latest = os.path.join(self.LOFARdataArchive, 'Latest')
        # ID of current scan
        self.scan_id = None
        # Path to folder that will contain scans:
        self.scanpath = os.path.join(self.LOFARdataArchive, 'Scans')
        # Path to folder that will contain the scan:
        self.scanpath_scdat = None  # Set under self.scanpath in setup_scan()
        self.exit_check = True
        self.halt_observingstate_when_finished = False
        # Initialize currently running beamctl_cmds & rcuctl_cmds
        self.rcusetup_cmds = []
        self.beamctl_cmds = []
        # and likewise the last executed commands
        self.last_rcusetup_cmds = []
        self.last_beamctl_cmds = []
        # Initialize septonconf setting; implies tile-off (tof) mode
        self.septonconf = None
        # Initialize scanresult
        ## Structure: {'rec': [], 'acc'|'bfs'|'bsx': ScanRecInfo,
        ##             'scan_id': str, 'scanpath_scdat': str}
        self.scanresult = {'rec': []}
        # Initialize beamstart time
        self.beamstart = None

    def is_inobservingstate(self):
        """Check if station is in main observing state for user.
        Returns True if it is, else False.
        """
        if not self.is_observingallowed():
            return False
        swlevel = self._lcu_interface.get_swlevel()
        if swlevel == 3:
            return True
        if swlevel == 2:
            # If in tile-off mode, count it as observing state
            if self.septonconf:
                return True
            else:
                return False
        else:
            return False

    def goto_observingstate(self, warmup=False):
        """Put station into the main observing state.

        Parameters
        ----------

        warmup: bool, optional
                Do a beam warmup after reaching swlevel 3.
        """
        if not self.is_observingallowed():
            raise RuntimeError('Observations not allowed')
        self._lcu_interface.cleanup()  # Could be leftovers from previous runs
        # Clean up locally also
        for ccu_cache in [self.ccu_cache_acc, self.ccu_cache_bsx]:
            for file in os.scandir(ccu_cache):
                os.remove(file.path)

        swlevel_changed = self._lcu_interface.set_swlevel(3)
        if swlevel_changed and warmup:
            # Dummy or hot beam start: (takes about 10sec)
            # This seems necessary: first beamctl after going to swlevel 3
            # seems to crash.
            print("Running warmup beam... @ {}".format(
                datetime.datetime.utcnow()))
            self.streambeams(modeparms.FreqSetup('10_90'),
                             '0.,1.5707963,AZELGEO')
            self._lcu_interface.stop_beam()
            print("Finished warmup beam... @ {}".format(
                datetime.datetime.utcnow()))

    def halt_observingstate(self):
        """Halt observing state on station."""
        if self.is_observingallowed():
            self._lcu_interface.set_swlevel(0)
            # Cleanup any data left on LCU.
            self._lcu_interface.cleanup()

    def __del__(self):
        """
        Delete this object

        May shutdown observation mode on station.
        """
        if self._lcu_interface:
            # If LCU is connected, do the following:
            # Stop any hanging beams running (can happen if an Exception occurs)
            self.stop_beam()
            if self.halt_observingstate_when_finished:
                self.halt_observingstate()
            elif self.exit_check:
                if self.is_observingallowed():
                    swlevel = self._lcu_interface.get_swlevel()
                    if swlevel != 0:
                        print("Warning: You are leaving station in swlevel {} != 0"
                              .format(swlevel))

    def movefromlcu(self, source, dest, recursive=False):
        """Move file(s) off LCU to DRU."""
        if not os.path.exists(dest):
            os.makedirs(dest)
        move_cmdline = ["scp", "-3"]
        if recursive:
            move_cmdline.append("-r")
        src_arg = self._lcu_interface.lcuURL + ":" + source
        move_cmdline.append(src_arg)
        dst_arg = self._dru_interface.url + ":" + dest
        move_cmdline.append(dst_arg)
        cmdprompt = "spawn on driver>"
        if self._lcu_interface.verbose:
            print("{} {}".format(cmdprompt, " ".join(move_cmdline)))
        proc = subprocess.Popen(move_cmdline)
        proc.wait()  # Since cleanup() after might zap data before completion
        self._lcu_interface._rm(source)

    def transfer_accs(self, update=10):
        """\
        Transfer new acc files to DRU as long as acc2stop in not set

        Parameters
        ----------
        update : int
            Look for file every update seconds.
        """
        acc_dur = 512
        wait = acc_dur
        while not self.acc2stop.is_set():
            time.sleep(wait)
            wait = update
            filetimestamps = self.get_datafiletimes(acc=True)
            if filetimestamps:
                wait = acc_dur
                accfilepath = os.path.join(self._lcu_interface.ACCsrcDir,
                                           filetimestamps[0]+'*.dat')
                self.movefromlcu(accfilepath, self.ccu_cache_acc)

    def get_datafiletimes(self, acc=False):
        """\
        Get filetime names of datafiles on LCU sort chronologically
        """
        dd_dir, acc_dir = self._lcu_interface.getdatalist()
        # Assumes files in datadump dir have
        # format YYYYmmdd_HHMMSS_[bsx]st_[rcu???|00[X|Y]].dat
        if not acc:
            the_dir = dd_dir
        else:
            the_dir = acc_dir
        _filetimenames = set()
        for filename in the_dir:
            obsdate, obstime, _obssuff = filename.split('_', 2)
            obsdatetime_stamp = obsdate+'_'+obstime
            _filetimenames.add(obsdatetime_stamp)
        filetimenames = sorted(list(_filetimenames))
        return filetimenames

    def get_stnid(self):
        """Return the station id that this StationDriver is managing."""
        return self._lcu_interface.stnid

    def get_laneports(self):
        """
        Return the UDP ports the 4 data lanes are sent to on the DRU.

        Returns
        -------
        port0: int
            1st port
        port1: int
            2nd port
        part2: int
            3rd port
        part3: int
            4th port

        Notes
        -----
        Seen as a string, the standard format for UDP ports are
        '1'+stnid+port_ordinal, where stnid is a three digit number
        and port_ordinal is 0,1,2 or 3.
        """
        if self.mockrun:
            return 4346, 4347, 4348, 4349  # Return default ports
        rspdriver_conf = self._lcu_interface.get_RSPDriver_conf()
        port0 = int(rspdriver_conf['RSPDriver']['LANE_00_DSTPORT'])
        port1 = int(rspdriver_conf['RSPDriver']['LANE_01_DSTPORT'])
        port2 = int(rspdriver_conf['RSPDriver']['LANE_02_DSTPORT'])
        port3 = int(rspdriver_conf['RSPDriver']['LANE_03_DSTPORT'])
        return port0, port1, port2, port3

    def get_lcuDumpDir(self):
        """Get LCU dump directory from LCU."""
        return self._lcu_interface.lcuDumpDir

    def start_acc_scan(self, duration_tot_req):
        """Start recording ACC"""
        self.scanresult['rec'].append('acc')
        self.scanresult['acc'] = data_io.ScanRecInfo()
        self.scanresult['acc'].set_stnid(self.get_stnid())
        self.scanresult['acc'].set_caltabinfos([])
        # Also duration of ACC sweep since each sb is 1 second.
        dur1acc = modeparms.TotNrOfsb  # Duration of one ACC
        interv2accs = 7  # time between end of 1 ACC and start of next one
        acc_cadence = dur1acc + interv2accs  # =519s between start of 2 ACC
        (nraccs, timrest) = divmod(duration_tot_req, acc_cadence)
        if timrest > dur1acc:
            nraccs += 1
        duration_tot = nraccs * acc_cadence - interv2accs
        if duration_tot != duration_tot_req:
            print("""Note: will use total duration {}s to fit with ACC
                              cadence.""".format(duration_tot))

        # Make sure swlevel=<2
        self._lcu_interface.set_swlevel(2)

        # Set CalServ.conf to dump ACCs:
        self._lcu_interface.acc_mode(enable=True)

        # Boot to swlevel 3 so the calserver service starts
        self._lcu_interface.set_swlevel(3)

        # Possibly make mock acc statistics:
        if self.mockrun:
            self._lcu_interface.mockstatistics('acc', 1.0, duration_tot)

        # Start thread to transfer new acc files
        self.acc2stop = threading.Event()
        self._transfer_acc_thread = threading.Thread(target=self.transfer_accs)
        self._transfer_acc_thread.start()

        return duration_tot

    def stop_acc_scan(self, duration_tot, freqsetup):
        """Stop ACC scan"""
        self.acc2stop.set()  # Send acc2stop event to transfer_accs
        self._transfer_acc_thread.join()  # Wait for transfer_accs

        self._lcu_interface.set_swlevel(2)
        self._lcu_interface.acc_mode(enable=False)
        self._lcu_interface.set_swlevel(3)

        # Create obsinfo each ACC file
        acc_files = os.listdir(self.ccu_cache_acc)
        for acc_file in acc_files:
            obsid, _ = acc_file.split('_acc_')
            rspctl_cmds = []  # ACC doesn't have any rspctl cmds
            ldatinfo_acc = data_io.LDatInfo('acc', self.last_rcusetup_cmds,
                                            self.last_beamctl_cmds, rspctl_cmds)
            ldatinfo_acc.filenametime = obsid
            self.scanresult['acc'].add_obs(ldatinfo_acc)

        # Set scanrecinfo
        pointing = ldatinfo_acc.direction
        acc_integration = 1.0
        self.scanresult['acc'].set_scanrecparms('acc', freqsetup.arg,
                                                duration_tot, pointing,
                                                acc_integration)
        self.scanresult['acc'].set_scanpath(self.scanpath_scdat)
        scanrecpath = self.scanresult['acc'].get_scanrecpath()

        # Create destination folder for data ot end storage
        if os.path.exists(scanrecpath):
            print("Dest directory exists already (will put data here)")
        else:
            print("Creating directory " + scanrecpath + " for ACC "
                  + str(duration_tot) + " s band=" + str(freqsetup.rcumodes[0]))
            os.makedirs(scanrecpath)

        # Move ACC data files and any metadata to end storage
        for _file in os.scandir(self.ccu_cache_acc):
            shutil.move(_file.path, scanrecpath)

    def set_caltable(self, which):
        """Select a calibration table on LCU to use."""
        self._lcu_interface.selectCalTable(which)

    def get_caltableinfos(self, rcumodes):
        """\
        Get Calibration Table Info.

        Parameters
        ----------
        rcumodes : list of rcumodes

        Returns
        -------
        caltabinfos : list of caltabinfos
        """
        caltabinfos = []
        if not self.septonconf:  # calservice not running when septonconf
            for rcumode in rcumodes:
                caltabinfo = self._lcu_interface.get_calinfo(rcumode)
                caltabinfos.append(caltabinfo)
        return caltabinfos

    def _load_user_rcu_disable_list(self, rcumode):
        """\
        Load user specified list of rcus to disable

        The disabled files are in the iLiSA config folder, under
        <stnid>/DISABLED/disabled-mode[3|5|7].txt.
        The text file is just a comma separated list of integers with no spaces.

        This method is analogous to the rcu_disable_list() method of LCUinterface,
        but is defined by user on DRU.

        Parameters
        ----------
        rcumode : int
            The selected spw.

        Returns
        -------
        disabledrcus : list
            List of RCUs the user wishes to disable in special file.
            If the disabled file does not exist, an list with a empty str is
            returned.
        """
        path2disableddir = ilisa.operations.user_conf_dir
        filename = os.path.join(path2disableddir, self._lcu_interface.stnid,
                                "DISABLED",
                                "disabled-mode{}.txt".format(rcumode))
        fp = Path(filename)
        if fp.is_file():
            filecontents = fp.read_text().rstrip()
        else:
            filecontents = ""
        disabledrcus = filecontents.split(',')
        return disabledrcus

    def allowed_rcus(self, rcumode):
        """\
        Get allowed rcus as a flag argument for rcumode

        Currently, the allowed RCUs are determined from the user defined
        list of RCUs to disable.

        Parameters
        ----------
        rcumode : int
            Selected rcumode.

        Returns
        -------
        enabledrcuflagstr : str
            String in CLI flag format of RCUs to enable.
        """
        nrofrcus = modeparms.nrofrcus
        # disabledrcu_lcu = self._lcu_interface.rcu_disable_list(rcumode)
        disabledrcu_usr = self._load_user_rcu_disable_list(rcumode)
        disabledrcu_tot = disabledrcu_usr
        if disabledrcu_tot[0] != '':
            disabledrcus = [int(rcustr) for rcustr in disabledrcu_tot]
            allrcus = range(nrofrcus)
            enabledrcus = [rcu for rcu in allrcus if rcu not in disabledrcus]
            enabledrcuflagstr = modeparms.list2seqarg(enabledrcus)
        else:
            enabledrcuflagstr = "0:{}".format(nrofrcus - 1)
        return enabledrcuflagstr

    def selected_rcus(self, desired_rcus, allowed_rcus):
        """\
        Selected RCUs based on desired and allowed

        Parameters
        ----------
        desired_rcus : str
            CLI flag formatted list of desired RCUs.
        allowed_rcus : str
            CLI flag formatted list of allowed RCUs.

        Returns
        -------
        rcusel : str
            CLI flag formatted str of selected RCUs.
        """
        rcus_allowed_set = set(modeparms.seqarg2list(allowed_rcus))
        rcus_desired_set = set(modeparms.seqarg2list(desired_rcus))
        rcu_list = list(rcus_desired_set.intersection(rcus_allowed_set))
        rcusel = modeparms.list2seqarg(rcu_list)
        return rcusel

    def _rcusetup(self, bits, attenuation, mode=None):
        """Setup RCUs on LCU."""
        rcusetup_cmds = self._lcu_interface.rcusetup(bits, attenuation,
                                                     mode=mode)
        self.rcusetup_cmds = rcusetup_cmds
        return rcusetup_cmds

    def _run_beamctl(self, beamlets, subbands, band, anadigdir, rcusel='all'):
        """Run beamctl command on LCU."""
        beamctl_cmd = self._lcu_interface.run_beamctl(beamlets, subbands,
                                                      band, anadigdir, rcusel)
        self.beamctl_cmds.append(beamctl_cmd)
        waittime = 0  # 11
        print("Waiting {}s for beam to settle...".format(waittime))
        time.sleep(waittime)  # Wait for beam to settle
        return beamctl_cmd

    def streambeams(self, freqbndobj, direction, attenuation=0,
                    dummywarmup=False):
        """\
        Form beams with station
        """
        if dummywarmup:
            print("Warning: warmup not currently implemented")
        bits = freqbndobj.bits
        rcuctl_cmds = self._rcusetup(bits, attenuation)
        beamctl_cmds = []
        bmlt_pntr = 0
        totnrbmlts = 0
        for bandbeamidx in range(len(freqbndobj.rcumodes)):
            _antset = freqbndobj.antsets[bandbeamidx]
            rcumode = freqbndobj.rcumodes[bandbeamidx]
            subbands = freqbndobj.subbands_spw[bandbeamidx]
            beamlets, bmlt_pntr, nrbmlts \
                = modeparms.alloc_beamlets(subbands, bmlt_pntr)
            totnrbmlts += nrbmlts  # TODO Check total # beamlets allowed
            # Select RCUs
            rcusel = self.selected_rcus(freqbndobj.rcusel[bandbeamidx],
                                        self.allowed_rcus(rcumode))
            # Run beamctl
            beamctl_main = self._run_beamctl(beamlets, subbands, rcumode,
                                             direction, rcusel)
            beamctl_cmds.append(beamctl_main)
        self.beamstart = datetime.datetime.utcnow()
        return rcuctl_cmds, beamctl_cmds

    def stop_beam(self):
        """Turn off beam."""
        self._lcu_interface.stop_beam()
        # Save last commands for e.g. ACC
        self.last_rcusetup_cmds = self.rcusetup_cmds
        self.last_beamctl_cmds = self.beamctl_cmds
        self.rcusetup_cmds = []
        self.beamctl_cmds = []

    def start_bsx_scan(self, bsxtype, freqsetup, duration_tot, integration=1.0,
                       duration_file=None):
        """\
        Start BSX scanrec

        Parameters
        ----------
        bsxtype : str
            One of either 'bst', 'sst' or 'xst'.
        freqsetup : FreqSetup
            The frequency setup to use.
        duration_tot : float
            Total duration of scan in seconds.
        integration : float
            Integration time in seconds.
        duration_file : float
            Duration of the bsx recorded file. If None and xst is not freq.
            swept, then it will be set to the total duration, duration_tot.
        """
        rcusetup_cmds = self.rcusetup_cmds
        beamctl_cmds = self.beamctl_cmds

        caltabinfos = []
        sweep_sbs = []  # list of subbands to be swept through
        if not duration_file:
            duration_file = duration_tot
        else:
            duration_file = int(duration_file)

        # Record statistic for duration_tot seconds
        if bsxtype == 'bst':
            caltabinfos = self.get_caltableinfos(freqsetup.rcumodes)
            sweep_sbs = [None]
        elif bsxtype == 'sst':
            sweep_sbs = [None]
        elif bsxtype == 'xst':
            for sb_rcumode in freqsetup.subbands_spw:
                sweep_sbs += modeparms.seqarg2list(sb_rcumode)
        nrsbs2sweep = len(sweep_sbs)
        if nrsbs2sweep > 1:
            duration_file = integration
        # TODO Consider that the desired duration is not the
        #  same as actual duration, if multiple calls to rspctl stats.
        #  xst step takes about 1s extra
        #  bst step takes about 7s extra
        (rep, _rst) = divmod(duration_tot, duration_file * nrsbs2sweep)
        rep = int(rep)
        if rep == 0:
            duration_tot = duration_file * nrsbs2sweep
            warnings.warn(
                "Total duration too short for 1 full repetition."
                "Increasing total duration to {}s.".format(duration_tot))
            rep = 1
        # Set scanrecinfo
        self.scanresult['rec'].append('bsx')
        self.scanresult['bsx'] = data_io.ScanRecInfo()
        self.scanresult['bsx'].set_stnid(self.get_stnid())
        self.scanresult['bsx'].set_caltabinfos(caltabinfos)
        self.scanresult['bsx'].set_scanrecparms(bsxtype,
                                                freqsetup.arg,
                                                duration_tot,
                                                self.pointing_spec['direction'],
                                                integration)
        # Repeat rep times (freq sweep)
        for _itr in range(rep):
            # Sweep through subbands the sweep_sbs list
            for xst_subband in sweep_sbs:
                # Record data
                rspctl_cmds = \
                    self._lcu_interface.run_rspctl_statistics(
                        bsxtype, integration, duration_file, xst_subband)
                ldatinfo = data_io.LDatInfo(bsxtype, rcusetup_cmds,
                                            beamctl_cmds, rspctl_cmds,
                                            septonconf=self.septonconf)
                ft_last = self.get_datafiletimes()[-1]
                ldatinfo.filenametime = ft_last
                self.movefromlcu(self.get_lcuDumpDir() + ldatinfo.filenametime
                                 + '*.dat', self.ccu_cache_bsx)
                yield ldatinfo
        yield None

    def stop_bsx_scan(self, ldatinfos):
        """\
        Stop BSX scan recording
        """
        for ldatinfo in ldatinfos:
            self.scanresult['bsx'].add_obs(ldatinfo)
        # Move data to archive
        self.scanresult['bsx'].set_scanpath(self.scanpath_scdat)
        scanrecpath = self.scanresult['bsx'].get_scanrecpath()
        if not os.path.exists(scanrecpath):
            os.makedirs(scanrecpath)
        for _file in os.listdir(self.ccu_cache_bsx):
            _src = os.path.join(self.ccu_cache_bsx, _file)
            shutil.move(_src, scanrecpath)

    def start_bfs_scan(self, starttime, freqsetup, duration_tot,
                       compress=False):
        """\
        Start recording BFS data

        Parameters
        ----------
        starttime : str
            Start date-time of this recording
        freqsetup : FreqSetup
            Frequency setup
        duration_tot : float
            Total duration in seconds of recording
        compress : bool
            Is recording to be compressed?

        Returns
        -------
        ldatinfo_bfs : LDatInfo
            LDatInfo for this recording.
        bfsdatapaths : list
            List of paths to the recorded BFS data.
        bfslogpaths : list
            List of paths to logs of the recorded BFS data.
        """
        caltabinfos = self.get_caltableinfos(freqsetup.rcumodes)
        rspctl_cmds = []  # BFS doesn't use rspctl cmds
        ldatinfo_bfs = data_io.LDatInfo('bfs', self.rcusetup_cmds,
                                        self.beamctl_cmds, rspctl_cmds)
        scanpath_bfdat = os.path.join(self.bf_data_dir, self.scan_id)
        lanesalloc = modeparms.getlanes(freqsetup.subbands_spw,
                                        freqsetup.bits, freqsetup.nrlanes)
        self.lanes = tuple(lanesalloc.keys())
        if self._lcu_interface.stnid == 'UK902':
            # FIXME
            self.lanes = (0, 1)  # UK902 only has 2 lanes
        # Select only ports for lanes to be used
        laneports = tuple(self.get_laneports()[i] for i in self.lanes)
        datafiles, logfiles = \
            self._dru_interface._rec_bf_proxy(laneports, duration_tot,
                                              scanpath_bfdat,
                                              starttime=starttime,
                                              compress=compress,
                                              band=freqsetup.rcubands[0],
                                              stnid=self.get_stnid())
        bfsdatapaths = []
        bfslogpaths = []
        for lane in self.lanes:
            datafileguess = datafiles.pop()
            dumplogname = logfiles.pop()
            if not datafileguess:
                _outdumpdir, _outarg, datafileguess, dumplogname = \
                    bfbackend.bfsfilepaths(lane, starttime,
                                           modeparms.band2rcumode(freqsetup.rcubands[0]),
                                           scanpath_bfdat,
                                           self.bf_port0,
                                           self.get_stnid())
            bfsdatapaths.append(datafileguess)
            bfslogpaths.append(dumplogname)
        # Duration of BFS not determinable via LCU commands so add this by hand
        ldatinfo_bfs.duration_subscan = duration_tot
        # Set scanrecinfo
        self.scanresult['rec'].append('bfs')
        self.scanresult['bfs'] = data_io.ScanRecInfo()
        self.scanresult['bfs'].set_stnid(self.get_stnid())
        # No integration for BFS
        self.scanresult['bfs'].set_scanrecparms('bfs', freqsetup.arg,
                                                duration_tot,
                                                ldatinfo_bfs.direction,
                                                integration=None)
        self.scanresult['bfs'].set_caltabinfos(caltabinfos)
        return ldatinfo_bfs, bfsdatapaths, bfslogpaths

    def stop_bfs_scan(self, ldatinfo_bfs, rectime, bfsdatapaths, bfslogpaths):
        """Stop BFS scan"""

        file_dt_name = rectime.strftime("%Y%m%d_%H%M%S")
        ldatinfo_bfs.filenametime = file_dt_name
        self.scanresult['bfs'].add_obs(ldatinfo_bfs)

        # Make a project folder for BFS data
        self.scanresult['bfs'].set_scanpath(self.scanpath_scdat)
        scanrecpath = self.scanresult['bfs'].get_scanrecpath()
        # Create BFS destination folder on DPU:
        os.makedirs(scanrecpath)
        if self._dru_interface.hostname == 'localhost':
            # Make soft links to actual BFS files and move logs to scanrec
            # folder
            for lane in self.lanes:
                if bfsdatapaths[lane] is not None:
                    _basename = os.path.basename(bfsdatapaths[lane])
                    _lnkname = os.path.join(scanrecpath, _basename)
                    os.symlink(bfsdatapaths[lane], _lnkname)
                if bfslogpaths[lane] is not None:
                    shutil.move(bfslogpaths[lane], scanrecpath)

    def _waittoboot(self, starttime, pause=0):
        """Before booting, wait until time given by starttime which includes
        a pause . """
        nw = datetime.datetime.utcnow()
        # st = datetime.datetime.strptime(starttime, "%Y-%m-%dT%H:%M:%S")

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
        self.septonconf = modeparms.elementMap2str(elemsOn)
        return self.septonconf

    def stop_tof(self):
        """Stop tiling off mode."""
        self._lcu_interface.set_swlevel(3)
        self.septonconf = None

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
        pointing = directions.std_pointings('Z')
        freqband = modeparms.FreqSetup(band)
        # FreqSetup obtained from band spec sets 8 bit mode,
        # so create a new FreqSetup object with only center frequency
        freqlo, freqhi = freqband.edgefreqs()
        freq0 = (freqlo+freqhi)/2.0
        actualfb = modeparms.FreqSetup(freq0)
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

    def setup_scan(self, bsx_stat=False, bfs=False, acc=False, beamstarted=None,
                   scanroot=None):
        """Setup a Scan"""
        self.goto_observingstate()
        self.scanresult = {'rec': []}
        if acc:
            self.scanresult['rec'].append('acc')
            self.scanresult['acc'] = data_io.ScanRecInfo()
        if bsx_stat:
            self.scanresult['rec'].append('bsx')
            self.scanresult['bsx'] = data_io.ScanRecInfo()
        if bfs:
            self.scanresult['rec'].append('bfs')
            self.scanresult['bfs'] = data_io.ScanRecInfo()
        if scanroot:
            self.scanpath = scanroot
        self.scan_id = self.get_scanid(beamstarted)
        self.scanresult['scan_id'] = self.scan_id
        # Work out where station-correlated data should be stored:
        self.scanpath_scdat = os.path.join(self.scanpath, self.scan_id)
        self.scanresult['scanpath_scdat'] = self.scanpath_scdat
        # and create the directory: (may not have ldat if no rec but will have
        # info files)
        os.makedirs(self.scanpath_scdat)
        # Make symbolic link from 'Latest' folder to latest scan folder
        if bsx_stat:
            os.remove(self.link2latest)
            os.symlink(self.ccu_cache_bsx, self.link2latest)

    def get_scanid(self, beamstarted=None):
        """
        Determine scan ID.
        
        The ID is the MJD of the data file stamp time or when the beam started
        (datetime).
        """
        try:
            scan_dt = datetime.datetime.strptime(self.get_datafiletimes()[-1],
                                                 "%Y%m%d_%H%M%S")
        except Exception:
            if beamstarted is None:
                beamstarted = datetime.datetime.utcnow()
            scan_dt = beamstarted
        scan_mjd_id = modeparms.dt2mjd(scan_dt)
        scan_id = "scan_{}".format(scan_mjd_id)
        return scan_id


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
    Wait until requested datetime starttime_req with a margin in seconds.

    Parameters
    ----------
    starttime_req : datetime
        Requested datetime to wait until. If datetime is 'NOW' or 'ASAP',
        then this is interpreted as current time or as soon as possible.
    margin : timedelta
        Margin of timedelta before actually returning at starttime_req.

    Returns
    -------
    starttime : datetime
        datetime this function waited until (excludes margin).
    """
    now = datetime.datetime.utcnow()
    if starttime_req == "NOW" or starttime_req == "ASAP":
        starttime = now
    else:
        starttime = starttime_req
    timeleft = (starttime - margin) - now
    secondsleft = int(timeleft.total_seconds())
    if secondsleft < 0:
        secondsleft = 0
    print("Waiting {}s until {}.".format(secondsleft, starttime))
    time.sleep(secondsleft)
    return starttime
