#!/usr/bin/python
"""Package that provides functions for setting up and running operations
via one LCUinterface instance and one DRUinterface instance.
This package knows about the data archive and
should not run anything directly on LCU."""
import shutil
import sys
import time
import datetime
import subprocess
import os
from pathlib import Path
import multiprocessing
import argparse
import logging

import ilisa.operations
import ilisa.operations.directions as directions
from ilisa.operations.lcuinterface import LCUInterface
from ilisa.operations.druinterface import DRUinterface
from ilisa.operations._rem_exec import ostimenow
import ilisa.operations.modeparms as modeparms
import ilisa.operations.data_io as data_io

_LOGGER = logging.getLogger(__name__)


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
                _LOGGER.warning("Station is not in stand-alone mode.")
                return False
        else:
            _LOGGER.warning(
                """Someone else ({}) is using LCU (You are running as {})"""
                .format(serviceuser, self._lcu_interface.user))
            return False

    def __init__(self, accessconf_lcu=None, accessconf_dru=None, mockrun=False):
        """\
        Initialize a StationDriver object, which has access to a station via
        a LCUInterface object configured with setting given by accessconf dict.

        Parameters
        ----------
        accessconf_lcu : dict
            LCU's access configuration metadata.
        accessconf_dru : dict
            DRU's access configuration metadata.
        mockrun : bool
            Whether this is mock run. Default False.

        Raises
        ------
        AssertionError
            If LCU is not set up correctly.
        ConnectionError
            If LCU could not be accessed.
        """

        self.mockrun = mockrun
        if not accessconf_lcu or not accessconf_dru:
            accessconf = ilisa.operations.default_access_lclstn_conf()
            if not accessconf_lcu:
                accessconf_lcu = accessconf['LCU']
            if not accessconf_dru:
                accessconf_dru = accessconf['DRU']
        try:
            self._lcu_interface = LCUInterface(accessconf_lcu)
        except (ConnectionError, AssertionError) as err:
            self._lcu_interface = None
            raise err
        bf_ports = self.get_laneports()
        self.bf_port0 = bf_ports[0]
        self._dru_interface = DRUinterface(accessconf_dru, bf_ports)
        if self.mockrun:
            self._lcu_interface.DryRun = self.mockrun

        _lofardatadir = accessconf_dru['LOFARdataArchive']
        # Make _lofardatadir relative to self._dru_root so it can be joined:
        if os.path.isabs(_lofardatadir):
            _lofardatadir = os.path.relpath(_lofardatadir, os.sep)
        self.tbbraw2h5cmd = accessconf_dru['TBBraw2h5Cmd']
        self.tbbh5dumpdir = accessconf_dru['TBBh5dumpDir']

        # Possibly setup sshfs
        self.use_sshfs = True
        if self._dru_interface.hostname == 'localhost':
            self.use_sshfs = False
        self._dru_root = os.sep  # DRU root directory on DRU
        if self.use_sshfs:
            # Setup sshfs to DRU
            ## Let self._dru_root be path to DRU root on CCU
            ##   set as ~/.cache/ilisa/DRU/<dru_hostname>
            self._dru_root = os.path.join(ilisa.operations.USER_CACHE_DIR,
                                          'DRU', self._dru_interface.hostname)
            # Mount it using sshfs
            subprocess.run(['sshfs', self._dru_interface.hostname + ':/',
                            self._dru_root])
        # Set the root path to where the lofar data should be stored
        self.dru_data_root = os.path.join(self._dru_root, _lofardatadir)
        # ID of current scan
        self.scan_id = None
        # Path to folder that will contain scans:
        self.scanpath = os.path.join(self.dru_data_root, 'Scans')
        # Path to folder that will contain the scan:
        self.scanpath_scdat = None  # Set under self.scanpath in init_scan()
        self.exit_check = True
        self.halt_observingstate_when_finished = False
        # Initialize currently running beamctl_cmds & rcuctl_cmds
        self.rcusetup_cmds = []
        self.beamctl_cmds = []
        # and likewise the last executed commands
        # Initialize septonconf setting; implies tile-off (tof) mode
        self.septonconf = None
        # Initialize beamstart time
        self.beamstart = None
        # Initialize field, name of field station pointing at
        self.field = ''

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

        swlevel_changed = self._lcu_interface.set_swlevel(3)
        if swlevel_changed and warmup:
            # Dummy or hot beam start: (takes about 10sec)
            # This seems necessary: first beamctl after going to swlevel 3
            # seems to crash.
            _LOGGER.info("Running warmup beam... @ {}".format(
                datetime.datetime.utcnow()))
            self.streambeams(modeparms.FreqSetup('10_90'),
                             '0.,1.5707963,AZELGEO')
            self._lcu_interface.stop_beam()
            _LOGGER.info("Finished warmup beam... @ {}".format(
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
                        _LOGGER.warning(
                            "You are leaving station in swlevel {} != 0"
                            .format(swlevel))

    def movefromlcu(self, source, dest, recursive=False):
        """Move file(s) off LCU to DRU."""
        move_cmdline = ["scp", "-3"]
        if recursive:
            move_cmdline.append("-r")
        src_arg = self._lcu_interface.url + ":" + source
        move_cmdline.append(src_arg)
        if self.use_sshfs or self._dru_interface.hostname == 'localhost':
            dst_arg = dest
            if not os.path.exists(dest):
                os.makedirs(dest)
        else:
            dst_arg = self._dru_interface.url + ":" + dest
            # TODO: mkdir on DRU
        move_cmdline.append(dst_arg)
        cmdprompt = "spawn on driver>"
        if self._lcu_interface.verbose:
            _LOGGER.info("{} {}".format(cmdprompt, " ".join(move_cmdline)))
        proc = subprocess.Popen(move_cmdline)
        proc.wait()  # Since cleanup() after might zap data before completion
        self._lcu_interface._rm(source)

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

    def get_ccu2dru(self):
        """Get path to DRU root from CCU"""
        return self._dru_root

    def get_bfsdatlogpaths(self):
        """\
        Get BFS data and log file paths
        """
        bfs_data_names, bfs_log_names = self._dru_interface.get_bfs_filenames(
            self.scanpath_bfdat)
        dru_bf_data_dir = self._dru_interface.bf_data_dir
        bfsdatapaths = []
        bfslogpaths = []
        for bfs_data_name in bfs_data_names:
            bfsdatapaths.append(os.path.join(dru_bf_data_dir, bfs_data_name))
        for bfs_log_name in bfs_log_names:
            bfslogpaths.append(os.path.join(dru_bf_data_dir, bfs_log_name))
        return bfsdatapaths, bfslogpaths

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
        path2disableddir = ilisa.operations.USER_CONF_DIR
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
        _LOGGER.info("Waiting {}s for beam to settle...".format(waittime))
        time.sleep(waittime)  # Wait for beam to settle
        return beamctl_cmd

    def streambeams(self, freqbndobj, direction, attenuation=0,
                    dummywarmup=False):
        """\
        Form beams with station
        """
        if dummywarmup:
            _LOGGER.warning("Warmup not currently implemented")
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
        self.rcusetup_cmds = []
        self.beamctl_cmds = []

    def start_acc_scan(self):
        """\
        Start recording ACC

        See also
        --------
        stop_acc_scan : Stops ACC scan.

        Notes
        -----
            ACC stands for Autocorrelation cubes. They are produced by the
            `CalServer` in swlevel 3 as soon as `beamctl` command is running.

        Yields
        ------
        None : if no new ACC dump file exists
        scanrecpath : for 1st ACC dump file
        ldatinfo_acc : once for every new ACC dump file
        """
        # Make sure swlevel=<2
        self._lcu_interface.set_swlevel(2)

        # Set CalServ.conf to dump ACCs:
        self._lcu_interface.acc_mode(enable=True)

        # Boot to swlevel 3 so the calserver service starts
        self._lcu_interface.set_swlevel(3)

        integration = 1.0
        dur1acc = modeparms.ACC_DUR
        firstacc = True
        continue_acc = True
        while continue_acc:
            filetimestamps = self.get_datafiletimes(acc=True)
            if not filetimestamps:
                # There is no ACC on LCU yet
                continue_acc = yield None
                continue
            accfilepath = os.path.join(self._lcu_interface.ACCsrcDir,
                                       filetimestamps[0] + '*.dat')
            # Create obsinfo each ACC file
            rspctl_cmds = []  # ACC doesn't have any rspctl cmds
            ldatinfo_acc = data_io.LDatInfo(
                'acc', self.rcusetup_cmds, self.beamctl_cmds,
                rspctl_cmds)
            ldatinfo_acc.filenametime = filetimestamps[0]
            if firstacc:
                obsfileinfo = {
                    'duration': dur1acc,
                    'filenametime': ldatinfo_acc.filenametime,
                    'integration': integration,
                    'spw': ldatinfo_acc.get_spw(),
                    'ldat_type': 'acc',
                    'pointing': self.field,
                    'sb': [],
                    'station_id': self.get_stnid()
                }
                accfilefolder = data_io.obsfileinfo2filefolder(obsfileinfo)
                scanrecpath = os.path.join(self.scanpath_scdat, accfilefolder)
            self.movefromlcu(accfilepath, scanrecpath)
            ldatinfo_acc.write_ldat_header(scanrecpath)
            if firstacc:
                yield scanrecpath
                firstacc = False
            continue_acc = yield ldatinfo_acc

    def stop_acc_scan(self):
        """\
        Stop ACC scan

        See also
        --------
        start_acc_scan : Starts ACC scan.
        """
        self._lcu_interface.set_swlevel(2)
        self._lcu_interface.acc_mode(enable=False)
        self._lcu_interface.set_swlevel(3)

    def rec_bsx_scan(self, bsxtype, freqsetup, duration_tot, integration=1.0,
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

        Yields
        ------
        ldatinfo : LDatInfo
            Metadata on the LOFAR data recorded in this subscan
        """
        rcusetup_cmds = self.rcusetup_cmds
        beamctl_cmds = self.beamctl_cmds

        sweep_sbs = []  # list of subbands to be swept through
        if not duration_file:
            duration_file = duration_tot
        else:
            duration_file = int(duration_file)

        # Record statistic for duration_tot seconds
        if bsxtype == 'bst':
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
            _LOGGER.warning(
                "Total duration too short for 1 full repetition."
                "Increasing total duration to {}s.".format(duration_tot))
            rep = 1
        filenametime_first = ''
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
                if not filenametime_first:
                    filenametime_first = ldatinfo.filenametime
                    obsfileinfo = {
                        'duration': duration_tot,
                        'filenametime': filenametime_first,
                        'integration': integration,
                        'spw': ldatinfo.get_spw(),
                        'ldat_type': bsxtype,
                        'pointing': self.field,
                        'sb': freqsetup.subbands_spw,
                        'station_id': self.get_stnid()
                    }
                    bsxfilefolder = data_io.obsfileinfo2filefolder(obsfileinfo)
                    scanrecpath = os.path.join(self.scanpath_scdat,
                                               bsxfilefolder)
                    continue_sub = yield scanrecpath
                self.movefromlcu(self.get_lcuDumpDir() + ldatinfo.filenametime
                                 + '*.dat', scanrecpath)
                ldatinfo.write_ldat_header(scanrecpath)
                continue_sub = yield ldatinfo

    def start_bfs_scan(self, starttime, freqsetup, duration_tot,
                       duration_file=None, compress=False):
        """\
        Start recording BFS data

        Parameters
        ----------
        starttime : str
            Start date-time of this recording
        freqsetup : FreqSetup
            Frequency setup
        duration_file : float
            Duration of dumpfiles in seconds. If set to None (default), then
            only make one dumpfile.
        duration_tot : float
            Total duration in seconds of recording
        compress : bool
            Is recording to be compressed?

        Returns
        -------
        ldatinfos_bfs : LDatInfo
            LDatInfo for this recording.
        bfsdatapaths : list
            List of paths to the recorded BFS data.
        bfslogpaths : list
            List of paths to logs of the recorded BFS data.
        """
        #caltabinfos = self.get_caltableinfos(freqsetup.rcumodes)
        rspctl_cmds = []  # BFS doesn't use rspctl cmds
        #ldatinfo_bfs = data_io.LDatInfo('bfs', self.rcusetup_cmds,
        #                                self.beamctl_cmds, rspctl_cmds)
        lanesalloc = modeparms.getlanes(freqsetup.subbands_spw,
                                        freqsetup.bits, freqsetup.nrlanes)
        self.lanes = tuple(lanesalloc.keys())
        if self._lcu_interface.stnid == 'UK902':
            # FIXME
            self.lanes = (0, 1)  # UK902 only has 2 lanes
        # Select only ports for lanes to be used
        laneports = tuple(self.get_laneports()[i] for i in self.lanes)
        datafiles, logfiles = \
            self._dru_interface.start_bf_rec(laneports, duration_tot,
                                             self.scanpath_bfdat,
                                             starttime=starttime,
                                             file_dur=duration_file,
                                             compress=compress,
                                             band=freqsetup.rcubands[0],
                                             stnid=self.get_stnid())
        # bf_rec started, now yield None to sync
        yield None
        firstbfs = True
        continue_bfs = True
        filetime_prev = None
        while continue_bfs:
            bfs_datfiles, bfs_logfiles = self._dru_interface.get_bfs_filenames(
                self.scanpath_bfdat
            )
            filetime_new = None
            if bfs_datfiles:
                port, hostname, startstr, ms, cmprss_suf = \
                    self._dru_interface.parse_bfs_filename(bfs_datfiles[-1])
                filetime_curr = modeparms.timestr2datetime(startstr).strftime(
                    "%Y%m%d_%H%M%S")
                if filetime_curr != filetime_prev:
                    filetime_new = filetime_curr
                    filetime_prev = filetime_curr
            if not filetime_new:
                # There is no new BFS on DRU yet
                continue_bfs = yield None
                continue
            # Create obsinfo each BFS file
            rspctl_cmds = []  # BFS doesn't have any rspctl cmds
            ldatinfo_bfs = data_io.LDatInfo(
                'bfs', self.rcusetup_cmds, self.beamctl_cmds,
                rspctl_cmds)
            # Duration of BFS not determinable via LCU commands
            # so add this by hand
            ldatinfo_bfs.duration_subscan = duration_file
            ldatinfo_bfs.filenametime = filetime_new
            if firstbfs:
                obsfileinfo = {
                    'duration': duration_tot,
                    'filenametime': ldatinfo_bfs.filenametime,
                    'integration': None,
                    'spw': ldatinfo_bfs.get_spw(),
                    'ldat_type': 'bfs',
                    'pointing': self.field,
                    'sb': freqsetup.subbands_spw,
                    'station_id': self.get_stnid()
                }
                # Create BFS destination folder on DRU:
                bfsfilefolder = data_io.obsfileinfo2filefolder(obsfileinfo)
                scanrecpath = os.path.join(self.scanpath_scdat, bfsfilefolder)
                os.makedirs(scanrecpath)
                yield scanrecpath
                firstbfs = False
            ldatinfo_bfs.write_ldat_header(scanrecpath)
            continue_bfs = yield ldatinfo_bfs
        #return ldatinfos_bfs, bfsdatapaths, bfslogpaths, scanrecpath

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
        _LOGGER.info("Will boot to observe state after {}s..."
                     .format(timeuntilboot))
        time.sleep(timeuntilboot)
        # From swlevel 0 it takes about 1:30min? to reach swlevel 3
        _LOGGER.info("Booting")

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

    def init_scan(self, scan_id, scanroot=None, destsubpath_bfs=None,
                  bsx_stat=False):
        # Work out where station-correlated data should be stored:
        if scanroot:
            self.scanpath = scanroot
        self.scanpath_scdat = os.path.join(self.scanpath, scan_id)
        # and create the directory: (may not have ldat if no rec but will have
        # info files)
        os.makedirs(self.scanpath_scdat)
        # Also setup for BFS data:
        bfs_scan_root = self._dru_interface.bf_data_dir
        if destsubpath_bfs:
            bfs_scan_root = os.path.join(bfs_scan_root, destsubpath_bfs)
        self.scanpath_bfdat = os.path.join(bfs_scan_root, scan_id)

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

    def _now_diff(self, remunit='LCU'):
        """\
        Measure difference between now time on LCU and stationdriver
        """
        if remunit=='LCU':
            url = self._lcu_interface.url
        elif remunit=='DRU':
            url = self._dru_interface.url
        stndrv_before = datetime.datetime.utcnow()
        remunit_dattim = ostimenow(url, remunit)
        stndrv_after = datetime.datetime.utcnow()
        # Compute time diff as difference between remote units time
        # and mean of local before and after time
        stndrv_mean = (stndrv_after - stndrv_before)/2+stndrv_before
        timdif = remunit_dattim - stndrv_mean
        return timdif

    def _time2startup_hint(self, whatservice='beam'):
        """\
        Give a guess as to how much time it will take for beam to start

        Parameters
        ----------
        whatservice : str
            The service to lookup the startup time for.
            Defined services are: 'beam' (default), 'boot', 'bst' and 'tof'.

        Returns
        -------
        service_time2startup : datetime.timedelta
            The time it takes for this service to startup in seconds.
        """
        # N.B. These times are valid for current scheme of iLiSA internal
        # commands issued from SE607's SDU to its LCU.
        stnid = self.get_stnid()
        sshcmd_delay = 20  # This is a maximum time e.g. IE613 from OSO
        boot_time = 138  # 77
        boot3_time = 4
        beam_time = 2
        beamsettle_time = 13  # Time it takes after beamctl to get real data
        idle_time = 41
        checkobs_time = 2
        bst_time = 7
        tof_time = 10
        bfs_time = beam_time + 18

        if stnid == 'SE607':
            sshcmd_delay = 0
        service_inittime = 0
        if whatservice == 'boot':
            service_inittime = boot_time
            if self._lcu_interface.get_swlevel() == 3:
                service_inittime = boot3_time
        elif whatservice == 'beam':
            service_inittime = beam_time
        elif whatservice == 'idle':
            service_inittime = idle_time
        elif whatservice == 'checkobs':
            service_inittime = checkobs_time
        elif whatservice == 'bst':
            service_inittime = bst_time
        elif whatservice == 'tof':
            service_inittime = tof_time
        service_time2startup = service_inittime + 0*sshcmd_delay
        return datetime.timedelta(seconds=service_time2startup)


def _is_sshfs_mounted(hostname):
    mounts_pnts = subprocess.run('mount', stdout=subprocess.PIPE
                                 ).stdout.decode().rstrip().split('\n')
    sshfs_mnt_pnts = filter(lambda l : 'fuse.sshfs' in l, mounts_pnts)
    hostname_mnt_pnts = list(filter(lambda l : l.startswith(hostname),
                                    sshfs_mnt_pnts))
    _LOGGER.debug("hostname_mnt_pnts={}".format(hostname_mnt_pnts))
    if hostname_mnt_pnts:
        return True
    else:
        return False


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
    time_at_return : datetime
        datetime this function returned.
    """
    now = datetime.datetime.utcnow()
    if starttime_req == "NOW" or starttime_req == "ASAP":
        starttime = now
    else:
        starttime = starttime_req
    timeleft = (starttime - margin) - now
    secondsleft = int(timeleft.total_seconds())
    if secondsleft < 0:
        lag = now - starttime
        _LOGGER.info("Lagging {}s past {} (margin={})".format(lag,
                                                              starttime,
                                                              margin))
    else:
        _LOGGER.info("Waiting {}s before {} (margin={})".format(secondsleft,
                                                                starttime,
                                                                margin))
        time.sleep(secondsleft)
    time_at_return = datetime.datetime.utcnow()
    return time_at_return


def boot(stndrv):
    """\
    Put station into ready to observe state

    Note: This typically takes 77 s.
    """
    stndrv.goto_observingstate()


def idle(stndrv):
    """\
    Put station into idle state

    Note: This typically takes 41 s.
    """
    stndrv.halt_observingstate()


def handback(stndrv):
    """Handback station to ILT control."""
    idle(stndrv)


def checkobs(stndrv):
    """Check if user can observe on LCU."""
    is_inobsstate = stndrv.is_inobservingstate()
    print("User can observe on station {} now: {}".format(stndrv.get_stnid(),
          is_inobsstate))
    if not is_inobsstate:
        obs_allowed = stndrv.is_observingallowed()
        if obs_allowed:
            reason = "swlevel not 3"
        else:
            reason = "Observing not allowed now"
        print("Reason: {}.".format(reason))


def main_cli():
    """iLiSA adm CLI"""
    cmdln_prsr = argparse.ArgumentParser()
    cmdln_prsr.add_argument('-t', '--time', type=str, default='ASAP',
                            help="Execute at time (format: YYYY-mm-ddTHH:MM:SS)"
                            )
    cmdln_prsr.add_argument('-s', '--station', type=str, default=None,
                            help="Station ID")
    cmdln_prsr.add_argument('-m', '--mockrun', action='store_true',
                            help="Mockrun")
    cmdln_prsr.add_argument('admcmd', help='Admin command')
    args = cmdln_prsr.parse_args(sys.argv[1:])
    accessconf = ilisa.operations.default_access_lclstn_conf(args.station)
    if accessconf==None:
        sys.exit("No default access config for local station {}"
                 .format(args.station))
    stndrv = StationDriver(accessconf['LCU'], accessconf['DRU'],
                           mockrun=args.mockrun)
    starttime = modeparms.timestr2datetime(args.time)
    waituntil(starttime, stndrv._time2startup_hint(args.admcmd))
    # Dispatch admin commands
    if args.admcmd == 'boot':
        boot(stndrv)
    elif args.admcmd == 'idle':
        idle(stndrv)
    elif args.admcmd == 'handback':
        handback(stndrv)
    elif args.admcmd == 'checkobs':
        checkobs(stndrv)


if __name__ == "__main__":
    main_cli()
