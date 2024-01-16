import argparse
import datetime
import os
import shutil
import time
import logging

import ilisa.operations
import ilisa.operations.directions as directions
import ilisa.operations.modeparms as modeparms
import ilisa.operations.data_io as data_io
from ilisa.operations.stationdriver import StationDriver, waituntil

_LOGGER = logging.getLogger(__name__)


class LScan:
    def __init__(self, stndrv, rec_type, freqsetup, pointing_spec, duration_tot,
                 integration, starttime, acc=False, bfs=False, destpath=None,
                 destpath_bfs=None, file_dur=None, scan_id=None):
        """\
        Record a scan of LOFAR station data

        Parameters
        ----------
        stndrv: StationDriver
            StationDriver object.
        rec_type: str
            Main type of LOFAR data to record.
            The possible data types are:
                'bst', 'sst', 'xst', 'tbb', 'dmp', None
            where:
                'bst': Beamlet STatistics data
                'sst': Subband STatistics data
                'xst': Xrosslet STatistics data
                'tbb': Transient Buffer Board data
                'dmp': Only record bfs, do not start a beam
                None: No recording of data.
        freqsetup: FreqSetup
            FreqSetup() instance.
        pointing_spec: dict
            Pointing source direction of scan. Can be a source name or a
            beamctl direction str. If None, then allsky image is implied.
        duration_tot: float
            Requested total duration of scan in seconds.
        integration: float
            Integration time in seconds.
        starttime: str
            The time at which scan should start.
        acc: bool
            If true, record autocovariance-cubes (ACC) files. The actual total
            duration will at most be duration_tot. (Usually it will be
            shorter so it fits within the cadence of whole ACC aquisition,
            which is 512+7=519 seconds). ACC files are the covariance of
            all array elements with each as a function of subband.
        bfs: bool
            If true, record 'bfs': Record BeamFormed Stream data.
        destpath: str
            Destination path for this recording.
        destpath_bfs: str
            Destination path on DRU to where bfs to be recorded.
        file_dur: float
            Duration of file of recorded (subscan) data.
        scanid: str
            ID for this scan.
        """
        self.stndrv = stndrv
        self.rec_type = rec_type
        self.freqsetup = freqsetup
        self.duration_tot = duration_tot
        self.pointing_spec = pointing_spec
        self.integration = integration
        self.starttime = starttime
        self.acc = acc
        self.bfs = bfs
        self.file_dur = file_dur if file_dur else duration_tot
        self.destpath = destpath if destpath else stndrv.scanpath
        self.destpath_bfs = destpath_bfs

        self.dir_bmctl = pointing_spec['direction']  # ilisa.operations.directions.normalizebeamctldir(pointing)
        if self.dir_bmctl:
            if not ilisa.operations.directions.check_directionstr(self.dir_bmctl):
                raise ValueError("Invalid pointing syntax: {}"
                                 .format(self.dir_bmctl))
        beam_needed = False
        if self.bfs or self.rec_type == 'bst':
            # Also ACC needs beam but not included here since it only leads
            # to blank ACC data and conceptually ought not to need a beam
            beam_needed = True
        if beam_needed and not self.dir_bmctl:
            raise ValueError("No pointing, but beam needed")
        if not self.dir_bmctl and self.acc:
            # If ACC set but no beam defined, start a dummy beam
            _LOGGER.info("Setting dummy beam for ACC data.")
            self.dir_bmctl = "0,0,AZELGEO"

        self.stndrv.goto_observingstate()
        # Initialize scanresult
        ## Structure: {'rec': [], 'acc'|'bfs'|'bsx': ScanRecInfo,
        ##             'scan_id': str, 'scanpath_scdat': str}
        self.scanresult = {}
        bsx_stat = modeparms._xtract_bsx(rec_type)
        caltabinfos = []
        if bsx_stat == 'bst' or bfs:
            caltabinfos = self.stndrv.get_caltableinfos(freqsetup.rcumodes)
        self.scanpath = self.destpath
        self.scan_id = scan_id if scan_id else self._gen_scan_id(self.starttime)
        self.scanresult['scan_id'] = self.scan_id

        if rec_type != 'tbb' and rec_type != 'dmp':
            # Initialize scan on station driver
            stndrv.init_scan(self.scan_id, scanroot=self.destpath,
                             destsubpath_bfs=self.destpath_bfs)
            self.scanresult['scanpath_scdat'] = self.stndrv.scanpath_scdat
            self.scanresult['rec'] = []

            if acc:
                self.scanresult['rec'].append('acc')
                self.scanresult['acc'] = data_io.ScanRecInfo(
                        self.stndrv.get_stnid(), caltabinfos, stndrv.mockrun)
                self.scanresult['acc'].set_caltabinfos([])
                self.scanresult['acc'].set_scanrecparms(
                    'acc', freqsetup.freq_arg, duration_tot,
                    self.pointing_spec['direction'], 1.0)
            if bsx_stat:
                self.scanresult['rec'].append('bsx')
                self.scanresult['bsx'] = data_io.ScanRecInfo(
                        self.stndrv.get_stnid(), caltabinfos, stndrv.mockrun)
                self.scanresult['bsx'].set_scanrecparms(
                    bsx_stat, freqsetup.freq_arg, duration_tot,
                    self.pointing_spec['direction'], integration)
            if bfs:
                self.scanresult['rec'].append('bfs')
                self.scanresult['bfs'] = data_io.ScanRecInfo(
                        self.stndrv.get_stnid(), caltabinfos, stndrv.mockrun)
                self.scanresult['bfs'].set_scanrecparms(
                    'bfs', freqsetup.freq_arg, duration_tot,
                    self.pointing_spec['direction'], None)

        self.ldatinfos = []
        self.ldatinfos_acc = []
        self.ldatinfos_bfs = []

    def __iter__(self):
        stndrv = self.stndrv
        rec_type = self.rec_type
        freqsetup = self.freqsetup
        duration_tot = self.duration_tot
        integration = self.integration
        starttime = self.starttime
        file_dur = self.file_dur

        bsx_type = modeparms._xtract_bsx(rec_type)
        stndrv.field = self.pointing_spec.get('source')

        if rec_type != 'tbb' and rec_type != 'dmp':
            if self.acc:
                # ACC needs to be enabled before beam
                # so initialize the subscan generator...
                acc_subscan = stndrv.start_acc_scan(duration_tot)
                # actually start the subscan...
                _ = next(acc_subscan)
                if not self.dir_bmctl:
                    _LOGGER.warning("ACC set but no beam run"
                                    " (ACC data will be blank)")
            if self.dir_bmctl:
                stndrv.streambeams(freqsetup, self.dir_bmctl)
            else:
                stndrv._rcusetup(freqsetup.bits, 0, mode=freqsetup.rcumodes[0])
                if freqsetup.rcumodes[0] > 4:
                    # No pointing for HBA implies tiles-off mode
                   stndrv.setup_tof()
            if self.bfs:
                bfs_subscan = stndrv.start_bfs_scan(starttime, freqsetup,
                                                    duration_tot,
                                                    duration_file=file_dur,
                                                    compress=True)
                bfs_yield = next(bfs_subscan)
            if bsx_type:
                bsx_subscan = stndrv.rec_bsx_scan(bsx_type, freqsetup,
                                                  duration_tot, integration,
                                                  duration_file=file_dur)

            # LCU is now generating data. Yield to signal ready for subscans
            continue_scan = yield
            # Enter subscan data collection loop
            normal_acc_wait = True
            bfs_wait = file_dur
            acc_yield = None
            while continue_scan:
                if bsx_type:
                    if self.ldatinfos == []:
                        self.scanresult['bsx'].scanrecpath = next(bsx_subscan)
                    try:
                        ldatinfo = bsx_subscan.send(continue_scan)
                    except StopIteration:
                        continue_scan = False
                    else:
                        self.ldatinfos.append(ldatinfo)
                else:
                    # Might have to wait since bsx_type is only blocking
                    # datataking mode...
                    # TODO: handle case when both ACC and BFS are running
                    if self.bfs and not self.acc:
                        if not bfs_yield:
                            # BFS does not need to block and
                            # so with no BSX have to throttle it...
                            _LOGGER.info("waiting {}s for BFS".format(bfs_wait))
                            time.sleep(bfs_wait)
                            if bfs_wait==file_dur:
                                bfs_wait = 1
                    if self.acc:
                        if not acc_yield:
                            # ACC does block and with no BSX have to throttle it
                            if normal_acc_wait:
                                wait_acc = modeparms.ACC_DUR
                            else:
                                wait_acc = 10
                            _LOGGER.info("waiting {}s for ACC".format(wait_acc))
                            time.sleep(wait_acc)
                            normal_acc_wait = False
                        else:
                            normal_acc_wait = True
                if self.bfs:
                    try:
                        bfs_yield = bfs_subscan.send(continue_scan)
                    except StopIteration:
                        continue_scan = False
                    if bfs_yield:
                        if self.ldatinfos_bfs == []:
                            # 1st yield is scanrecpath
                            self.scanresult['bfs'].scanrecpath = bfs_yield
                            # subsequent not None yields are ldatinfos
                            bfs_yield = bfs_subscan.send(continue_scan)
                        self.ldatinfos_bfs.append(bfs_yield)
                if self.acc:
                    acc_yield = acc_subscan.send(continue_scan)
                    if acc_yield:
                        if self.ldatinfos_acc == []:
                            # 1st yield is scanrecpath
                            self.scanresult['acc'].scanrecpath = acc_yield
                            # subsequent not None yields are ldatinfos
                            acc_yield = acc_subscan.send(continue_scan)
                        self.ldatinfos_acc.append(acc_yield)
                continue_top_scan = yield continue_scan
                continue_scan = continue_scan and continue_top_scan

        elif rec_type == 'tbb':
            stndrv.do_tbb(duration_tot, freqsetup.rcubands[0])
        elif rec_type == 'dmp':
            stndrv.halt_observingstate_when_finished = False
            stndrv.exit_check = False
            rectime = starttime
            lanes = modeparms.getlanes(freqsetup.subbands_spw,
                                       freqsetup.bits,
                                       freqsetup.nrlanes)
            band = freqsetup.rcubands[0]
            scanpath_bfdat = stndrv.bf_data_dir
            stnid = stndrv.get_stnid()
            ports = (stndrv.bf_port0, stndrv.bf_port0+1,
                     stndrv.bf_port0+2, stndrv.bf_port0+3)
            _datafiles, _logfiles = stndrv._dru_interface._rec_bf_proxy(
                ports, duration_tot, scanpath_bfdat, starttime=rectime,
                band=band, stnid=stnid)

    def close(self):
        rec_type = self.rec_type
        stndrv = self.stndrv
        acc = self.acc
        bfs =  self.bfs
        pointing_spec = self.pointing_spec
        ldatinfos = self.ldatinfos
        ldatinfo_bfs = self.ldatinfos_bfs

        bsx_type = modeparms._xtract_bsx(rec_type)

        if rec_type != 'tbb' and rec_type != 'dmp':
            if self.dir_bmctl:
                stndrv.stop_beam()
                del(self.dir_bmctl)
            elif stndrv.septonconf:
                # No pointing and tiles-off mode, so stop tiles-off mode
                stndrv.stop_tof()
            if acc:
                self.scanresult['acc'].sourcename = pointing_spec['source']
                self.scanresult['acc']._pointing = pointing_spec['pointing']
                self._stop_acc_scan()
            if bsx_type:
                self.scanresult['bsx'].sourcename = pointing_spec['source']
                self.scanresult['bsx']._pointing = pointing_spec['pointing']
                self._stop_bsx_scan(ldatinfos)
            if bfs:
                self.scanresult['bfs'].sourcename = pointing_spec['source']
                self.scanresult['bfs']._pointing = pointing_spec['pointing']
                self._stop_bfs_scan(ldatinfo_bfs)
        elif rec_type == 'tbb':
            pass
        elif rec_type == 'dmp':
            pass

        self._write_scanrecs()

    def _stop_acc_scan(self):
        """\
        Stop ACC scan recording
        """
        self.stndrv.stop_acc_scan()
        # Check if there are left over ACCs on LCU that were missed
        # now that acc scan stopped.
        leftacc_filtims = self.stndrv.get_datafiletimes(acc=True)
        flush_acc = False
        if len(leftacc_filtims):
            flush_acc = True
        if flush_acc:
            # Move left-over ACCs to DRU
            _LOGGER.info("Flushing leftover ACCs")
            _scanrecpath = self.scanresult['acc'].scanrecpath
            self.stndrv.flush_acc(_scanrecpath)
            # Add their ldatinfos and write them to file
            for filtim in leftacc_filtims:
                self.ldatinfos_acc.append(self.ldatinfos_acc[0])
                self.ldatinfos_acc[-1].filenametime = filtim
                self.ldatinfos_acc[-1].write_ldat_header(_scanrecpath)
        # Create obsinfo each ACC file
        for ldatinfo in self.ldatinfos_acc:
            self.scanresult['acc'].add_obs(ldatinfo)

    def _stop_bsx_scan(self, ldatinfos):
        """\
        Stop BSX scan recording
        """
        for ldatinfo in ldatinfos:
            self.scanresult['bsx'].add_obs(ldatinfo)

    def _stop_bfs_scan(self, ldatinfos_bfs):
        """\
        Stop BFS scan
        """
        def _DRU2CCUpath(drupath):
            ccu2dru = self.stndrv.get_ccu2dru()
            # Convert drupath from abs path on DRU to rel path from root:
            drupath = os.path.relpath(drupath, os.sep)
            # Determine path to bfs datafiles on DRU via CCU:
            ccupath = os.path.join(ccu2dru, drupath)
            return ccupath

        bfsdatapaths, bfslogpaths = self.stndrv.get_bfsdatlogpaths()
        for ldatinfo in ldatinfos_bfs:
            self.scanresult['bfs'].add_obs(ldatinfo)

        # Get path under Scans folder that will contain BFS data
        scanrecpath = self.scanresult['bfs'].scanrecpath
        #if self.stndrv._dru_interface.hostname == 'localhost':
        # Make hard links between actual BFS files and scanrec folder:
        for abfsdatapath in bfsdatapaths:
            _drusrc = _DRU2CCUpath(abfsdatapath)
            # Create new link path under scanrecpath on DRU via CCU:
            _basename = os.path.basename(abfsdatapath)
            _lnkname = os.path.join(scanrecpath, _basename)
            try:
                os.link(_drusrc, _lnkname)
            except:
                os.symlink(_drusrc, _lnkname)
        # Move logs to scanrec folder
        for abfslogpath in bfslogpaths:
            _drusrc = _DRU2CCUpath(abfslogpath)
            shutil.move(_drusrc, scanrecpath)

    def _gen_scan_id(self, scan_started):
        """\
        Generate Scan ID
        """
        scan_dt = modeparms.timestr2datetime(scan_started)
        scan_mjd_id = modeparms.dt2mjd(scan_dt)
        scan_id = str(scan_mjd_id)
        return scan_id

    def _write_scanrecs(self):
        """\
        Write the scanrec for each recorded ldat
        """
        for ldat in self.scanresult['rec']:
            scanrecinfo = self.scanresult.get(ldat)
            scanrecpath = scanrecinfo.scanrecpath
            if scanrecpath:
                scanrecinfo.write_scanrec(scanrecpath)

    def describe_scan(self):
        """\
        User friendly description of the scan
        """
        rec = self.rec_type
        id = self.scanresult['scan_id']
        pnting = self.pointing_spec.get('source')
        if not pnting:
            pnting = self.pointing_spec.get('pointing')
            if not pnting:
                pnting = self.pointing_spec.get('direction')
        freqarg = self.freqsetup.freq_arg
        descriptor = "id: {}, rec: {}, pnting: {}, freq: {}"\
                     .format(id, rec, pnting, freqarg)
        return descriptor


def still_time_fun(stoptime):
    """\
    Construct a function that says if there is still time before 'stoptime'

    Parameters
    ----------
    stoptime : datetime
        The UT datetime after which this function will be False.

    Returns
    -------
    passed_time : function
        A function that is True if current UT time is before the stoptime,
        otherwise False.
    """
    def still_time():
        now = datetime.datetime.utcnow()
        timeleft = stoptime - now
        secondsleft = int(timeleft.total_seconds())
        _LOGGER.info('Stop condition: TIME LEFT {}'.format(secondsleft))
        if secondsleft > 0:
            return True
        else:
            _LOGGER.info("Stop condition: TIME'S UP")
            return False
    return still_time


def do_nominal_scan(args):
    """\
    Do a nominal scan on a station

    Parameters
    ----------
    args : argparse.Namespace
        Arguments object. See main_cli().
    """
    accessconf = ilisa.operations.default_access_lclstn_conf()
    stndrv = StationDriver(accessconf['LCU'], accessconf['DRU'],
                           mockrun=args.mockrun)
    stndrv._lcu_interface._fake_slow_conn = 0  # Fake slow connection for test
    sesspath = stndrv.dru_data_root  # accessconf['DRU']['LOFARdataArchive']
    rec_type = args.ldat_type
    if rec_type == 'None':
        rec_type = None
    if args.acc:
        sesspath = os.path.join(sesspath, 'acc')
    elif args.bfs:
        sesspath = os.path.join(sesspath, 'bfs')
    else:
        sesspath = os.path.join(sesspath, args.ldat_type)
    freqsetup = modeparms.FreqSetup(args.freqspec)
    # Postprocess beam to get direction
    pointing = args.pointing
    direction = directions.normalizebeamctldir(pointing)
    source = 'None'
    if not directions.pointing_str2tuple(pointing):
        source = pointing
    try:
        duration_tot = float(eval(str(args.duration_tot)))
    except NameError:
        _LOGGER.error("Cannot understand duration='{}'.".format(args.duration_tot))
        raise
    # Start criteria: Time

    starttime = modeparms.as_asapdatetime(args.timestart)
    _nowtime = waituntil(starttime, datetime.timedelta(seconds=2))
    stoptime = (_nowtime + datetime.timedelta(seconds=int(duration_tot))
                + datetime.timedelta(seconds=0))
    _LOGGER.info('Expected scan stop @ {}'.format(stoptime))
    stop_cond = lambda : True  # Always True, thus wait for LCU procs to finish
    # Initialize LScan
    _pointing_spec = {'pointing': pointing, 'direction': direction,
                      'source': source}
    lscan = LScan(stndrv, rec_type, freqsetup, _pointing_spec, duration_tot,
                  args.integration, starttime, args.acc, args.bfs,
                  destpath=sesspath, destpath_bfs='Scans')
    subscan = iter(lscan)
    # Start the subscan
    next(subscan)
    while stop_cond():
        try:
            _LOGGER.debug('do_nominal_scan: IN SUBSCAN LOOP')
            _ = subscan.send(stop_cond())
        except StopIteration:
            break
    _LOGGER.debug('do_nominal_scan: END SUBSCAN LOOP')
    lscan.close()
    for res in lscan.scanresult['rec']:
        _LOGGER.info("Saved {} scanrec here: {}"
                     .format(res, lscan.scanresult[res].scanrecpath))
    if not lscan.scanresult['rec']:
        _LOGGER.info("No data recorded ('None' selected)")


def main_cli():
    """
    Record LOFAR station data via CLI

    Entry_point for ilisa_rec.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mockrun', help="Run mock rec",
                        action='store_true')
    parser.add_argument('-t', '--timestart',
                        help="Start time: YYYY-mm-ddTHH:MM:SS [UT] or 'ASAP'",
                        type=str, default='ASAP')
    parser.add_argument('-i', '--integration',
                        help="Integration time [s]",
                        type=float, default=modeparms.MIN_STATS_INTG)
    parser.add_argument('-a', '--acc', help="Enable ACC",
                        action='store_true')
    parser.add_argument('-b', '--bfs', help="Record also BFS",
                        action='store_true')
    parser.add_argument('ldat_type',
                        help="lofar data type to record."
                             " Choose from 'bst', 'sst', 'tbb', 'xst', 'dmp'"
                             " or 'None'.")
    parser.add_argument('freqspec',
                        help='Frequency spec in Hz.')
    parser.add_argument('duration_tot',
                        help='Duration in seconds. '
                             '(Can be an arithmetic expression)',
                        type=str)
    parser.add_argument('pointing', nargs='?', default=None,
                        help='Direction in az,el,ref (radians) or source name.')
    args = parser.parse_args()

    do_nominal_scan(args)

if __name__ == "__main__":
    main_cli()
