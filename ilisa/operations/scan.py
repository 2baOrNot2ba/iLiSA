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


class LScan:
    def __init__(self, stndrv, rec_type, freqsetup, duration_tot,
                 pointing_spec, integration, starttime, acc=False, bfs=False,
                 destpath=None, destpath_bfs=None, file_dur=None, scan_id=None):
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
        duration_tot: float
            Requested total duration of scan in seconds.
        pointing_spec: dict
            Pointing source direction of scan. Can be a source name or a
            beamctl direction str. If None, then allsky image is implied.
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

        pointing = pointing_spec['pointing']
        self.dir_bmctl = ilisa.operations.directions.normalizebeamctldir(pointing)
        if pointing and not self.dir_bmctl:
            raise ValueError("Invalid pointing syntax: {}".format(pointing))
        beam_needed = bfs or rec_type == 'bst'
        if beam_needed and not pointing:
            raise ValueError("No pointing, but beam needed")

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
                self.scanresult['acc'] = data_io.ScanRecInfo()
                self.scanresult['acc'].set_stnid(self.stndrv.get_stnid())
                self.scanresult['acc'].set_caltabinfos([])
                self.scanresult['acc'].set_scanrecparms(
                    'acc', freqsetup.arg, duration_tot,
                    self.pointing_spec['direction'], 1.0)
            if bsx_stat:
                self.scanresult['rec'].append('bsx')
                self.scanresult['bsx'] = data_io.ScanRecInfo()
                self.scanresult['bsx'].set_stnid(self.stndrv.get_stnid())
                self.scanresult['bsx'].set_caltabinfos(caltabinfos)
                self.scanresult['bsx'].set_scanrecparms(
                    bsx_stat, freqsetup.arg, duration_tot,
                    self.pointing_spec['direction'], integration)
            if bfs:
                self.scanresult['rec'].append('bfs')
                self.scanresult['bfs'] = data_io.ScanRecInfo()
                self.scanresult['bfs'].set_stnid(self.stndrv.get_stnid())
                self.scanresult['bfs'].set_caltabinfos(caltabinfos)
                self.scanresult['bfs'].set_scanrecparms(
                    'bfs', freqsetup.arg, duration_tot,
                    self.pointing_spec['direction'], None)

        self.ldatinfos = []
        self.ldatinfos_acc = []
        self.ldatinfos_bfs = []

    def __iter__(self):
        stndrv = self.stndrv
        rec_type = self.rec_type
        freqsetup = self.freqsetup
        duration_tot = self.duration_tot
        pointing_spec = self.pointing_spec
        integration = self.integration
        starttime = self.starttime
        file_dur = self.file_dur

        pointing = pointing_spec['pointing']
        bsx_type = modeparms._xtract_bsx(rec_type)

        if rec_type != 'tbb' and rec_type != 'dmp':
            if pointing:
                stndrv.field = pointing_spec['source']
                if self.acc:
                    # ACC needs to be enabled before beam
                    # so initialize the subscan generator...
                    acc_subscan = stndrv.start_acc_scan()
                    # actually start the subscan...
                    _ = next(acc_subscan)
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
            acc_yield = None
            # bfs_yield = None
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
                    if self.bfs:
                        if not bfs_yield:
                            # BFS does not need to block and
                            # so with no BSX have to throttle it...
                            logging.info(f'waiting {file_dur}s for BFS')
                            time.sleep(file_dur)
                    if self.acc:
                        if not acc_yield:
                            # ACC does block and with no BSX have to throttle it
                            if normal_acc_wait:
                                wait_acc = modeparms.ACC_DUR
                            else:
                                wait_acc = 10
                            logging.info(f'waiting {wait_acc}s for ACC')
                            time.sleep(wait_acc)
                            normal_acc_wait = False
                        else:
                            normal_acc_wait = True
                if self.bfs:
                    bfs_yield = bfs_subscan.send(continue_scan)
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

        pointing = pointing_spec['pointing']
        bsx_type = modeparms._xtract_bsx(rec_type)

        if rec_type != 'tbb' and rec_type != 'dmp':
            if pointing:
                stndrv.stop_beam()
                if acc:
                    self.scanresult['acc'].sourcename = pointing_spec['source']
                    self.scanresult['acc']._pointing = pointing_spec['pointing']
                    self._stop_acc_scan()
            elif stndrv.septonconf:
                # No pointing and tiles-off mode, so stop tiles-off mode
                stndrv.stop_tof()
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
        bfsdatapaths, bfslogpaths = self.stndrv.get_bfsdatlogpaths()
        for ldatinfo in ldatinfos_bfs:
            self.scanresult['bfs'].add_obs(ldatinfo)

        # Make a project folder for BFS data
        scanrecpath = self.scanresult['bfs'].scanrecpath
        # Create BFS destination folder on DPU:
        if self.stndrv._dru_interface.hostname == 'localhost':
            # Make soft links to actual BFS files and move logs to scanrec
            # folder
            for abfsdatapath in bfsdatapaths:
                _basename = os.path.basename(abfsdatapath)
                _lnkname = os.path.join(scanrecpath, _basename)
                os.symlink(abfsdatapath, _lnkname)
            for abfslogpath in bfslogpaths:
                shutil.move(abfslogpath, scanrecpath)

    def _gen_scan_id(self, scan_dt):
        scan_mjd_id = modeparms.dt2mjd(scan_dt)
        scan_id = "scan_{}".format(scan_mjd_id)
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
        freqarg = self.freqsetup.arg
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
        logging.info('Stop condition: TIME LEFT {}'.format(secondsleft))
        if secondsleft > 0:
            return True
        else:
            logging.info("Stop condition: TIME'S UP")
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
        logging.error("Cannot understand duration='{}'.".format(args.duration_tot))
        raise
    # Start criteria: Time
    args.starttime = modeparms.timestr2datetime(args.starttime)
    starttime = waituntil(args.starttime, datetime.timedelta(seconds=2))
    stoptime = (starttime + datetime.timedelta(seconds=int(duration_tot))
                + datetime.timedelta(seconds=0))
    logging.info('Will stop @{}'.format(stoptime))
    stop_cond = still_time_fun(stoptime)
    # Initialize LScan
    lscan = LScan(stndrv, rec_type, freqsetup, duration_tot,
                  {'pointing': pointing, 'direction': direction,
                   'source': source}, args.integration, starttime,
                  args.acc, args.bfs, destpath=sesspath, destpath_bfs='Scans')
    subscan = iter(lscan)
    # Start the subscan
    next(subscan)
    while stop_cond():
        try:
            logging.debug('do_nominal_scan: IN SUBSCAN LOOP')
            _ = subscan.send(stop_cond())
        except StopIteration:
            break
    logging.debug('do_nominal_scan: END SUBSCAN LOOP')
    lscan.close()
    for res in lscan.scanresult['rec']:
        logging.info("Saved {} scanrec here: {}"
            .format(res, lscan.scanresult[res].scanrecpath))
    if not lscan.scanresult['rec']:
        logging.info("No data recorded ('None' selected)")


def main_cli():
    """
    Record LOFAR station data via CLI

    Entry_point for ilisa_rec.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mockrun', help="Run mock rec",
                        action='store_true')
    parser.add_argument('-s', '--starttime',
                        help="Start Time (format: YYYY-mm-ddTHH:MM:SS)",
                        type=str, default='ASAP')
    parser.add_argument('-i', '--integration',
                        help="Integration time [s]",
                        type=float, default=modeparms.MIN_STATS_INTG)
    parser.add_argument('-a', '--acc', help="Enable ACC",
                        action='store_true')
    parser.add_argument('-b', '--bfs', help="Record also BST",
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
