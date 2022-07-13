"""Backends for beamformed data streams."""
import sys
import os
import subprocess
import multiprocessing
import platform
import argparse
import datetime

from . import __version__ as pipeline_version
import ilisa.pipelines.rec_bf_streams_py as rec_bf_streams_py
DATETIMESTRFMT = '%Y-%m-%dT%H:%M:%S'  # Same as in ilisa.operations

def timestr2datetime(timestr):
    # Note: this is the same as in modeparms, but to avoid pipeline package's
    # dependence on anything I copied it here.
    """\
    Convert time string into a python datetime object

    Parameters
    ----------
    timestr: str
        Date-Time string in ISO-like format '%Y-%m-%dT%H:%M:%S'
        OR 'NOW' or 'ASAP', which imply the current UT datetime.

    Returns
    -------
    dattim: datetime.datetime
        Python datetime object corresponding to input.
    """
    if timestr == 'NOW' or timestr == 'ASAP':
        # Set time to nearest rounded second from now:
        dattim = datetime.datetime.utcnow()
        dattim = dattim.replace(microsecond=0)
        dattim += datetime.timedelta(seconds=1)
    else:
        try:
            dattim = datetime.datetime.strptime(timestr, DATETIMESTRFMT)
        except:
            raise RuntimeError("Wrong datetime format.")
    return dattim

# DUMPERNAME is name of binary executable on DRU which is run by
# PL_REC_WRAPPER when capturing UDP packets with LOFAR beamformed voltages data.
# (Currently requires manually install putting it in DRU user's PATH env var)
DUMPERNAME = 'dump_udp_ow'  # Alias to local version
#pathtodumper = os.path.dirname(ilisa.pipelines.__file__)
DUMPERCMD = DUMPERNAME  # Assume dumper is in user's PATH
# DUMPERCMD = 'echo'  # For testing purposes
PL_REC_WRAPPER = 'pl_rec'


def bfsfilepaths(lane, starttime, rcumode, bf_data_dir, port0, stnid,
                 compress=True):
    """Generate paths and name for BFS recording.

    Parameters
    ----------
    lane : int
        Lane number 0,1,2, or 3.
    starttime : str
        The datetime string when the BF stream started.
    rcumode :
        The band name for the BF stream.
    bf_data_dir : str
        Template for BF data lane dump directory. Should have format:
            <pre_bf_dir>?<pst_bf_dir>
        where '?' will be replaced by the lane number.
    port0 :
        The port number of lane 0.
    stnid :
        Station ID.
    compress: bool
        Whether or not compression is used.

    Returns
    -------
    outdumpdir : str
        Directory in which to dump bf data. Has format:
            <outdumpdir>/udp_<stnid>
        where
            <outdumpdir> := <rootdir>/lane<lanenr>/<pst_bf_dir>
    outarg : str
        Argument 'out' passed to dumper CLI.
    datafileguess : str
        Path to data file. Has format:
            _<port>.start.<%Y-%m-%dT%H:%M:%S>.000
    dumplogname :
        Name of dumper's logfile.
    """
    port = port0 + lane
    pre_bf_dir, pst_bf_dir = bf_data_dir.split('?')
    outdumpdir = pre_bf_dir + str(lane) + pst_bf_dir
    outfilepre = "udp_" + stnid
    outarg = os.path.join(outdumpdir, outfilepre)
    dumplogname = os.path.join(outdumpdir,
                               '{}_lane{}_rcu{}.log'.format(DUMPERNAME, lane,
                                                            rcumode))
    datafileguess = outarg + '_' + str(port) + platform.node()\
        + starttime.strftime("%Y-%m-%dT%H:%M:%S") + '.000'
    if compress:
        datafileguess += '.zst'
    return outdumpdir, outarg, datafileguess, dumplogname


def _startlanerec(lane, starttime, duration, file_dur, rcumode, bf_data_dir,
                  port0, stnid, compress=True, threadqueue=None):
    """Start recording a lane using an external dumper process.
    """
    if compress:
        # Compress stored data using zstd.
        compress_flag = ' --compress'
    else:
        compress_flag = ''
    port = port0 + lane
    outdumpdir, outarg, datafileguess, dumplogname = \
        bfsfilepaths(lane, starttime, rcumode, bf_data_dir,
                     port0, stnid, compress)
    if not os.path.exists(outdumpdir):
        os.mkdir(outdumpdir)
    dur_flagarg = ''
    if duration:
        dur_flagarg = ' --duration ' + str(int(duration))
    maxfilesz_arg = ''
    if file_dur:
        datarate_perlane = 95.5E6  # Bytes per second per lane approx
        maxfilesz = int(file_dur * datarate_perlane)
        maxfilesz_arg = ' --Maxfilesize ' + str(maxfilesz)
    cmdline_full_shell = (DUMPERCMD + ' --ports ' + str(port) + ' --check '
                          + ' --Start ' +
                          starttime.strftime("%Y-%m-%dT%H:%M:%S")
                          + dur_flagarg
                          + maxfilesz_arg
                          + ' --timeout 9999'
                          + compress_flag
                          + ' --out ' + outarg
                          + ' > ' + dumplogname)
    print("Running: {}".format(cmdline_full_shell))
    subprocess.run(cmdline_full_shell, shell=True)
    print("{}".format(datafileguess))
    print("{}".format(dumplogname))
    if threadqueue:
        threadqueue.put((datafileguess, dumplogname))
    # return datafileguess, dumplogname


def rec_bfs_lanes(starttime, duration, file_duration, lanes, rcumode,
                  bf_data_dir, port0, stnid, compress):
    retvalq = multiprocessing.Queue()
    laneprocs = []
    for lane in lanes:
        laneproc = multiprocessing.Process(target=_startlanerec,
                                           args=(lane, starttime, duration,
                                                 file_duration, rcumode,
                                                 bf_data_dir, port0, stnid,
                                                 compress, retvalq))
        laneproc.start()
        laneprocs.append(laneproc)
    datafiles = []
    logfiles = []
    for laneproc in laneprocs:
        laneproc.join()
        datafileguess, dumplogname = retvalq.get()
        datafiles.append(datafileguess)
        logfiles.append(dumplogname)
    return datafiles, logfiles


def rec_bf_streams(starttime, duration, file_duration, lanes, rcumode,
                   bf_data_dir, port0, stnid, compress):
    """
    Wrapper that runs a beamformed data stream capture command for each lane

    Uses either fork or a multiprocess Process.
    Invokes either udp_dump_ow or iLiSA python recorder.
    """
    usefork = False
    use_python_recorder = False
    if usefork:
        child_pids = []
        _child_lanes = []
        for lane in lanes:
            newpid = os.fork()
            if newpid == 0:
                _startlanerec(lane, starttime, duration, rcumode, bf_data_dir,
                              port0, stnid)
                sys.exit(0)
            else:
                child_pids.append(newpid)
        for lane in lanes:
            _pid, status = os.waitpid(child_pids[lane], 0)
            print("lane {} finished with status {}.".format(lane, status))
        datafiles, logfiles = [None]*len(lanes), [None]*len(lanes)
    else:
        if not use_python_recorder:
            retvalq = multiprocessing.Queue()
            datafiles = []
            logfiles = []
            laneprocs = []
            for lane in lanes:
                laneproc = multiprocessing.Process(target=_startlanerec,
                                                   args=(lane, starttime,
                                                         duration,
                                                         file_duration,
                                                         rcumode,
                                                         bf_data_dir, port0,
                                                         stnid, compress,
                                                         retvalq))
                laneproc.start()
                laneprocs.append(laneproc)
            for laneproc in laneprocs:
                laneproc.join()
                datafileguess, dumplogname = retvalq.get()
                datafiles.append(datafileguess)
                logfiles.append(dumplogname)
        else:
            datafiles, logfiles = rec_bf_streams_py.main(port0, bf_data_dir,
                                                         duration=duration)
    return datafiles, logfiles


def bfsrec_main_cli():
    # Entry point for pl_rec
    parser = argparse.ArgumentParser()
    parser.add_argument('-m', '--mockrun', help="Run mock rec",
                        action='store_true')
    parser.add_argument('-t', '--starttime',
                        type=str, default='NOW',
                        help = "Start-time: (iso format) YYYY-mm-ddTHH:MM:SS"
                        )
    parser.add_argument('-p', '--ports',
                        type=str, default='4346',
                        help="List of port number(s)"
                        )
    parser.add_argument('-b', '--bfdatadir',
                        type=str, default='/mnt/lane?/',
                        help="Template directory for BF data"
                        )
    parser.add_argument('-d', '--duration',
                        type=float, default=None,
                        help="Duration of recording in seconds"
                        )
    parser.add_argument('-f', '--file_duration',
                        type=float, default=None,
                        help="Duration of dumped files in seconds"
                        )
    parser.add_argument('-w', '--which',
                        type=str, default='ow',
                        help="Which backend recorder: ow or py",
                        )
    parser.add_argument('-r', '--rcumode',
                        type=str, default='5',
                        help="rcumode or spectral window",
                        )
    parser.add_argument('-s', '--stnid',
                        type=str, default='SE607',
                        help="station id",
                        )
    parser.add_argument('-c', '--compress', action="store_true",
                        help="Compress recorded data")
    parser.add_argument('-v', '--version', action="store_true",
                        help="Print version of module")
    args = parser.parse_args()
    starttime = timestr2datetime(args.starttime)
    args.ports = [int(portstr) for portstr in args.ports.split(',')]
    if args.version:
        print(pipeline_version)
        return
    if not args.mockrun:
        if args.which == 'py':
            rec_bf_streams_py.main(args.ports[0], args.bfdatadir, args.duration)
        else:
            lanes = range(len(args.ports))
            rec_bfs_lanes(starttime, args.duration, args.file_duration, lanes,
                          args.rcumode, args.bfdatadir, args.ports[0],
                          args.stnid, args.compress)
    else:
        print('MOCKRUN. Arguments:')
        print(args)


if __name__ == '__main__':
    bfsrec_main_cli()
