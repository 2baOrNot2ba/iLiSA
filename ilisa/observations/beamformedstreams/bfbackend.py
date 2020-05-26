"""Backends for beamformed data streams."""
import sys
import os
import subprocess
import multiprocessing

import ilisa.observations.modeparms
import ilisa.observations.beamformedstreams
import ilisa.observations.beamformedstreams.rec_bf_streams_py as rec_bf_streams_py

dumpername = 'dump_udp_ow_12'
pathtodumper = os.path.dirname(ilisa.observations.beamformedstreams.__file__)
dumpercmd = os.path.join(pathtodumper, dumpername)
# dumpercmd = 'echo'  # For testing purposes


def bfsfilepaths(lane, starttime, band, bf_data_dir, port0, stnid, compress=True):
    """Generate paths and name for BFS recording.

    Parameters
    ----------
    lane : int
        Lane number 0,1,2, or 3.
    starttime : str
        The datetime string when the BF stream started.
    band :
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
    rcumode = ilisa.observations.modeparms.band2rcumode(band)
    outarg = os.path.join(outdumpdir, outfilepre)
    dumplogname = os.path.join(outdumpdir, '{}_lane{}_rcu{}.log'.format(dumpername,
                                                                        lane,
                                                                        rcumode))
    datafileguess = outarg + '_' + str(port) + '.start.'\
                    + starttime.strftime("%Y-%m-%dT%H:%M:%S") + '.000'
    if compress:
        datafileguess += '.zst'
    return outdumpdir, outarg, datafileguess, dumplogname


def _startlanerec(lane, starttime, duration, band, bf_data_dir, port0, stnid,
                  compress=True, threadqueue=None):
    """Start recording a lane using an external dumper process.
    """
    if compress:
        # Compress stored data using zstd.
        compress_flag = ' --compress'
    else:
        compress_flag = ''
    port = port0 + lane
    outdumpdir, outarg, datafileguess, dumplogname = bfsfilepaths(lane, starttime,
                                                                  band, bf_data_dir,
                                                                  port0, stnid, compress)
    if not os.path.exists(outdumpdir):
        os.mkdir(outdumpdir)
    cmdline_full_shell = (dumpercmd + ' --ports ' + str(port) + ' --check '
                         + ' --Start ' + starttime.strftime("%Y-%m-%dT%H:%M:%S")
                         + ' --duration ' + str(duration)
                         + ' --timeout 9999'
                         + compress_flag
                         + ' --out ' + outarg
                         + ' > ' + dumplogname)
    print("Running: {}".format(cmdline_full_shell))
    subprocess.call(cmdline_full_shell, shell=True)
    print("{}".format(datafileguess))
    print("{}".format(dumplogname))
    if threadqueue:
        threadqueue.put((datafileguess, dumplogname))
    #return datafileguess, dumplogname


def rec_bf_streams(starttime, duration, lanes, band, bf_data_dir, port0, stnid):
    """Basically a wrapper that runs dump_udp processes to capture beamformed data streams.
    It sets up multiple forked processes that record one lane each."""
    usefork = False
    use_python_recorder = False
    if usefork:
        child_pids = []
        child_lanes = []
        for lane in lanes:
            newpid = os.fork()
            if newpid == 0:
                _startlanerec(lane, starttime, duration, band, bf_data_dir, port0, stnid)
                sys.exit(0)
            else:
                child_pids.append(newpid)
        for lane in lanes:
            pid, status = os.waitpid(child_pids[lane], 0)
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
                                                   args=(lane, starttime, duration, band,
                                                         bf_data_dir, port0, stnid, True,
                                                         retvalq)
                                                   )
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
