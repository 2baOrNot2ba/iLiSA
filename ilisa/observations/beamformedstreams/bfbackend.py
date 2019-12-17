"""Backends for beamformed data streams."""
import sys
import os
import subprocess

import ilisa.observations.modeparms
import ilisa.observations.beamformedstreams


dumpername = 'dump_udp_ow_12'
pathtodumper = os.path.dirname(ilisa.observations.beamformedstreams.__file__)
dumpercmd = os.path.join(pathtodumper, dumpername)
# dumpercmd = 'echo'  # For testing purposes

def bfsfilepaths(lane, starttime, band, bf_data_dir, port0, stnid):
    """Generate paths and name for BFS recording."""
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
    return outdumpdir, outarg, datafileguess, dumplogname

def _startlanerec(lane, starttime, duration, band, bf_data_dir, port0, stnid,
                  compress=True):
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
                                                                  port0, stnid)
    if not os.path.exists(outdumpdir):
        os.mkdir(outdumpdir)
    subprocess.call(dumpercmd + ' --ports ' + str(port) + ' --check '
                    + ' --Start ' + starttime.strftime("%Y-%m-%dT%H:%M:%S")
                    + ' --duration ' + str(duration)
                    + ' --timeout 9999'
                    + compress_flag
                    + ' --out ' + outarg
                    + ' > ' + dumplogname,
                    shell=True
                    )
    print("dumper log written in {}".format(dumplogname))
    print("Hopefully created dumpfile: ".format(datafileguess))


def rec_bf_streams(starttime, duration, lanes, band, bf_data_dir, port0, stnid):
    """Basically a wrapper that runs dump_udp processes to capture beamformed data streams.
    It sets up multiple forked processes that record one lane each."""
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
