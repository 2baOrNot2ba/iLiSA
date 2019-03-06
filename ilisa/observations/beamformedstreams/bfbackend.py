"""Backends for beamformed data streams."""
import sys
import os
import subprocess

import ilisa.observations.modeparms
import ilisa.observations.beamformedstreams


dumpername = 'dump_udp_ow_4'
pathtodumper = os.path.dirname(ilisa.observations.beamformedstreams.__file__)
dumpercmd = os.path.join(pathtodumper, dumpername)
# dumpercmd = 'echo'  # For testing purposes


def _startlanerec(lane, starttimestr, duration, band, bf_data_dir, port0, stnid):
    """Start recording a lane using an external dumper process.
    """
    port = port0 + lane
    pre_bf_dir, pst_bf_dir = bf_data_dir.split('?')
    outdumpdir = pre_bf_dir + str(lane) + pst_bf_dir
    outfilepre = "udp_" + stnid
    rcumode = ilisa.observations.modeparms.band2rcumode(band)
    dumplogname = outdumpdir + dumpername \
                  + '_lane' + str(lane) \
                  + '_rcu' + rcumode \
                  + '.log'
    subprocess.call(dumpercmd + ' --ports ' + str(port) + ' --check '
                    + ' --Start ' + starttimestr
                    + ' --duration ' + str(duration)
                    + ' --timeout 9999'
                    + ' --out ' + outdumpdir + outfilepre
                    + ' > ' + dumplogname,
                    shell=True
                    )
    print("dumper log written in {}".format(dumplogname))
    datafileguess = outdumpdir + outfilepre + '_' + str(
        port) + '.start.' + starttimestr + '.000'
    print("Hopefully created dumpfile: ".format(datafileguess))


def rec_bf_streams(starttimestr, duration, lanes, band, bf_data_dir, port0, stnid):
    """Basically a wrapper that runs dump_udp processes to capture beamformed data streams.
    It sets up multiple forked processes that record one lane each."""
    child_pids = []
    child_lanes = []
    for lane in lanes:
        newpid = os.fork()
        if newpid == 0:
            _startlanerec(lane, starttimestr, duration, band, bf_data_dir, port0, stnid)
            sys.exit(0)
        else:
            child_pids.append(newpid)
    for lane in lanes:
        pid, status = os.waitpid(child_pids[lane], 0)
        print("lane {} finished with status {}.".format(lane, status))
