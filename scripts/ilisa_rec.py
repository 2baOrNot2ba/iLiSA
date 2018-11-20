#!/usr/bin/python
"""Start observation on station.
Types of observation are:
* ACC
* BST
* SST
* XST
* TBB

"""

import sys
import os
import math
import subprocess
import time
import datetime
import argparse
import ilisa.observations.beamformedstreams
import ilisa.observations.observing as observing
import ilisa.observations.stationcontrol as stationcontrol
import ilisa.observations.dataIO as dataIO


dumpername = 'dump_udp_ow_4'
pathtodumper = os.path.dirname(ilisa.observations.beamformedstreams.__file__)
dumpercmd = os.path.join(pathtodumper, dumpername)
#dumpercmd = 'echo'  # For testing purposes


def startlanerec(lane, starttimestr, duration, band, bf_data_dir, port0, stnid):
    """Start recording a lane.
    """
    port = port0+lane
    pre_bf_dir, pst_bf_dir = bf_data_dir.split('?')
    outdumpdir = pre_bf_dir+str(lane)+pst_bf_dir
    outfilepre = "udp_"+stnid
    rcumode = stationcontrol.band2rcumode(band)
    dumplogname=outdumpdir+dumpername \
                          +'_lane'+str(lane) \
                          +'_rcu'+rcumode \
                          +'.log'
    subprocess.call(dumpercmd+' --ports '+str(port)+' --check '
                          +' --Start '+starttimestr
                          +' --duration '+str(duration)
                          +' --timeout 9999'
                          +' --out '+outdumpdir+outfilepre
                          +' > '+dumplogname,
                         shell=True
                   )
    print("dumper log written in {}".format(dumplogname))
    datafileguess=outdumpdir+outfilepre+'_'+str(port)+'.start.'+starttimestr+'.000'
    print("Hopefully created dumpfile: ".format(datafileguess))
#    return datafileguess


def setuprecording(starttimestr, duration, lanes, band, bf_data_dir, port0, stnid):
    """Setup multiple forked processes to record lanes."""
    child_pids = []
    child_lanes = []
    for lane in lanes:
        newpid = os.fork()
        if newpid==0:
            startlanerec(lane, starttimestr, duration, band, bf_data_dir, port0, stnid)
            sys.exit(0)
        else:
            child_pids.append(newpid)
    for lane in lanes:
        pid, status = os.waitpid(child_pids[lane], 0)
        print("lane {} finished with status {}.".format(lane,status))


def waittoboot(starttimestr, pause):
    """Before booting, wait until time given by starttimestr. This includes
     a dummy beam warmup."""
    nw = datetime.datetime.utcnow()
    st = datetime.datetime.strptime(starttimestr, "%Y-%m-%dT%H:%M:%S")

    maxminsetuptime = datetime.timedelta(seconds=105 + pause)  # Longest minimal time
    # before observation
    # start to set up
    d = (st - maxminsetuptime) - nw
    timeuntilboot = d.total_seconds()
    if timeuntilboot < 0.:
        timeuntilboot = 0
    print("Will boot to observe state after " + str(timeuntilboot) + " seconds...")
    time.sleep(timeuntilboot)
    return st


def do_bfs(args):
    """Record BeamFormed Streams (BFS)."""
    starttimestr = args.starttimestr
    duration = eval(args.duration)  # Duration in seconds
    band = args.band  # RCUMODE={LBA : 10_90 3, HBAlo : 110_190 5, HBAhi : 210_250 7}

    ###
    bits = 8  # 8
    attenuation = None
    # Subbands allocation
    if band == '10_90' or band == '30_90':
        # LBA
        lanes = (0, 1)  # (0,1)
        beamletIDs = '0:243'  # '0:243'
        subbandNrs = '164:407'  # '164:407'
    elif band == '110_190':
        # HBAlo
        lanes = (0, 1, 2, 3)  # Normally (0,1,2,3) for all 4 lanes.
        beamletIDs = '0:487'
        subbandNrs = '12:499'
    elif band == '210_250':
        # HBAhi
        lanes = (0, 1)
        beamletIDs = '0:243'
        subbandNrs = '12:255'
    else:
        raise ValueError, \
            "Wrong band: should be 10_90 (LBA), 110_190 (HBAlo) or 210_250 (HBAhi)."
    pointing = observing.normalizebeamctldir(args.pointsrc)

    # Wait until it is time to start
    pause = 5  # Sufficient?
    st = waittoboot(starttimestr, pause)

    # From swlevel 0 it takes about 1:30min? to reach swlevel 3
    print("Booting @ {}".format(datetime.datetime.utcnow()))

    # Necessary since fork creates multiple instances of myobs and each one
    # will call it's __del__ on completion and __del__ shutdown...
    myobs.halt_observingstate_when_finished = False
    myobs.exit_check = False

    # BEGIN Dummy or hot beam start: (takes about 10sec)
    # TODO: This seems necessary, otherwise beamctl will not start up next time,
    #       although it should not have to necessary.)
    print("Running warmup beam... @ {}".format(datetime.datetime.utcnow()))
    myobs.stationcontroller.runbeamctl(beamletIDs, subbandNrs, band, pointing)
    myobs.stationcontroller.rcusetup(bits,
                                        attenuation)  # setting bits also seems necessary.
    myobs.stationcontroller.stopBeam()
    # END Dummy or hot start

    print("Pause {}s after boot.".format(pause))
    time.sleep(pause)

    # Real beam start:
    print("Now running real beam... @ {}".format(datetime.datetime.utcnow()))
    beamctl_CMD = myobs.stationcontroller.runbeamctl(beamletIDs, subbandNrs,
                                                        band, pointing)
    rcu_setup_CMD = myobs.stationcontroller.rcusetup(bits, attenuation)
    nw = datetime.datetime.utcnow()
    timeleft = st - nw
    if timeleft.total_seconds() < 0.:
        starttimestr = nw.strftime("%Y-%m-%dT%H:%M:%S")
    print("(Beam started) Time left before recording: {}".format(
        timeleft.total_seconds()))

    REC = True
    if REC == True:
        bf_data_dir, port0, stnid = myobs.bf_data_dir, myobs.bf_port0, \
                                    myobs.stationcontroller.stnid
        setuprecording(starttimestr, duration, lanes, band, bf_data_dir, port0, stnid)
    else:
        print("Not recording")
    sys.stdout.flush()
    myobs.stationcontroller.stopBeam()
    headertime = datetime.datetime.strptime(starttimestr, "%Y-%m-%dT%H:%M:%S"
                                            ).strftime("%Y%m%d_%H%M%S")

    obsinfo = dataIO.ObsInfo(myobs.stationcontroller.stnid,
                             myobs.project, myobs.observer)
    obsinfo.setobsinfo_fromparams('bfs', headertime, beamctl_CMD, rcu_setup_CMD, "")
    bsxSTobsEpoch, datapath = obsinfo.getobsdatapath(myobs.LOFARdataArchive)
    print "Putting header in", datapath
    #obsinfo.create_LOFARst_header('bf', '.', headertime, "", beamctl_CMD, rcu_setup_CMD, "")
    obsinfo.create_LOFARst_header(datapath)
    myobs.halt_observingstate_when_finished = args.shutdown  # Necessary due to forking


def do_acc(args):
    """Record acc data for one of the LOFAR bands over a duration.
    """
    duration = int(eval(args.duration))
    accdestdir = myobs.do_acc(args.band, duration, args.pointsrc)
    print("Saved ACC data in folder: {}".format(accdestdir))


def _do_bsx(statistic, args):
    """Records bst,sst,xst data in one of the LOFAR bands and creates a header file
    with observational settings.
    """
    duration = int(math.ceil(eval(args.duration)))
    myobs.bits = 16
    # TODO the freqbnd is not fully functional yet. Implement it.
    if statistic == 'sst':
        if args.band == '10_90' or  args.band == '30_90':
            args.freqbnd = '50e6'
        elif args.band == '110_190' or  args.band == '170_230':
            args.freqbnd = '150e6'
        elif args.band == '210_250':
            args.freqbnd = '250e6'
    freqlo = float(args.freqbnd)
    if args.allsky and freqlo > 100.0e6 :
        myobs.do_SEPTON(statistic, freqlo, args.integration, duration)
    else:
        myobs.bsxST(statistic, freqlo, args.integration, duration,
                    args.pointsrc)


def do_bst(args):
    """Records BST data in one of the LOFAR bands and creates a header file
    with observational settings.
    """
    _do_bsx('bst', args)


def do_sst(args):
    """Records SST data in one of the LOFAR bands and creates a header file
    with observational settings.
    """
    _do_bsx('sst', args)


def do_xst(args):
    """Records XST data in one of the LOFAR bands and creates a header file
    with observational settings.
    """
    _do_bsx('xst', args)


def do_tbb(args):
    """Record Transient Buffer Board (TBB) data from one of the LOFAR bands for
    duration seconds.
    """
    myobs.do_tbb(args.duration, args.band)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--allsky', help="Set allsky FoV", action='store_true')
    parser.add_argument('-s', '--shutdown', help="Shutdown observing state when finished", action='store_true')
    subparsers = parser.add_subparsers(title='Observation mode',
                                       description='Select a type of data to record.',
                                       help='Type of datataking:')

    # Specify common parameter args:
    arg_rcuband_kwargs = {'type': str,
                       'help': """(RCU) Band to use: 10_90, 30_90, 110_190, \
                               170_230 or 210_250."""}
    arg_freqband_kwargs = {'type': str,
                           'help': """\
                            Frequency band spec in Hz.\
                            Format: freqct | freqlo:freqhi | freqlo:freqstp:freqhi\
                            """}
    arg_integration_kwargs ={'type': int,
                             'help': "Integration time in s"}
    arg_duration_kwargs = {'type': str,
                           'help': "Duration of calibration obs. in seconds. \
                                    Can be an arithmetic formula e.g. 24*60*60."}
    arg_pointsrc_kwargs = {'type': str, 'nargs': '?', 'default': 'Z',
                           'help': "Pointing [format: RA,DEC,REF with RA,DEC in radians, \
                                    REF=J2000] or a source name. Default 'Z' stands for \
                                    zenith."}

    # ACC data
    parser_acc = subparsers.add_parser('acc',
                                       help="Make an ACC observation.")
    parser_acc.set_defaults(func=do_acc)
    parser_acc.add_argument('band', **arg_rcuband_kwargs)
    parser_acc.add_argument('duration', **arg_duration_kwargs)
    parser_acc.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # BFS data
    parser_bfs = subparsers.add_parser('bfs',
                                       help="Make an BFS observation.")
    parser_bfs.set_defaults(func=do_bfs)
    parser_bfs.add_argument('starttimestr',
                 help="Start-time (format: YYYY-mm-ddTHH:MM:SS)")
    parser_bfs.add_argument('band', **arg_rcuband_kwargs)
    parser_bfs.add_argument('duration',**arg_duration_kwargs)
    parser_bfs.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # BST data
    parser_bst = subparsers.add_parser('bst',
                                       help="Make a BST observation")
    parser_bst.set_defaults(func=do_bst)
    parser_bst.add_argument('freqbnd', **arg_freqband_kwargs)
    parser_bst.add_argument('integration', **arg_integration_kwargs)
    parser_bst.add_argument('duration',**arg_duration_kwargs)
    parser_bst.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # SST data
    parser_sst = subparsers.add_parser('sst',
                                       help="Make a SST observation")
    parser_sst.set_defaults(func=do_sst)
    parser_sst.add_argument('band', **arg_rcuband_kwargs)
    parser_sst.add_argument('integration',**arg_integration_kwargs)
    parser_sst.add_argument('duration',**arg_duration_kwargs)
    parser_sst.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # XST data
    parser_xst = subparsers.add_parser('xst',
                                       help="Make a XST observation")
    parser_xst.set_defaults(func=do_xst)
    parser_xst.add_argument('freqbnd', **arg_freqband_kwargs)
    parser_xst.add_argument('integration',**arg_integration_kwargs)
    parser_xst.add_argument('duration',**arg_duration_kwargs)
    parser_xst.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # TBB data
    parser_tbb = subparsers.add_parser('tbb',
                                       help="Make a TBB observation")
    parser_tbb.set_defaults(func=do_tbb)
    parser_tbb.add_argument('band', **arg_rcuband_kwargs)
    parser_tbb.add_argument('duration', **arg_duration_kwargs)

    args = parser.parse_args()

    myobs = observing.Session()
    myobs.halt_observingstate_when_finished = args.shutdown

    args.func(args)
