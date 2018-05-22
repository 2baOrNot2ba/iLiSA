#!/usr/bin/python
"""A script to record beamformed stream data."""
#TobiaC 2018-04-26 (2018-04-26)

import sys
import os
import subprocess
import time
import datetime
import argparse
import ilisa.observations.beamformedstreams
import ilisa.observations.observing as observing
import ilisa.observations.stationcontrol as stationcontrol

BLOCK = False
dumpername = 'dump_udp_ow_4'
pathtodumper = os.path.dirname(ilisa.observations.beamformedstreams.__file__)
dumpercmd = os.path.join(pathtodumper, dumpername)
#dumpercmd = 'echo' #For testing purposes


def startlanerec(lane, starttimestr, duration, band):
    """Start recording a lane.
    """
    # myObsSes refers to object in global namespace
    bf_data_dir, port0, stnid = myObsSes.bf_data_dir, myObsSes.bf_port0, \
                                myObsSes.stationcontroller.stnid
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


def setuprecording(starttimestr, duration, lanes, band):
    """Setup multiple forked processes to record lanes."""
    child_pids = []
    child_lanes = []
    for lane in lanes:
        newpid = os.fork()
        if newpid==0:
            startlanerec(lane, starttimestr, duration, band)
            sys.exit(0)
        else:
            child_pids.append(newpid)
    for lane in lanes:
        pid, status = os.waitpid(child_pids[lane], 0)
        print("lane {} finished with status {}.".format(lane,status))


if __name__=="__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("starttimestr",
                 help='Start-time (format: YYYY-mm-ddTHH:MM:SS)')
    parser.add_argument("duration",
                        help='Duration of observation in seconds')
    parser.add_argument("band",
                        help='Band: 30_90, 110_190, 210_250')
    parser.add_argument("pointing",
                        help='Pointing direction: eg CasA')
    parser.add_argument("--halt", action="store_true",
                        help="Halt observing state on completion")
    args = parser.parse_args()
    
    starttimestr = args.starttimestr
    duration = eval(args.duration) # Duration in seconds
    band = args.band # RCUMODE={LBA : 10_90 3, HBAlo : 110_190 5, HBAhi : 210_250 7}
    
    ###
    bits = 8 #8
    attenuation = None
    # Subbands allocation
    if band == '10_90' or band == '30_90':
        # LBA
        lanes = (0,1) #(0,1)
        beamletIDs = '0:243' #'0:243'
        subbandNrs = '164:407' #'164:407'
    elif band == '110_190':
        # HBAlo
        lanes = (0,1,2,3) # Normally (0,1,2,3) for all 4 lanes.
        beamletIDs = '0:487'
        subbandNrs = '12:499'
    elif band == '210_250':
        # HBAhi
        lanes = (0,1)
        beamletIDs = '0:243'
        subbandNrs = '12:255'
    else :
        raise ValueError, \
        "Wrong band: should be 10_90 (LBA), 110_190 (HBAlo) or 210_250 (HBAhi)."
    pointing = observing.stdPointings(args.pointing)

    # Wait until it is time to start
    nw=datetime.datetime.utcnow()
    st=datetime.datetime.strptime(starttimestr,"%Y-%m-%dT%H:%M:%S")

    pause = 5 # Sufficient?
    maxminsetuptime=datetime.timedelta(seconds=105+pause) # Longest minimal time
                                                          # before observation
                                                          # start to set up
    d= (st-maxminsetuptime)-nw
    timeuntilboot = d.total_seconds()
    if timeuntilboot < 0.:
        timeuntilboot = 0
    print("Will boot to observe state after "+str(timeuntilboot)+" seconds...")
    time.sleep(timeuntilboot)
    print("Booting @ {}".format(datetime.datetime.utcnow()))
    # From swlevel 0 it takes about 1:30min? to reach swlevel 3
    myObsSes = observing.Session()
    
    # Necessary since fork creates multiple instances of myObsSes and each one
    # will call it's __del__ on completion and __del__ shutdown...
    myObsSes.halt_observingstate_when_finished = False
    myObsSes.exit_check = False
    
    # BEGIN Dummy or hot beam start: (takes about 10sec) 
    # TODO: This seems necessary, otherwise beamctl will not start up next time,
    #       although it should not have to necessary.)
    print("Running warmup beam... @ {}".format(datetime.datetime.utcnow()))
    myObsSes.stationcontroller.runbeamctl(beamletIDs, subbandNrs, band, pointing)
    myObsSes.stationcontroller.rcusetup(bits, attenuation) #setting bits also seems necessary.
    myObsSes.stationcontroller.stopBeam()
    # END Dummy or hot start

    print("Pause {}s after boot.".format(pause))
    time.sleep(pause)
    
    # Real beam start:
    print("Now running real beam... @ {}".format(datetime.datetime.utcnow()))
    beamctl_CMD = myObsSes.stationcontroller.runbeamctl(beamletIDs, subbandNrs,
                                                                 band, pointing)
    rcu_setup_CMD = myObsSes.stationcontroller.rcusetup(bits, attenuation)
    nw=datetime.datetime.utcnow()
    timeleft = st-nw
    if timeleft.total_seconds() < 0.:
        starttimestr = nw.strftime("%Y-%m-%dT%H:%M:%S")
    print("(Beam started) Time left before recording: {}".format(
                                                     timeleft.total_seconds()))
    
    REC = True
    if REC == True:
        setuprecording(starttimestr, duration, lanes, band)
    else:
        print("Not recording")
    sys.stdout.flush()
    myObsSes.stationcontroller.stopBeam()
    headertime = datetime.datetime.strptime(starttimestr,"%Y-%m-%dT%H:%M:%S"
                                                     ).strftime("%Y%m%d_%H%M%S")
    myObsSes.create_LOFARst_header('bf', '.', headertime, "", beamctl_CMD,
                                                              rcu_setup_CMD, "")
    myObsSes.halt_observingstate_when_finished = args.halt # This modifies
                                                           # Session destruction

