#!/usr/bin/python
import sys
import multiprocessing
import subprocess
import socket
import datetime
import time
import os
import signal

integrate_step = 10
nr_lanes = 4
port = 4346
#hostIPs = ['10.211.7.2', '10.212.7.2', '10.213.7.2', '10.214.7.2']
logFileName = "/mnt/lane0/logs/BFdump.log"
datarec_dir = "/mnt/lane"
filename_base = "SE607-BFS"
filetag4dumpongoing = "ONGOING"
ext = ".lbp" #lofar beamlet packets
REC_SET = True
DATA_REDUCE = False
bufsize = 1048576


def log_obs(ObsDataFile, bc_args, OneBeamletRec):
    if not OneBeamletRec:
        ObsDataFileBase,ObsDataFileExt=ObsDataFile.split('.')
        logfile = "/mnt/lane0/BF/"+ObsDataFileBase+".log"
    else:
        # beamctl log for beamlets
        logfile = "/mnt/lane0/logs/Beamlet.log"
    lfh=open(logfile,"a")
    bcline=" ".join(bc_args.split())
    bclines="\\\n".join(bcline.split('\\'))
    if OneBeamletRec:
        lfh.write(ObsDataFile+" => "+"BFlanes "+bclines+"\n")
    else:
        lfh.write(ObsDataFile+' =>\n'+'"'+bclines+'"\n')
    lfh.close()
    print("LCL> created BF logfile: "+logfile)


def get_databurst(lane):
    try:
        clientsock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    except socket.error as msg:
        raise msg
    clientsock.bind(('', port+lane))  # clientsock.bind((hostIPs[lane], port+lane))

    # Wait until I get some data:
    print("Lane %d waiting..." % lane)
    recv_msg, addr = clientsock.recvfrom(bufsize)
    start_time = datetime.datetime.utcnow()
    print("...got some data at {}".format(start_time))
    datapath = None
    if REC_SET:
        print("Recording until data-stream times-out.")
        datapath = datarec_dir + str(lane) + "/BF/" + filename_base\
                   + filetag4dumpongoing + "_" + str(lane) + ext
        f = open(datapath, "wb")
        f.write(recv_msg)
    clientsock.settimeout(1.0)
    while True:
        try:
            recv_msg, addr = clientsock.recvfrom(bufsize)
        except socket.timeout:
            clientsock.close()
            if REC_SET:
                f.close()
            break
        if REC_SET:
            f.write(recv_msg)
    return start_time, datapath


def startlane(lane, threadqueue = None):
    start_time, tmpdatapath = get_databurst(lane)
    stop_time = datetime.datetime.utcnow()
    stop_time_near_sec = stop_time.replace(microsecond=0)

    obs_file_name = filename_base + "_" + stop_time_near_sec.isoformat() + "_"\
                    + str(lane) + ext
    datapath = datarec_dir + str(lane) + "/BF/" + obs_file_name
    # Let Mother process find ONGOING file...
    os.rename(tmpdatapath, datapath)
    print("{}: Created file: {}".format(os.path.basename(__file__), datapath))
    if threadqueue:
        threadqueue.put((lane, obs_file_name))
    if DATA_REDUCE:
        subprocess.call(
           ["/home/tobia/lofarBFdataProcess/LofBFraw2bin",
            datapath, str(integrate_step), "o"]
        )


def main(lanes, block=True):
    child_pids=[]
    for lane in lanes:
        newpid=os.fork()
        if newpid == 0:
            startlane(lane)
            exit()
        else:
            child_pids.append(newpid)
    if block:
        try :
           while True:
                time.sleep(3600)
        except KeyboardInterrupt:
            for lane in range(len(lanes)):
                os.kill(child_pids[lane],signal.SIGTERM)
                print("Stopping lane {}".format(lanes[lane]))


def main_multiproc(lanes):
    laneprocs = []
    thelanes = []
    threadqueue = multiprocessing.Queue()
    for lane in lanes:
        laneprocs.append(multiprocessing.Process(target=startlane,
                                                 args=(lane, threadqueue)))
        laneprocs[-1].start()
        thelanes.append(lane)


if __name__=="__main__":
    parallelmode = 'nfork'
    if len(sys.argv) == 1:
        lanes=range(nr_lanes)
    else:
        lanes=map(int, sys.argv[1:])
    if parallelmode=='fork':
        main(lanes)
    else:
        main_multiproc(lanes)
