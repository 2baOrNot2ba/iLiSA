#!/usr/bin/python
import sys
import multiprocessing
import subprocess
import socket
import datetime
import time
import os
import os.path
import signal

integrate_step = 10
nr_lanes = 4
port = 4346
datadir_template = "/mnt/lane?/BF/SE607/"
basedir = "Scans"
filename_base = ""
ext = "lbp"  # lofar beamlet packets
REC_SET = True
DATA_REDUCE = False
bufsize = 1048576


def startlane(lane, threadqueue = None):
    """Monitor a streaming lane of beamformed data."""

    def get_databurst():
        try:
            clientsock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        except socket.error as msg:
            raise msg
        clientsock.bind(('', port + lane))  # clientsock.bind((hostIPs[lane], port+lane))

        # Wait until I get some data:
        print("Lane %d waiting..." % lane)
        recv_msg, addr = clientsock.recvfrom(bufsize)
        start_time = datetime.datetime.utcnow()
        start_time_near_sec = start_time.replace(microsecond=0)
        print("...got some data at {}".format(start_time))
        datapath = None
        if REC_SET:
            print("Recording until data-stream times-out.")
            maindatadir = datadir_template.replace('?', str(lane))
            datapath = os.path.join(maindatadir, basedir,
                                    "{}{}_{}.{}".format(filename_base,
                                                        start_time_near_sec.isoformat(),
                                                        lane, ext))
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

    start_time, recdatapath = get_databurst()
    stop_time = datetime.datetime.utcnow()

    print("{}: Created file: {}".format(os.path.basename(__file__), recdatapath))
    if threadqueue:
        threadqueue.put((lane, recdatapath))
    if DATA_REDUCE:
        subprocess.call(
           ["LofBFraw2bin",
            recdatapath, str(integrate_step), "o"]
        )


def main_forking(lanes, block=True):
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
        main_forking(lanes)
    else:
        main_multiproc(lanes)
