#!/usr/bin/python3
"""Python based recorder for beamformed data streams."""
import multiprocessing
import subprocess
import socket
import datetime
import os
import os.path
import argparse


integrate_step = 10
nr_lanes = 4
# port = 4346
# datadir_template = "/mnt/lane?/BF/SE607/"
# basedir = "Scans"
filename_base = ""
ext = "lbp"  # lofar beamlet packets
REC_SET = True
DATA_REDUCE = False
bufsize = 1048576


def startlane(lane, port, dumppath, dur=None, threadqueue=None):
    """Monitor a streaming lane of beamformed data."""

    def get_databurst():
        try:
            clientsock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        except socket.error as msg:
            raise msg
        clientsock.bind(('', port))  # clientsock.bind((hostIPs[lane], port+lane))

        # Wait until I get some data:
        print("Lane %d waiting..." % lane)
        recv_msg, addr = clientsock.recvfrom(bufsize)
        start_time = datetime.datetime.utcnow()
        start_time_near_sec = start_time.replace(microsecond=0)
        print("...got some data at {}".format(start_time))
        datapath = None
        if REC_SET:
            print("Recording until data-stream times-out.")
            datapath = os.path.join(dumppath,
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
            elapsedtime = datetime.datetime.utcnow() - start_time
            if dur:
                if elapsedtime > datetime.timedelta(seconds=int(dur)):
                    break
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


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--port0',
                        type=int, default=4346,
                        help = "Port number for lane 0"
                        )
    parser.add_argument('-b', '--bfdatadir',
                        type=str, default='/mnt/lane?/BF/SE607/',
                        help="Template directory for BF data"
                        )
    parser.add_argument('-d', '--dur',
                        type=int, default=None,
                        help="Duration of recording in seconds"
                        )
    args = parser.parse_args()
    main(args.port0, args.bfdatadir, args.dur)


def main(port0, datadir, duration):
    lanes = range(nr_lanes)
    retvalq = multiprocessing.Queue()
    datafiles = [None]*nr_lanes
    logfiles = [None]*nr_lanes
    laneprocs = []
    for lane in lanes:
        dumppath = datadir.replace('?', str(lane))
        if not os.path.exists(dumppath):
            os.mkdir(dumppath)
        port = port0 + lane
        laneproc = multiprocessing.Process(target=startlane, args=(lane, port, dumppath,
                                                                   duration, retvalq))
        laneproc.start()
        laneprocs.append(laneproc)
    for laneproc in laneprocs:
        laneproc.join()
        lane_ret, datafileguess = retvalq.get()
        datafiles[lane_ret] = datafileguess
        dumplogname = None
        logfiles.append(dumplogname)
    return datafiles, logfiles


if __name__=="__main__":
    cli_main()
