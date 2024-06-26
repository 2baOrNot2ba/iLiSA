#!/usr/bin/python3
"""Record beamformed streams."""
import multiprocessing
import subprocess
import socket
import datetime
import os
import os.path
import argparse
import platform


integrate_step = 10
nr_lanes = 4
# port = 4346
# datadir_template = "/mnt/lane?/BF/SE607/"
# basedir = "Scans"
filename_base = "udp"
ext = "lbp"  # lofar beamlet packets
REC_SET = True
DATA_REDUCE = False
bufsize = 1048576
DRU_NAME = platform.node()


def startlane(lane, port, dumppath, dur=None, stnid='', threadqueue=None):
    """
    Python based recorder for beamformed data streams

    Monitor a streaming lane of beamformed data and save data to file.
    """

    def get_databurst():
        try:
            clientsock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
        except socket.error as msg:
            raise msg
        clientsock.bind(('', port))  # clientsock.bind((hostIPs[lane], port+lane))

        # Wait until I get some data:
        print("Lane %d waiting..." % lane)
        recv_msg, _addr = clientsock.recvfrom(bufsize)
        start_time = datetime.datetime.utcnow()
        start_time_near_sec = start_time.replace(microsecond=0)
        print("...got some data at {}".format(start_time))
        datapath = None
        if REC_SET:
            if not dur:
                print("Recording until data-stream times-out.")
            else:
                print("Recording {} s of data.".format(dur))
            datapath = os.path.join(dumppath,
                "{}_{}_{}.{}.{}.{}".format(filename_base, stnid, port, DRU_NAME,
                                           start_time_near_sec.isoformat(), ext))
            f = open(datapath, "wb")
            f.write(recv_msg)
        clientsock.settimeout(1.0)
        while True:
            try:
                recv_msg, _addr = clientsock.recvfrom(bufsize)
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

    _start_time, recdatapath = get_databurst()
    _stop_time = datetime.datetime.utcnow()

    print("{}: Created file: {}".format(os.path.basename(__file__), recdatapath))
    if threadqueue:
        threadqueue.put((lane, recdatapath))
    if DATA_REDUCE:
        subprocess.call(
           ["LofBFraw2bin",
            recdatapath, str(integrate_step), "o"]
        )


def main(ports, datadir, duration, stnid=None):
    _nrlanes = len(ports)
    lanes = range(_nrlanes)
    retvalq = multiprocessing.Queue()
    datafiles = [None]*_nrlanes
    logfiles = [None]*_nrlanes
    laneprocs = []
    for lane in lanes:
        dumppath = datadir.replace('?', str(lane))
        if not os.path.exists(dumppath):
            os.mkdir(dumppath)
        port = ports[lane]
        laneproc = multiprocessing.Process(target=startlane,
                                           args=(lane, port, dumppath,
                                                 duration, stnid, retvalq))
        laneproc.start()
        laneprocs.append(laneproc)
    for laneproc in laneprocs:
        laneproc.join()
        lane_ret, datafileguess = retvalq.get()
        datafiles[lane_ret] = datafileguess
        dumplogname = None
        logfiles.append(dumplogname)
    return datafiles, logfiles


def cli_main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--starttime',
                        type=str, default='NOW',
                        help = "Start-time,: (iso format) YYYY-mm-ddTHH:MM:SS"
                        )
    parser.add_argument('-p', '--port0',
                        type=int, default=4346,
                        help = "Port number for lane 0"
                        )
    parser.add_argument('-b', '--bfdatadir',
                        type=str, default='/mnt/lane?/BF/SE607/Scans/',
                        help="Template directory for BF data"
                        )
    parser.add_argument('-d', '--dur',
                        type=int, default=None,
                        help="Duration of recording in seconds"
                        )
    parser.add_argument('-s', '--stnid',
                        type=str, default='',
                        help="Station ID",
                        )
    args = parser.parse_args()
    main(args.port0, args.bfdatadir, args.dur, args.stnid)


if __name__=="__main__":
    cli_main()
