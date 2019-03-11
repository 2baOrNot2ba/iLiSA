#!/usr/bin/python
import sys
import multiprocessing
import subprocess
import socket
import datetime
import time
import os
import signal

integrateStep = 10
lane = 0  #0,1,2,3
port = 4346
#hostIPs = ['10.211.7.2', '10.212.7.2', '10.213.7.2', '10.214.7.2']
#host = ''
#seqNr=15
logFileName="/mnt/lane0/logs/BFdump.log"
dataDir = "/mnt/lane"
filenameBase = "SE607-Beam_Obs"
filetag4dumpongoing="ONGOING"
ext = ".lbp" #lofar beamlet packets
REC_SET = True
DATA_REDUCE = False
bufsize=1048576
#datapath=dataDir+filenameBase+ext


def logObs(ObsDataFile, bc_args, OneBeamletRec):
    if not OneBeamletRec:
        ObsDataFileBase,ObsDataFileExt=ObsDataFile.split('.')
        logfile="/mnt/lane0/BF/"+ObsDataFileBase+".log"
    else:
        #beamctl log for beamlets
        logfile="/mnt/lane0/logs/Beamlet.log"
    lfh=open(logfile,"a");
    bcline=" ".join(bc_args.split())
    bclines="\\\n".join(bcline.split('\\'))
    if OneBeamletRec:
        lfh.write(ObsDataFile+" => "+"BFlanes "+bclines+"\n")
    else:
        lfh.write(ObsDataFile+' =>\n'+'"'+bclines+'"\n')
    lfh.close()
    print "LCL> created BF logfile: "+logfile


def getAdataBurst(lane):
    try:
        clientsock = socket.socket(socket.AF_INET, socket.SOCK_DGRAM)
    except socket.error as msg:
        print msg 
    try:
        #clientsock.setsockopt(socket.SOL_SOCKET, IN.SO_BINDTODEVICE,"eth10.201"+str(lane+1)+'\0')
        pass
    except socket.error as msg:
        print msg
    #clientsock.bind((hostIPs[lane], port+lane))
    clientsock.bind(('', port+lane))
    #print hostIPs[lane]

#Wait until I get some data:
    print "Lane %d waiting..." % lane
    recv_msg, addr = clientsock.recvfrom(bufsize)
    startTime=datetime.datetime.utcnow()
    print "...got some data at ",startTime,
    if REC_SET==True:
        print "Recording until data-stream times-out."
        #datapath=dataDir+str(lane)+"/"+filenameBase+str(seqNr)+"_"+str(lane)+ext
        datapath=dataDir+str(lane)+"/BF/"+filenameBase+filetag4dumpongoing+"_"+str(lane)+ext
        f=open(datapath, "wb")
        f.write(recv_msg)
    clientsock.settimeout(1.0)
    while True :
        try:
            recv_msg, addr = clientsock.recvfrom(bufsize)
        except socket.timeout:
            clientsock.close()
            if REC_SET==True:
                f.close()
            break
        if REC_SET==True:
            f.write(recv_msg)
    return startTime,datapath


def startlane(lane, threadqueue):
    startTime,tmpdatapath=getAdataBurst(lane)
    stopTime=datetime.datetime.utcnow()
    stopTimeNearSec=stopTime.replace(microsecond=0)
    #lfh=open(logFileName,"r")
    #fLineLst=lfh.readlines()
    #lfh.close()
    #LastEntries=fLineLst[-1].split(",",3)
    #LastSeqNr=int(LastEntries[0])
    #LastLane=int(LastEntries[1])
    #LastStartTime=datetime.datetime.strptime(LastEntries[2],"%Y-%m-%d %H:%M:%S.%f")
    #LastDur_time=time.strptime(LastEntries[3].strip(),"%H:%M:%S.%f")
    #LastDur=datetime.timedelta(seconds=LastDur_time.tm_sec, minutes=LastDur_time.tm_min, hours=LastDur_time.tm_hour)
    #LastStopTime=LastStartTime+LastDur
    #if (stopTime-LastStopTime).total_seconds()>3:
    #  seqNr=LastSeqNr+1
    #else:
    #  seqNr=LastSeqNr
    #lfh=open(logFileName,"a")
    #lfh.write(str(seqNr)+","+str(lane)+","+str(startTime)+","+str(stopTime-startTime)+'\n')
    #lfh.close()

    ObsFileName=filenameBase+"_"+stopTimeNearSec.isoformat()+"_"+str(lane)+ext
    datapath=dataDir+str(lane)+"/BF/"+ObsFileName
    #datapath=dataDir+str(lane)+"/BF/"+filenameBase+str(seqNr)+"_"+str(lane)+ext
    time.sleep(0.0) #Let Mother process find ONGOING file... FIX
    os.rename(tmpdatapath,datapath)
    print "recordBFlanes.py: Created file "+datapath
    threadqueue.put((lane, ObsFileName))
    if DATA_REDUCE:
        subprocess.call(
           ["/home/tobia/lofarBFdataProcess/LofBFraw2bin",
            datapath,str(integrateStep),"o"]
        )
        #subprocess.call(["scp", datapath, "tobia@bele:/data/tobia/xjobb/"], shell=True) 
        #os.remove(datapath)


def main(lanes,BLOCK=True):
    child_pids=[]
    child_lanes=[]
    
    for lane in lanes:
        newpid=os.fork()
        if newpid==0:
            startlane(lane)
            exit()
        else:
            child_pids.append(newpid)
    if BLOCK:
        try :
           while True:
                time.sleep(3600)
        except KeyboardInterrupt:
            for laneNr in range(len(lanes)):
                os.kill(child_pids[laneNr],signal.SIGTERM)
                print "Stopping lane "+str(lanes[laneNr])


laneprocs = []
thelanes = []
threadqueue = multiprocessing.Queue()
def mainmultiproc(lanes):
    for lane in lanes:
        laneprocs.append( multiprocessing.Process(
                              target=startlane, args=(lane,threadqueue)))
        laneprocs[-1].start()
        thelanes.append(lane)


def waitrec():
    ObsFileNames = [None, None, None, None]
    for lane in thelanes:
        (lane, ObsFileName)=threadqueue.get()
        ObsFileNames[lane]=ObsFileName
        #laneproc.join()

    return ObsFileNames


if __name__=="__main__":
    parallelmode = 'nfork'
    if len(sys.argv) == 1:
        lanes=range(4)
    else:
        lanes=map(int, sys.argv[1:])
    if parallelmode=='fork':
        main(lanes)
    else:
        mainmultiproc(lanes)

