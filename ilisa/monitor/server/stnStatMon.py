#!/usr/bin/python

## Script for aggregating station status information which is broadcast on UDP.
## Inspired by (and partially borrowed from) A. Horneffers "isEcStatus.py" code.
## T. Carozzi 30 Sep 2013 (25 Jan 2012) (TobiaC)
## 
from stnStatMon_config import *

GW_PC_UDP_PORT=6070 #Port for istn service (If you edit this you will need to update port on clients)

#Paths to various lofar commands:
OPERATIONSPATH="/opt/operations/bin/"
#OPERATIONSPATH="/usr/local/bin/"
LOFARBINPATH="/opt/lofar/bin/"

import sys
import socket
import struct
import time
from subprocess import Popen, PIPE
import datetime

VERSION = '2.2' # version of this script

status={}


def pathtoISSTATUS():
    """Determine path to ISSTATUS script."""
    lcuSWver = get_lofar_sw_ver()
    change1 = (1, 16, 0)
    change2 = (2, 17, 6)
    changelatest = change2
    if   lcuSWver >= changelatest:
        STATIONTESTPATH='/opt/lofar/sbin/'
        ISSTATUSPY='status.py'
    elif lcuSWver >= change1 and lcuSWver < change2:
        STATIONTESTPATH='/opt/lofar/sbin/'
        ISSTATUSPY='isStatus.py'
    elif lcuSWver <  change1:
        STATIONTESTPATH='/opt/stationtest/test/envcontroltest/'
        ISSTATUSPY='isStatus.py'
    else:
        print("Cannot determine version of path")
        exit(1)
        
    return STATIONTESTPATH+ISSTATUSPY


def get_isStatus():
    ISSTATUSscript = pathtoISSTATUS()
    #Environmental control status:
    ECstatOut=Popen(
           ISSTATUSscript,
           stdout=PIPE).communicate()[0]
    ECstatOutLns= ECstatOut.splitlines()
    status['station']=ECstatOutLns[1].split()[0]
    status['version']=ECstatOutLns[1].split()[2][1:-1]
    status['time']=time.mktime(time.strptime(ECstatOutLns[2].lstrip().rstrip()))
    key=''
    for lineNr in range(4,10):
        (desc,val)=ECstatOutLns[lineNr].split('=')
        description=desc.rstrip()
        value=val.lstrip()
        # The output labels from isStatus.py were changed 21 June 2013.
        # Effectively the string "cab3" was removed. I changed the
        # matching condition to the output, but not the subsequent
        # naming so as to not have to change anything downstream from here.
        if description == 'temperature':
           key='cab3temp'
           value=float(value)
        elif description == 'humidity':
           key='cab3hum'
           value=float(value)
        elif description == 'heater state':
           key='heater'
        elif description == 'power 48V state':
           key='48V'
        elif description == 'power LCU state':
           key='LCU'
        elif description == 'lightning state':
           key='lightning'
        status[key]=value


def get_lofar_sw_ver():
    ps_out=Popen(
           [LOFARBINPATH+'swlevel','-V'],
           stdout=PIPE).communicate()[0]
    verstr=ps_out.split('-')[-1]
    ver_maj, ver_min, ver_pat = [int(ver.strip()) for ver in verstr.split('_')]
    return ver_maj, ver_min, ver_pat


def aggregateInfo():
    """Aggregate information from various sources into a useful iStnMon status.
       The result is in the 'status' global variable."""
    #Get the isStatus status.
    get_isStatus()
    
    if True:
      #Station switch status:
      StnSwtchstatOut=Popen(
           #['sudo', OPERATIONSPATH+'stationswitch','-s'],
           [OPERATIONSPATH+'getstationmode'],
           stdout=PIPE).communicate()[0]
      StnSwtchstatOutLns= StnSwtchstatOut.splitlines()
      status['switch']=StnSwtchstatOutLns[0].split()[-1]
    if True:
      #Software level:
      swlstatOut=Popen(
           [LOFARBINPATH+'swlevel','-S'],
           stdout=PIPE).communicate()[0]
      swlstatOutLns= swlstatOut.splitlines()
      status['softwareLevel']=int(swlstatOutLns[0].split()[0])
    if CHECK_BC_USER:
      #beamctl user:
      bc_user = who_beamctl()
      status['beamctl']=bc_user


def who_beamctl():
    ps_out=Popen(
           ['/bin/ps', '-Cbeamctl', '--no-headers', '-ouser'],
           stdout=PIPE).communicate()[0]
    ps_out_lns = ps_out.splitlines()
    if len(ps_out_lns) == 0:
        bc_user = 'None'
    else:
        bc_user = ps_out_lns[0]
    return bc_user


def printInfo():
    print status['station']
    print status['version']
    print status['time']
    print status['cab3temp']
    print status['cab3hum']
    print status['heater']
    print status['48V']
    print status['LCU']
    print status['lightning']
    if CHECK_BC_USER:
      print status['beamctl']


def sendstatus(isUDP=True,isSendTest=False,isLogged=True):
        date = time.localtime(status['time'])
        outstring = "LOFAR_STN_STATUS (version): %s" % VERSION
        outstring += "\n"
        outstring += "%4d-%02d-%02d-%02d:%02d:%02d, "%(
            date.tm_year, date.tm_mon, date.tm_mday,
            date.tm_hour,  date.tm_min, date.tm_sec)
        #Environmental control status:
        outstring += "Station: %s, "%status['station']
        outstring += "ECvers: %s, "%status['version']
        outstring += "Cab3 Temp: %.2fC, "%status['cab3temp']
        outstring += "Cab3 Hum: %.2f%%, "%status['cab3hum']
        outstring += "Heater: %s, "%status['heater']
        outstring += "48V: %s, "%status['48V']
        outstring += "LCU: %s, "%status['LCU']
        outstring += "Lightning: %s"%status['lightning']
        #The switch status:
        outstring += "\n"
        outstring += "Switch: %s" % status['switch']
        #The switch status:
        outstring += "\n"
        outstring += "Software Level: %s" % status['softwareLevel']
        if CHECK_BC_USER:
          #Who is using beamctl:
          outstring += "\n"
          outstring += "beamctl User: %s" % status['beamctl']

        if isSendTest:
           outstring="TEST "+outstring
   
        if isUDP:
           UDP_IP=GW_PC_UDP_IP
           UDP_PORT=GW_PC_UDP_PORT

           sock = socket.socket( socket.AF_INET, # Internet
                      socket.SOCK_DGRAM ) # UDP
           sock.sendto( outstring, (UDP_IP, UDP_PORT) )
           if isLogged:
              print outstring 
        else:
           print outstring


from optparse import OptionParser
if __name__ == "__main__":
   parser = OptionParser()
   parser.add_option("-p", "--print",dest="prntout",
                  action="store_true",default=False,
                  help="just print status messages to stdout (no UDP)")
   (options, args) = parser.parse_args()

   aggregateInfo()
   #printInfo()
   sendstatus(isUDP= not( options.prntout))
   bc_user=who_beamctl()
