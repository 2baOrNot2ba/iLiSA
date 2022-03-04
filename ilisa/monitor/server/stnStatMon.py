#!/usr/bin/python3

# Script for aggregating station status information which is broadcast on UDP.
# Inspired by (and partially borrowed from) A. Horneffers "isEcStatus.py" code.
# T. Carozzi 30 Sep 2013 (25 Jan 2012) (TobiaC)
# Modified by A. Horneffer 2013--2017

from stnStatMon_config import *

GW_PC_UDP_PORT = 6070  # Port for istn service
                       # (If you edit this you will need to update port on
                       # clients)

# Paths to various lofar commands:
OPERATIONSPATH = "/opt/operations/bin/"
# OPERATIONSPATH="/usr/local/bin/"
LOFARBINPATH = "/opt/lofar/bin/"

import sys
import socket
import time
from subprocess import Popen, PIPE
from parse_rspctl import parse_rspctl_status, parse_rspctl_rcu
from parse_tbbctl import parse_tbbctl_status

VERSION = '2.4'  # version of this script

STATUS = {}  # Status of LCU


def pathtoISSTATUS():
    """Determine path to ISSTATUS script."""
    lcuSWver = get_lofar_sw_ver()
    change1 = (1, 16, 0)
    change2 = (2, 17, 6)
    changelatest = change2
    if lcuSWver >= changelatest:
        STATIONTESTPATH = '/opt/lofar/sbin/'
        ISSTATUSPY = 'status.py'
    elif lcuSWver >= change1 and lcuSWver < change2:
        STATIONTESTPATH = '/opt/lofar/sbin/'
        ISSTATUSPY='isStatus.py'
    elif lcuSWver < change1:
        STATIONTESTPATH = '/opt/stationtest/test/envcontroltest/'
        ISSTATUSPY = 'isStatus.py'
    else:
        print("Cannot determine version of path")
        sys.exit(1)
        
    return STATIONTESTPATH+ISSTATUSPY


def get_isStatus():
    ISSTATUSscript = pathtoISSTATUS()
    # Environmental control status:
    ECstatOut = Popen(ISSTATUSscript,
                      stdout=PIPE).communicate()[0].decode('UTF8')
    ECstatOutLns = ECstatOut.splitlines()
    STATUS['station'] = ECstatOutLns[1].split()[0]
    STATUS['version'] = ECstatOutLns[1].split()[2][1:-1]
    STATUS['time'] = time.mktime(time.strptime(
        ECstatOutLns[2].lstrip().rstrip()))
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
            key = 'cab3temp'
            value = float(value)
        elif description == 'humidity':
            key = 'cab3hum'
            value = float(value)
        elif description == 'heater state':
            key = 'heater'
        elif description == 'power 48V state':
            key = '48V'
        elif description == 'power LCU state':
            key = 'LCU'
        elif description == 'lightning state':
            key = 'lightning'
        STATUS[key] = value


def get_lofar_sw_ver():
    ps_out=Popen([LOFARBINPATH+'swlevel','-V'],
                 stdout=PIPE).communicate()[0].decode('UTF8')
    verstr=ps_out.split('-')[-1]
    ver_maj, ver_min, ver_pat = [int(ver.strip()) for ver in verstr.split('_')]
    return ver_maj, ver_min, ver_pat


def aggregateInfo():
    """Aggregate information from various sources into a useful iStnMon status.
       The result is in the 'status' global variable."""
    # Get the isStatus status.
    get_isStatus()
    
    if True:
        # Station switch status:
        StnSwtchstatOut = Popen(
           #['sudo', OPERATIONSPATH+'stationswitch','-s'],
           [OPERATIONSPATH+'getstationmode'],
           stdout=PIPE).communicate()[0].decode('UTF8')
        StnSwtchstatOutLns = StnSwtchstatOut.splitlines()
        STATUS['switch'] = StnSwtchstatOutLns[0].split()[-1]
    if True:
        # Software level:
        swlstatOut = Popen([LOFARBINPATH+'swlevel','-S'],
                           stdout=PIPE).communicate()[0].decode('UTF8')
        swlstatOutLns = swlstatOut.splitlines()
        STATUS['softwareLevel'] = int(swlstatOutLns[0].split()[0])
    if CHECK_BC_USER:
        # beamctl user:
        bc_user = who_beamctl()
        STATUS['beamctl']=bc_user
    if CHECK_NUM_LOGINS:
        # number of logged-in users
        woutput = Popen('w | wc', shell=True,
                        stdout=PIPE).communicate()[0].split()[0].decode('UTF8')
        STATUS['allUsers'] = int(woutput)-2
        woutputuser = Popen('w | grep user[0-9] | wc', shell=True,
                          stdout=PIPE
                          ).communicate()[0].split()[0].decode('UTF8')
        STATUS['localUsers'] = int(woutputuser)
    if CHECK_RSP_STATS and (STATUS['softwareLevel'] >= 2):
        # status of the RSP boards
        statoutput = Popen(LOFARBINPATH+'rspctl --status', shell=True,
                           stdout=PIPE,stderr=PIPE
                           ).communicate()[0].decode('UTF8')
        STATUS['rspstat'] = parse_rspctl_status(statoutput)
        rcuoutput = Popen(LOFARBINPATH+'rspctl --rcu', shell=True,
                          stdout=PIPE,stderr=PIPE
                          ).communicate()[0].decode('UTF8')
        STATUS['rcumodes'] = parse_rspctl_rcu(rcuoutput)
    if CHECK_TBB_STATS and (STATUS['softwareLevel'] >= 2):
        # status of the TBBs
        tbbstatoutput = Popen('/opt/lofar/bin/tbbctl --status', shell=True,
                              stdout=PIPE,stderr=PIPE
                              ).communicate()[0].decode('UTF8')
        STATUS['tbbstat'] = parse_tbbctl_status(tbbstatoutput)


def who_beamctl():
    ps_out = Popen(['/bin/ps', '-Cbeamctl', '--no-headers', '-ouser'],
                   stdout=PIPE).communicate()[0].decode('UTF8')
    ps_out_lns = ps_out.splitlines()
    if len(ps_out_lns) == 0:
        bc_user = 'None'
    else:
        bc_user = ps_out_lns[0]
    return bc_user


def printInfo():
    print(STATUS['station'])
    print(STATUS['version'])
    print(STATUS['time'])
    print(STATUS['cab3temp'])
    print(STATUS['cab3hum'])
    print(STATUS['heater'])
    print(STATUS['48V'])
    print(STATUS['LCU'])
    print(STATUS['lightning'])
    if CHECK_BC_USER:
        print(STATUS['beamctl'])


def sendstatus(isUDP=True, isSendTest=False, isLogged=True):
    date = time.localtime(STATUS['time'])
    outstring = "LOFAR_STN_STATUS (version): %s" % VERSION
    outstring += "\n"
    outstring += "%4d-%02d-%02d-%02d:%02d:%02d, "%(
        date.tm_year, date.tm_mon, date.tm_mday,
        date.tm_hour,  date.tm_min, date.tm_sec)
    # Environmental control status:
    outstring += "Station: %s, " % STATUS['station']
    outstring += "ECvers: %s, " % STATUS['version']
    outstring += "Cab3 Temp: %.2fC, " % STATUS['cab3temp']
    outstring += "Cab3 Hum: %.2f%%, " % STATUS['cab3hum']
    outstring += "Heater: %s, " % STATUS['heater']
    outstring += "48V: %s, " % STATUS['48V']
    outstring += "LCU: %s, " % STATUS['LCU']
    outstring += "Lightning: %s" % STATUS['lightning']
    # The switch status:
    outstring += "\n"
    outstring += "Switch: %s" % STATUS['switch']
    # The switch status:
    outstring += "\n"
    outstring += "Software Level: %s" % STATUS['softwareLevel']
    if CHECK_BC_USER:
        # Who is using beamctl:
        outstring += "\n"
        outstring += "beamctl User: %s" % STATUS['beamctl']
    if CHECK_NUM_LOGINS:
        # Logged in users:
        outstring += "\n"
        outstring += "All Users: %s " % STATUS['allUsers']
        outstring += "Local Users: %s " % STATUS['localUsers']
    if 'rcumodes' in STATUS:
        # rcumodes
        outstring += "\n"
        outstring += "RCUmodes -1:%d 0:%d 1:%d 2:%d 3:%d 4:%d 5:%d 6:%d 7:%d" \
                     % (STATUS['rcumodes']['-1'], STATUS['rcumodes']['0'],
                        STATUS['rcumodes']['1'], STATUS['rcumodes']['2'],
                        STATUS['rcumodes']['3'], STATUS['rcumodes']['4'],
                        STATUS['rcumodes']['5'], STATUS['rcumodes']['6'],
                        STATUS['rcumodes']['7'],)
    if 'rspstat' in STATUS:
        outstring += "\n"
        outstring += "RSPtemps PCBmean: %.2fC, BPmean: %.2fC, APmean: %.2fC\n" \
                     % (STATUS['rspstat']['PCBtempmean'],
                        STATUS['rspstat']['BPtempmean'],
                        STATUS['rspstat']['APtempmean'])
        outstring += "RSPtemps PCBmax: %.2fC, BPmax: %.2fC, APmax: %.2fC\n" \
                     % (STATUS['rspstat']['PCBtempmax'],
                        STATUS['rspstat']['BPtempmax'],
                        STATUS['rspstat']['APtempmax'])
        outstring += "RSPtemps PCBmin: %.2fC, BPmin: %.2fC, APmin: %.2fC\n" \
                     % (STATUS['rspstat']['PCBtempmin'],
                        STATUS['rspstat']['BPtempmin'],
                        STATUS['rspstat']['APtempmin'])
        outstring += "RSPvolt V1.2: %.2f, V2.5: %.2f, V3.3: %.2f\n" \
                     % (STATUS['rspstat']['volt12mean'],
                        STATUS['rspstat']['volt25mean'],
                        STATUS['rspstat']['volt33mean'])
        outstring += "RSPvoltMax V1.2: %.2f, V2.5: %.2f, V3.3: %.2f\n" \
                     % (STATUS['rspstat']['volt12max'],
                        STATUS['rspstat']['volt25max'],
                        STATUS['rspstat']['volt33max'])
        outstring += "RSPvoltMin V1.2: %.2f, V2.5: %.2f, V3.3: %.2f" \
                     % (STATUS['rspstat']['volt12min'],
                        STATUS['rspstat']['volt25min'],
                        STATUS['rspstat']['volt33min'])
    if 'tbbstat' in STATUS:
        outstring += "\n"
        outstring += "Bad-TBBs: %d, Good-TBBs: %d\n" \
                     % (STATUS['tbbstat']['badTBBs'],
                        STATUS['tbbstat']['goodlines'])
        outstring += "TBBtemps PCBmean: %.2fC, TPmean: %.2fC, MPmean: %.2fC\n" \
                     % (STATUS['tbbstat']['PCBtempmean'],
                        STATUS['tbbstat']['TPtempmean'],
                        STATUS['tbbstat']['MPtempmean'])
        outstring += "TBBtemps PCBmax: %.2fC, TPmax: %.2fC, MPmax: %.2fC\n" \
                     % (STATUS['tbbstat']['PCBtempmax'],
                        STATUS['tbbstat']['TPtempmax'],
                        STATUS['tbbstat']['MPtempmax'])
        outstring += "TBBtemps PCBmin: %.2fC, TPmin: %.2fC, MPmin: %.2fC\n" \
                     % (STATUS['tbbstat']['PCBtempmin'],
                        STATUS['tbbstat']['TPtempmin'],
                        STATUS['tbbstat']['MPtempmin'])
        outstring += "TBBvolt V1.2: %.2f, V2.5: %.2f, V3.3: %.2f\n" \
                     % (STATUS['tbbstat']['volt12mean'],
                        STATUS['tbbstat']['volt25mean'],
                        STATUS['tbbstat']['volt33mean'])
        outstring += "TBBvoltMax V1.2: %.2f, V2.5: %.2f, V3.3: %.2f\n" \
                     % (STATUS['tbbstat']['volt12max'],
                        STATUS['tbbstat']['volt25max'],
                        STATUS['tbbstat']['volt33max'])
        outstring += "TBBvoltMin V1.2: %.2f, V2.5: %.2f, V3.3: %.2f" \
                     % (STATUS['tbbstat']['volt12min'],
                        STATUS['tbbstat']['volt25min'],
                        STATUS['tbbstat']['volt33min'])

    if isSendTest:
        outstring = "TEST "+outstring

    if isUDP:
        UDP_IP = GW_PC_UDP_IP
        UDP_PORT = GW_PC_UDP_PORT

        sock = socket.socket(socket.AF_INET, # Internet
                             socket.SOCK_DGRAM ) # UDP
        sock.sendto(outstring.encode('UTF8'), (UDP_IP, UDP_PORT))
        if isLogged:
            print(outstring)
    else:
        print(outstring)


from optparse import OptionParser


if __name__ == "__main__":
   parser = OptionParser()
   parser.add_option("-p", "--print", dest="prntout",
                  action="store_true", default=False,
                  help="just print status messages to stdout (no UDP)")
   (options, args) = parser.parse_args()

   aggregateInfo()
   # printInfo()
   sendstatus(isUDP=not(options.prntout))
   bc_user=who_beamctl()
