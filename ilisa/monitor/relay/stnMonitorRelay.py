#!/usr/bin/python3

# This script runs on the Gateway PC and
# relays UDP packets from LCU on one ethernet port
# to an external computer.
#
# Tobia Carozzi 2012-10-30

# Import configuration parameters
from relay_config import *

MULTICAST_GROUP = (MULTICAST_ADDRESS, MULTICAST_PORT)

import socket
import select
import struct
import datetime


# struct SHAME_BLOCK {
#   char name[16];
#   int size;
#   int version;
#   int timestamp;
#   int flag;
#   struct MyData myData;
# } block;


# Inbound
sockin = socket.socket(socket.AF_INET,      # Internet
                       socket.SOCK_DGRAM )  # UDP
sockin.bind( (IPin,UDP_PORT) )
sockin.setblocking(0)
output_rate = 1  # output status once every output_rate seconds

# outbound
if RELAYSTNSTAT:
    sockout = socket.socket(socket.AF_INET,      # Internet
                            socket.SOCK_DGRAM ) # UDP
# multicast-out
if SHAMECAST:
    sockmout = socket.socket(socket.AF_INET,  # Internet
                             socket.SOCK_DGRAM, socket.IPPROTO_UDP)  # UDP
    # Set the time-to-live for messages to 1 so they do not go past the
    # local network segment.
    ttl = struct.pack('b', 1)
    sockmout.setsockopt(socket.IPPROTO_IP, socket.IP_MULTICAST_TTL, ttl)

    #mreq = struct.pack("4sl", socket.inet_aton(MULTICAST_ADDRESS), socket.INADDR_ANY)
    #sockmout.setsockopt(socket.IPPROTO_IP, socket.IP_ADD_MEMBERSHIP, mreq)


def stnstat2dict(message):
    """Convert intl station status in the form of an string (as used as message
    between server and cacti client) to a python dict."""
    if message != '':
        # Process status lines
        status_lns = message.splitlines()
        # Get iStnStatus version
        (version_heading, versionstr) = status_lns[0].split(': ')
        version = versionstr.split('.')
        if version >= ['2', '0']:
           (environ_ln, switch_ln, swlevel_ln, bmctluser_ln
            ) = status_lns[1:5]
        elif version >= ['1', '0']:
            (environ_ln, switch_ln, swlevel_ln
             ) = status_lns[1:4]
            bmctluser_ln = ''
        else:
            raise ValueError("Error: unknown version: "+version)
        # Process environment state line
        (timestamp, station_keyval, ECvers_keyval, cab3temp_keyval, cab3hum_keyval
         , heater_keyval, fourty8V_keyval, LCU_keyval, lightning_keyval
         ) = environ_ln.split(',')
        
        # Start creating status dict
        status = {}
        keyvalsep = ': '
        # LOFAR iStnStatus version
        status['version'] = versionstr
        # Environment
        status['time'] = timestamp
        (station_key,  status['station']) =   station_keyval.split(keyvalsep)
        (ECvers_key,   status['ECvers']) =    ECvers_keyval.split(keyvalsep)
        (cab3temp_key, status['cab3temp']) =  cab3temp_keyval.split(keyvalsep)
        (cab3hum_key,  status['cab3hum']) =   cab3hum_keyval.split(keyvalsep)
        (heater_key,   status['heater']) =    heater_keyval.split(keyvalsep)
        (fourty8V_key, status['48V']) =       fourty8V_keyval.split(keyvalsep)
        (LCU_key,      status['LCU']) =       LCU_keyval.split(keyvalsep)
        (lightning_key,status['lightning']) = lightning_keyval.split(keyvalsep)
        # Switch state
        (switch_key, status['switch']) = switch_ln.split(keyvalsep)
        # SW level
        (swlevel_key, status['swlevel']) = swlevel_ln.split(keyvalsep)
        # Beamctl User
        (bmctluser_key, status['beamctl']) = bmctluser_ln.split(keyvalsep)
    else:
        status = ''
    return status

def stnstat2shamecast(status):
    """Convert intl station status in dict form to a shame cast (OSO event
    monitor broadcasting framework) block."""
    timenow = datetime.datetime.now()
    epoch = datetime.datetime.utcfromtimestamp(0)
    secSince1970 = int((timenow-epoch).total_seconds())
    try:
        cab3temp=float(status['cab3temp'].split("C",1)[0])
        cab3hum=float(status['cab3hum'].split("%",1)[0])
    except:
        cab3temp=-1.0
        cab3hum=-1.0
    SHAME_BLOCK_HEAD = '16siiii'
    name = 'LOFAR_SE607'      # =11 char name[16];
    size = 40                 #  int size;
    version = 1               #  int version;
    timestamp = secSince1970  #  int timestamp;
    flag = 0                  #  int flag;
    SHAME_BLOCK_DATA = 'ff'
    endian = '<'  # little
    lofar_scb = struct.pack(endian+SHAME_BLOCK_HEAD+SHAME_BLOCK_DATA,
                            name, size, version, timestamp, flag,
                            cab3temp, cab3hum)
    return lofar_scb

message = ''

while True:
    ready = select.select([sockin], [], [], output_rate)
    if ready[0]:
        message, addr = sockin.recvfrom(1*1024)  # buffer size is 1024 bytes
    if isLogging:
        if logfilename == '':
            stnstatdict = stnstat2dict(message)
            print(stnstatdict)
        else:
            f = open(logfilename, 'w')
            f.write("Latest message:\n")
            f.write(message)
            f.close()
    if message[0:4] == "TEST":
        print("Testing:\n", message)
    if RELAYSTNSTAT:
        sockout.sendto(message, (IPto, UDP_PORT))
    if SHAMECAST:
        # Multicast
        stnstatdict = stnstat2dict(message)
        stnstatshm = stnstat2shamecast(stnstatdict)
        sockmout.sendto(stnstatshm, MULTICAST_GROUP)
