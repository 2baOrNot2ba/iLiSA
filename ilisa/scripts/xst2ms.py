#!/usr/bin/env python

from __future__ import print_function
import sys
from lofarstation.stationdata import XSTData
from casacore.measures import measures
import ilisa.antennameta.calibrationtables as calibrationtables
import ilisa.operations.data_io as dataIO

XSTfilepath = sys.argv[1]
stnid = "SE607"

obsinfo = dataIO.parse_xstextfilepath(XSTfilepath)
(RA,DEC,ref) = obsinfo['pointing']
pnting = measures().direction(ref, RA+'rad', DEC+'rad')

sd = XSTData(XSTfilepath, subband=int(obsinfo['subband']),
             rcu_mode=int(obsinfo['rcumode']), station_name=stnid,
             integration_time=float(obsinfo['integration']), direction=pnting)

print(sd.time[0])
print("{} MHz".format(sd.frequency / 1e6))

caltabpath = calibrationtables.findcaltabpath(obsinfo['rcumode'], stnid)
print(caltabpath)
sd.set_station_cal(caltabpath)
sd.write_ms(XSTfilepath[:-4]+".ms")
