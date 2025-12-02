import sys
from .antennafieldlib import antset_lonlat, list_stations

"""
Provides a basic CLI for antennameta package

Run with

bash$ python -m ilisa.antennameta
output> SE607, ..., etc

is the list of all known ILT stations, or

bash$ python -m ilisa.antennameta SE607
output> 11.930896027303548deg,57.398761763677555deg,41.35527306050062m

which is the longitude, latitude, and height of the SE607 station.
"""

if len(sys.argv) > 1:
    stnid=sys.argv[1]
    antset = 'HBA'
    lon_lat, noa = antset_lonlat(stnid, antset)
    print(f"{lon_lat[0]}deg,{lon_lat[1]}deg,{lon_lat[2]}m")
else:
    print(*list_stations())