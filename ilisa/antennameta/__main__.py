import sys
from .antennafieldlib import antset_lonlat, list_stations

"""
Provides a basic CLI for antennameta package

Run with

bash$ python -m ilisa.antennameta
output> SE607, ..., etc

is the list of all known ILT stations, or

bash$ python -m ilisa.antennameta SE607
output> WGS84 HBA 11.930896027303548deg,57.398761763677555deg,41.35527306050062m
output> Bearing(N->E): 4.27 deg

which is the WGS84 longitude, latitude, and height followed by bearing of
the HBA antennaset of the SE607 station.
"""

if len(sys.argv) > 1:
    stnid = sys.argv[1]
    antset = 'HBA'
    if len(sys.argv) > 2:
        antset = sys.argv[2]
    lon_lat, bearing = antset_lonlat(stnid, antset)
    print(f"WGS84 {antset}: {lon_lat[0]}deg,{lon_lat[1]}deg,{lon_lat[2]}m")
    print("Bearing(N->E):", round(bearing/3.14*180,2), "deg")
else:
    print(*list_stations())