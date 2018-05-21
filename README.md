# iLiSA - internatal LOFAR in Stand-Alone mode #

iLiSA provides a package to make observations on international LOFAR stations in
stand-alone mode called observations and a package to handle metadata about
the configuration of the LOFAR station's antenna metadata called antennameta.

### What is this repository for? ###

iLiSA is of interest to people that have access to a LOFAR station or if you wish
to read-in data taken with a LOFAR station in stand-alone mode. LOFAR is a low-frequency
radio telescope array located in Europe.

The "observations" package allows you to make and record standard stand-alone
observations such as:

* BST - Beamlet statics (Temporal power of beamlets)
* SST - Subband statics (Dynamic spectra of each receiving unit)
* XST - Cross-correlated statistic (Array visibilities for one frequency)
* ACC - Autocovariance cubes (XST stepped through all frequencies)
* TBB - Transient Buffer Board (Direct sampling of each receiving unit)
* BFS - Beam-formed streams (Streamed complex voltages of beamlets)

Observation using iLiSA are made with science-level python functions which issue
the hardware-level commands to the station's Monitoring And Control interface.

The "antennameta" is of interest if you wish access a LOFAR station's metadata.
LOFAR is a low-frequency radio telescope array located in Europe.

The metadata of a LOFAR station includes:

* AntennaFields - Positions of a station's antenna elements 
* Delays - Cable group delay per antenna element
* CalTables - Complex gain factors to apply to visiblities from station

This python package is ideal basis for LOFAR station tasks such as
processing stand-alone data (e.g. imaging xst data, handling TBB data).

### How do I get set up? ###

* Use setup.py script to install
* (If you want to actually access the station) Setup password-less ssh access to
  LCU via tunnel on the LCU gateway node.
* Requires: numpy (for some features also python-casacore)
* Optional features might require additional software...

### Status ###

Please note that iLiSA is still in developement. It has been used heavily on the
Swedish LOFAR station. Features such as TBB handling require DAL software.

