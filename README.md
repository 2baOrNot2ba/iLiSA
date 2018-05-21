# LiSA - LOFAR in Stand-Alone mode #

LiSA provides a package to make observations on LOFAR stations in stand-alone
called standalone and a package to handle LOFAR station metadata.

### What is this repository for? ###

LiSA is of interest to people that have access to a LOFAR station or if you wish
to read-in data taken with a LOFAR station in stand-alone mode. LOFAR is a low-frequency
radio telescope array located in Europe.

The "standalone" package allows you to make and record standard stand-alone
observations such as:

* BST - Beamlet statics (Temporal power of beamlets)
* SST - Subband statics (Dynamic spectra of each receiving unit)
* XST - Cross-correlated statistic (Array visibilities for one frequency)
* ACC - Autocovariance cubes (XST stepped through all frequencies)
* TBB - Transient Buffer Board (Direct sampling of each receiving unit)
* BFS - Beam-formed streams (Streamed complex voltages of beamlets)

Observation using LiSA are made with science-level python functions which issue
the hardware-level commands to the station's Monitoring And Control interface.

The "stationmeta" is of interest if you wish access a LOFAR station's metadata.
LOFAR is a low-frequency radio telescope array located in Europe.

The metadata of a LOFAR station includes:

* AntennaFields - Positions of a station's antenna elements 
* Delays - Cable group delay per antenna element
* CalTables - Complex gain factors to apply to visiblities from station

This python package is ideal basis for LOFAR station tasks such as
processing stand-alone data (e.g. imaging xst data, handling TBB data).

### How do I get set up? ###

* Use setup.py script to install
* Setup password-less ssh access to LCU via tunnel on the LCU gateway node.
* Requires: python-casacore, numpy
* Optional features require additional software.

### Status ###

Please note that LiSA is still in developement. It has been used heavily on the
Swedish LOFAR station. Several features require additional software.

