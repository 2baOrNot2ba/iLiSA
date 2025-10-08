# iLiSA - international LOFAR in Stand-Alone mode #

iLiSA provides a package to make observations on international LOFAR stations in
stand-alone mode called `observations` and a package to handle metadata about
the configuration of the LOFAR station's antenna metadata called `antennameta`.
LOFAR is a low-frequency radio telescope array located in Europe.

### What is this repository for? ###

iLiSA provides a package to make observations on international LOFAR stations in
stand-alone mode called `operations` and a package to handle metadata about
the configuration of the LOFAR station's antenna metadata called `antennameta`.
LOFAR is a low-frequency radio telescope array located in Europe.

### What is this repository for? ###

iLiSA is of interest to people that have access to a LOFAR station or if you
wish to process data taken with a LOFAR station in stand-alone mode.

The `operations` subpackage allows you to make and record standard stand-alone
observations such as:

* BST - Beamlet statics (Temporal power of beamlets)
* SST - Subband statics (Dynamic spectra of each receiving unit)
* XST - Cross-correlated statistic (Array visibilities for one frequency)
* ACC - Autocovariance cubes (XST stepped through all frequencies)
* TBB - Transient Buffer Board (Direct sampling of each receiving unit)
* BFS - Beam-formed streams (Streamed complex voltages of beamlets)

The package also features special observational modes:

* **Multi-band BST**:
  BST data recorded over multiple bands, i.e. LBA, HBA-lo and HBA-hi
* **SEPTON**:
  Gives HBA an allsky field-of-view. Data taken in this mode is either SST or XST
* **Multi-subband XST**:
  Records XST data in over a sequence of subbands 
 
Observation using iLiSA are made with science-level python functions which issue
the hardware-level commands to the station's Monitoring And Control interface.

The `antennameta` subpackage is of interest if you wish access a LOFAR station's
metadata concerning antennas. This metadata includes:

* AntennaFields - Positions of a station's antenna elements 
* Delays - Cable group delay per antenna element
* CalTables - Complex, per antenna, gain factors to apply on station array
  visibilities

This metadata is useful for postprocessing stand-alone data, e.g. imaging XST
data, handling TBB data or calibration.

The `monitor` subpackage provides tools for collecting information on the state
of a station so that operators can remotely monitor the state of the telescopes.

The `pipelines` subpackage enables the setting up of pipelines that process
data recorded during station operations. It also contains CLI support for
handling LOFAR BFS data files, allowing e.g. conversion to numpy etc.

`scripts` contain some high level scripts such as a scheduling tool that can
start observing sessions at specific times on specific stations and start
pipelines.

Finally, the `calim` subpackage provides functionality to do calibration,
imaging and modeling of radio interferometric observations, with e.g. routines
for handling visibilities and global sky models.

### How do I get set up? ###

* Use setup.py script to install
* (If you want to actually access the station) Setup password-less ssh access to
  LCU via tunnel on the LCU gateway node.
* Requires: `numpy` (for `calim` package also `python-casacore`)
* Optional features might require additional software...

### Status ###

