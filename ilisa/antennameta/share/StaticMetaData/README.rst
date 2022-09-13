================================
README for StaticMetaData folder
================================

Files in this directory are three different types of station configuration metadata:

* `<stnid>-AntennaField.conf`, are station `<stnid>`'s LBA and HBA element positions
* `<stnid>-iHBADeltas.conf`, are station `<stnid>`'s HBA positions in a tile
* `<stnid>-CableDelays.conf`, are station `<stnid>`'s cable delays.

They are *not* meant to be edited. They are part of iLiSA's package_data.
The reference for these files is the URL:
https://svn.astron.nl/viewvc/LOFAR/trunk/MAC/Deployment/data/StaticMetaData/
in the three folders: `AntennaFields/`, `iHBADeltas/` and `CableDelays/`.
