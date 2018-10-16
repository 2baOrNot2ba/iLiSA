"""Provides an instance of your local station. Use it by starting ipython
and then typing:
[in]: execfile('stationprompt.py')

The myobs object then provides interactive interface the local LOFAR station.
"""

#TobiaC (2018-03-09)

import ilisa.observations.observing as observing
import ilisa.observations.dataIO as dataIO

myObsSes = observing.Session(goto_observingstate_when_starting=False)
myObsSes.halt_observingstate_when_finished = False
