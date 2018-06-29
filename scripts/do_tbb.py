#!/usr/bin/python
"""Record Transient Buffer Board (TBB) data from one of the LOFAR bands for
duration seconds.
"""


import argparse
import ilisa.observations.observing as observing
import ilisa.observations.stationcontrol as stationcontrol


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=
"""Record Transient Buffer Board (TBB) data from one of the LOFAR bands for
duration seconds.""")
    parser.add_argument("band", type=str,
                        help="a LOFAR band: {}.".format(
                        ', '.join(stationcontrol.bandnames)))
    parser.add_argument("duration", type=float,
                        help="duration of tbb dump in seconds.")
    args = parser.parse_args()
    myobs = observing.Session()
    myobs.halt_observingstate_when_finished = False
    myobs.do_tbb(args.duration, args.band)

