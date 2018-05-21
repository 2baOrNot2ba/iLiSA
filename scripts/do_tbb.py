#!/usr/bin/python
"""Record Transient Buffer Board (TBB) data from one of the LOFAR bands for
duration seconds.
"""


import argparse
import ilisa.observations.observing as observing


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("duration",type=float)
    parser.add_argument("band", type=str)
    args = parser.parse_args()
    myobs = observing.Session()
    myobs.do_tbb(args.duration, args.band)

