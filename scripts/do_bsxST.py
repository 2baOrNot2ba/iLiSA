#!/usr/bin/python
# Records bst,sst,xst data in one of the LOFAR bands and creates a header file
# with observational settings.


import argparse
import ilisa.observations.observing as observing


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("statistic")
    parser.add_argument("frequency", type=float)
    parser.add_argument("integration",type=int)
    parser.add_argument("duration",type=str)
    parser.add_argument("pointSrc", type=str, nargs='?', default='Z')
    args = parser.parse_args()
    duration = int(eval(args.duration))
    myobs = observing.Session()
    # Norminally when the Session is over (i.e. when this script exits), the
    # LCU observing state is stopped. If you want to do further observations
    # you might not want to do this (to avoid the time it takes to boot to
    # an observing state). In this case you should uncomment next line:
    #myobs.halt_observingstate_when_finished = False
    myobs.bits = 16
    myobs.bsxST(args.statistic, args.frequency, args.integration, duration,
         args.pointSrc)

