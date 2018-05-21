#!/usr/bin/python
"""(NEW, replaces do_stnCalObs.py) Collect data for calibrating the LOFAR
station. Can also be used to sweep through all 512 subbands of an rcumode
provided in one ACC file."""

#TobiaC (2018-01-19)


import argparse
import ilisa.observations.observing as observing


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("band", help="band is 10_90, 110_190, or 210_250")
    parser.add_argument("duration", type=str,
                        help="Duration of calibration obs. in seconds. \
                              Can be an arithmetic formula e.g. 24*60*60.")
    parser.add_argument("pointSrc", type=str, nargs='?', default='Z',
                        help="Pointing direction name. Default 'Z' stands for \
                              zenith.")
    args = parser.parse_args()
    duration = int(eval(args.duration))
    
    myobs = observing.Session()
    myobs.do_acc(args.band, duration, args.pointSrc, exit_obsstate=True)

