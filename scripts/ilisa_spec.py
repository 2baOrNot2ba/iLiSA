#!/usr/bin/python
"""Spec observations.
"""

import argparse
import ilisa.observations.session as session


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('datatype',
                        help='lofar datatype rec')
    parser.add_argument('freqbndarg',
                        help='')
    parser.add_argument('duration_tot',
                        help='')
    parser.add_argument('pointsrc',
                        help='')
    args = parser.parse_args()

    myses = session.Session(halt_observingstate_when_finished = False)
    starttime = 'NOW'
    if args.datatype == 'bst':
        obsfun = 'do_bstnew'
    elif args.datatype == 'sst':
        obsfun = 'do_sstnew'
    elif args.datatype == 'xst':
        obsfun = 'do_xstnew'
    obsargs = {'freqbndarg': args.freqbndarg, 'duration_tot': eval(args.duration_tot),
               'integration':1, 'pointsrc': args.pointsrc}
    schedspec = [{'starttime': starttime, 'stations':'*', 'obsfunname': obsfun,
                  'obsargs': obsargs}]

    myses.implement_scanschedule(schedspec)
