#!/usr/bin/python
"""Spec observations.
"""

import argparse
import ilisa.observations.session as session


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('datatype',
                        help='lofar datatype rec')
    parser.add_argument('freqspec',
                        help='')
    parser.add_argument('duration_tot',
                        help='')
    parser.add_argument('pointing',
                        help='')
    args = parser.parse_args()

    myses = session.Session(halt_observingstate_when_finished = False)
    starttime = 'NOW'
    if args.datatype == 'bst':
        obsprog = 'do_bstnew'
    elif args.datatype == 'sst':
        obsprog = 'do_sstnew'
    elif args.datatype == 'xst':
        obsprog = 'do_xstnew'
    beam = {'freqspec': args.freqspec, 'pointing': args.pointing}
    rec_stat = {'duration_tot': eval(args.duration_tot), 'integration': 1}
    sessionsched = {'projectid': '0',
                    'stations': 'ALL',
                    'scans': [{'starttime': starttime,
                               'obsprog': obsprog,
                               'beam': beam,
                               'rec_stat': rec_stat
                               }]}
    myses.implement_scanschedule(sessionsched)
