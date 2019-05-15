#!/usr/bin/python
"""Spec observations.
"""

import argparse
import ilisa.observations.session as session
import ilisa.observations.modeparms as modeparms


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--allsky', help="Set allsky FoV", action='store_true')
    parser.add_argument('-s', '--starttime',
                        help="Start Time (format: YYYY-mm-ddTHH:MM:SS)",
                        type=str, default='NOW')
    parser.add_argument('-p', '--projectid',
                        help="Project ID",
                        type=str, default='0')
    parser.add_argument('-i', '--integration',
                        help="Integration time [s]",
                        type=float, default=modeparms.MIN_STATS_INTG)

    parser.add_argument('datatype',
                        help='lofar data type to record')
    parser.add_argument('freqspec',
                        help='Frequency spec in Hz.')
    parser.add_argument('duration_tot',
                        help='Duration in seconds. (Can be an arithmetic expression)')
    parser.add_argument('pointing', nargs='?',
                        help='A direction in az,el,ref (radians) or a source name.')
    args = parser.parse_args()

    myses = session.Session(halt_observingstate_when_finished = False)
    beam = {'freqspec': args.freqspec, 'pointing': args.pointing, 'allsky': args.allsky}
    rec_stat = {'type': args.datatype, 'integration': args.integration}
    sessionsched = {'projectid': args.projectid,
                    'stations': 'ALL',
                    'scans': [{'starttime': args.starttime,
                               'beam': beam,
                               'duration_tot': args.duration_tot,
                               'rec_stat': rec_stat
                               }]}
    myses.implement_scanschedule(sessionsched)
