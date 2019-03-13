#!/usr/bin/python
"""Spec observations.
"""

import argparse
import ilisa.observations.session as session


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('file',
                        help='schedule specification file')

    args = parser.parse_args()
    myses = session.Session(
                            halt_observingstate_when_finished = False)
    schedspec = []
    with open(args.file, 'r') as f:
        for l in f:
            if l[0] == '#': continue  # A comment line starting with '#'
            starttimearg, stnsarg, obsfun, obsargs = l.split(' ', 3)
            # if starttimearg == 'NOW':
            #     starttime = datetime.datetime.utcnow()
            # else:
            #     starttime = datetime.datetime.strptime(starttimearg, '%Y%m%dT%H%M%S')
            starttime = starttimearg
            schedspec.append({'starttime': starttime, 'stations': stnsarg,
                              'obsfunname': obsfun, 'obsargs': obsargs})

    myses.implement_scanschedule(schedspec)
