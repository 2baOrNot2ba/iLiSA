#!/usr/bin/python
"""Spec observations.
"""

import argparse
import yaml
import ilisa.observations.session as session


if __name__ == "__main__":
    parser = argparse.ArgumentParser()

    parser.add_argument('file',
                        help='schedule specification file')

    args = parser.parse_args()

    with open(args.file, 'r') as f:
        schedspec = yaml.load(f)
    myses = session.Session(projectid=schedspec['projectid'],
                            halt_observingstate_when_finished = False)
    myses.implement_scanschedule(schedspec)