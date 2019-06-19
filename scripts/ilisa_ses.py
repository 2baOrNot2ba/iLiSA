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
    try:
        projectid = str(schedspec['projectid'])
    except AttributeError:
        raise RuntimeError("'projectid' field not set.")
    try:
        mockrun = schedspec['mockrun']
    except KeyError:
        mockrun = False
    myses = session.Session(projectid=projectid, mockrun=mockrun,
                            halt_observingstate_when_finished = False)
    myses.implement_scanschedule(schedspec)
