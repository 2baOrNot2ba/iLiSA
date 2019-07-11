#!/usr/bin/python
"""Spec observations.
"""

import os
import argparse
import yaml
import ilisa
from ilisa.observations.stationsession import StationSession

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file',
                        help='station schedule specification file')
    args = parser.parse_args()

    userilisadir = ilisa.user_conf_dir
    accessconffile = os.path.join(userilisadir, 'accessconf_localstation.yml')
    with open(accessconffile) as cfigfilep:
        accessconf = yaml.load(cfigfilep)
    with open(args.file, 'r') as f:
        stn_ses_sched_in = yaml.load(f)
    try:
        mockrun = stn_ses_sched_in['mockrun']
    except KeyError:
        mockrun = False
    stnsess = StationSession(accessconf, mockrun=mockrun)
    stnsess.run_lcl_sched(stn_ses_sched_in)
