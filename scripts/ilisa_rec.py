#!/usr/bin/python
"""Spec observations.
"""

import os
import argparse
import yaml
import ilisa
import ilisa.observations.stationdriver as stationdriver

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file',
                        help='station schedule specification file')
    args = parser.parse_args()

    userilisadir = ilisa.user_conf_dir
    accessconffile = os.path.join(userilisadir, 'access_config.yml')
    with open(accessconffile) as cfigfilep:
        accessconf = yaml.load(cfigfilep)
    with open(args.file, 'r') as f:
        stn_ses_sched_in = yaml.load(f)
    try:
        mockrun = stn_ses_sched_in['mockrun']
    except KeyError:
        mockrun = False
    stndrv = stationdriver.StationDriver(accessconf, mockrun)
    stndrv.run_lcl_session(stn_ses_sched_in)
