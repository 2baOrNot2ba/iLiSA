#!/usr/bin/python
"""Spec observations.
"""

import os
import argparse
import yaml
import ilisa
from ilisa.observations.stationsession import StationSession, projid2meta

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file',
                        help='station schedule specification file')
    args = parser.parse_args()

    with open(args.file, 'r') as f:
        stn_ses_sched_in = yaml.load(f)
    try:
        mockrun = stn_ses_sched_in['mockrun']
    except KeyError:
        mockrun = False
    try:
        projectid = stn_ses_sched_in['projectid']
    except KeyError:
        projectid = None

    projectmeta, accessfiles = projid2meta(projectid)

    # Just first element in list since single station cntrl:
    acf_lcu, acf_dru = accessfiles[0]['LCU'], accessfiles[0]['DRU']
    userilisadir = ilisa.user_conf_dir
    ac_lcuf = os.path.join(userilisadir, acf_lcu)
    with open(ac_lcuf) as cfigfilep:
        ac_lcu = yaml.load(cfigfilep)
    ac_druf = os.path.join(userilisadir, acf_dru)
    with open(ac_druf) as cfigfilep:
        ac_dru = yaml.load(cfigfilep)

    stnsess = StationSession(ac_lcu, ac_dru, mockrun=mockrun)
    stnsess.run_lcl_sched(stn_ses_sched_in)
