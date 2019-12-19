#!/usr/bin/python
"""Spec observations.
"""

import os
import argparse
import yaml
import ilisa.observations
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
    acf_lclstn_name = list(accessfiles.pop().values()).pop()
    userilisadir = ilisa.observations.user_conf_dir
    acf_lclstn_path = os.path.join(userilisadir, acf_lclstn_name)
    with open(acf_lclstn_path) as acffp:
        acf = yaml.load(acffp)
    stnsess = StationSession(acf['LCU'], acf['DRU'], mockrun=mockrun,
                             halt_observingstate_when_finished=False)
    stnsess.run_lcl_sched(stn_ses_sched_in)
