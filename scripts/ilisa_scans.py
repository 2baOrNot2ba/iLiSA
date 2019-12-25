#!/usr/bin/python
"""Spec observations.
"""

import os
import argparse
import yaml
import ilisa.observations
from ilisa.observations.scansession import ScanSession, projid2meta
import ilisa.observations.stationdriver as stationdriver


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('file',
                        help='Scan-session file')
    args = parser.parse_args()

    with open(args.file, 'r') as f:
        scansess_in = yaml.load(f)
    try:
        mockrun = scansess_in['mockrun']
    except KeyError:
        mockrun = False
    try:
        projectid = scansess_in['projectid']
    except KeyError:
        projectid = None

    projectmeta, accessfiles = projid2meta(projectid)

    # Just first element in list since single station cntrl:
    acf_lclstn_name = list(accessfiles.pop().values()).pop()
    userilisadir = ilisa.observations.user_conf_dir
    acf_lclstn_path = os.path.join(userilisadir, acf_lclstn_name)
    with open(acf_lclstn_path) as acffp:
        acf = yaml.load(acffp)
    # Initialize stationdriver :
    stndrv = stationdriver.StationDriver(acf['LCU'], acf['DRU'], mockrun=mockrun,
                                              goto_observingstate_when_starting=False)
    stndrv.halt_observingstate_when_finished = False
    scnsess = ScanSession(stndrv)
    scnsess.run_scansess(scansess_in)
