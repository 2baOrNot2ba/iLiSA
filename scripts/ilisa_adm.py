#!/usr/bin/python

import os
import argparse
import yaml
import ilisa
import ilisa.observations.stationdriver as stationdriver
from ilisa.observations.scansession import projid2meta

def getswlevel(args):
    current_swl = stndrv.lcu_interface.get_swlevel()
    print("The current swlevel of {} is {}".format(stndrv.get_stnid(), current_swl))

def halt(args):
    stndrv.halt_observingstate()

def up(args):
    stndrv.goto_observingstate()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title='Admin commands',
                                       description='Select a type of LCU admin command.',
                                       help='LCU admin commands:')
    parser_swl = subparsers.add_parser('swlevel',
                                       help="Get LCU swlevel.")
    parser_swl.set_defaults(func=getswlevel)
    parser_hlt = subparsers.add_parser('halt',
                                       help="Shutdown LCU observation state.")
    parser_hlt.set_defaults(func=halt)
    parser_up = subparsers.add_parser('up',
                                       help="Put LCU into observation state.")
    parser_up.set_defaults(func=up)

    args = parser.parse_args()
    projectmeta, accessfiles = projid2meta('0')

    # Just first element in list since single station cntrl:
    acf_lclstn_name = accessfiles.values().pop()
    userilisadir = ilisa.observations.user_conf_dir
    acf_lclstn_path = os.path.join(userilisadir, acf_lclstn_name)
    with open(acf_lclstn_path) as acffp:
        acf = yaml.load(acffp)
    # Initialize stationdriver :
    stndrv = stationdriver.StationDriver(acf['LCU'], acf['DRU'])
    args.func(args)
