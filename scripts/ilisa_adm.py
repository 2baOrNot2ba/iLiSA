#!/usr/bin/python

import argparse
from ilisa.observations.session import Session

def getswlevel(args):
    myobs.halt_observingstate_when_finished = False
    for stnid in myobs.stationdrivers:
        current_swl = myobs.stationdrivers[stnid].stationcontroller.get_swlevel()
        print("The current swlevel of {} is {}".format(stnid, current_swl))

def halt(args):
    myobs.set_halt_observingstate()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(title='Admin commands',
                                       description='Select a type of LCU admin command.',
                                       help='LCU admin commands:')
    parser_swl = subparsers.add_parser('swlevel',
                                       help="Get LCU swlevel.")
    parser_swl.set_defaults(func=getswlevel)
    parser_swl = subparsers.add_parser('halt',
                                       help="Shutdown LCU observation state.")
    parser_swl.set_defaults(func=halt)

    args = parser.parse_args()

    myobs = Session(halt_observingstate_when_finished=False)

    args.func(args)
