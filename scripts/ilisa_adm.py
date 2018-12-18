#!/usr/bin/python

import argparse
import ilisa.observations.stationdriver as observing

def getswlevel(args):
    myobs.halt_observingstate_when_finished = False
    current_swl = myobs.stationcontroller.getswlevel()
    print("The current swlevel is {}".format(current_swl))

def halt(args):
    myobs.halt_observingstate_when_finished = True

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

    myobs = observing.Session(goto_observingstate_when_starting=False)

    args.func(args)
