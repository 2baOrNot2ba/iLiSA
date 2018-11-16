#!/usr/bin/python
"""Start observation on station.
Types of observation are:
* ACC
* BST
* SST
* XST
* TBB

"""

import math
import argparse
import ilisa.observations.observing as observing


def do_acc(args):
    """Record acc data for one of the LOFAR bands over a duration.
    """
    duration = int(eval(args.duration))
    accdestdir = myobs.do_acc(args.band, duration, args.pointsrc, exit_obsstate=False)
    print("Saved ACC data in folder: {}".format(accdestdir))


def _do_bsx(statistic, args):
    """Records bst,sst,xst data in one of the LOFAR bands and creates a header file
    with observational settings.
    """
    duration = int(math.ceil(eval(args.duration)))
    myobs.bits = 16
    if args.allsky and args.freqlo > 100.0e6 :
        myobs.do_SEPTON(statistic, args.freqlo, args.integration, duration)
    else:
        myobs.bsxST(statistic, args.freqlo, args.integration, duration,
                    args.pointsrc)


def do_bst(args):
    """Records BST data in one of the LOFAR bands and creates a header file
    with observational settings.
    """
    _do_bsx('bst', args)


def do_sst(args):
    """Records SST data in one of the LOFAR bands and creates a header file
    with observational settings.
    """
    _do_bsx('sst', args)


def do_xst(args):
    """Records XST data in one of the LOFAR bands and creates a header file
    with observational settings.
    """
    _do_bsx('xst', args)


def do_tbb(args):
    """Record Transient Buffer Board (TBB) data from one of the LOFAR bands for
    duration seconds.
    """
    myobs.do_tbb(args.duration, args.band)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--allsky', help="Set allsky FoV", action='store_true')
    subparsers = parser.add_subparsers(title='Observation mode',
                                       description='Select a type of data to record.',
                                       help='Type of datataking:')

    # Specify common parameter args:
    arg_band_kwargs = {'type': str,
                       'help': "Band to use: 10_90, 110_190, or 210_250."}
    arg_freqlo_kwargs = {'type': float,
                         'help': "Frequency in Hz (lower edge)"}
    arg_integration_kwargs ={'type': int,
                             'help': "Integration time in s"}
    arg_duration_kwargs = {'type': str,
                           'help': "Duration of calibration obs. in seconds. \
                                    Can be an arithmetic formula e.g. 24*60*60."}
    arg_pointsrc_kwargs = {'type': str, 'nargs': '?', 'default': 'Z',
                           'help': "Pointing [format: RA,DEC,REF with RA,DEC in radians, \
                                    REF=J2000] or a source name. Default 'Z' stands for \
                                    zenith."}

    # ACC data
    parser_acc = subparsers.add_parser('acc',
                                       help="Make an ACC observation.")
    parser_acc.set_defaults(func=do_acc)
    parser_acc.add_argument('band', **arg_band_kwargs)
    parser_acc.add_argument('duration', **arg_duration_kwargs)
    parser_acc.add_argument("pointsrc", **arg_pointsrc_kwargs)

    # BST data
    parser_bst = subparsers.add_parser('bst',
                                       help="Make a BST observation")
    parser_bst.set_defaults(func=do_bst)
    parser_bst.add_argument('freqlo', **arg_freqlo_kwargs)
    parser_bst.add_argument('integration', **arg_integration_kwargs)
    parser_bst.add_argument('duration',**arg_duration_kwargs)
    parser_bst.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # SST data
    parser_sst = subparsers.add_parser('sst',
                                       help="Make a SST observation")
    parser_sst.set_defaults(func=do_sst)
    parser_sst.add_argument('freqlo', **arg_freqlo_kwargs)
    parser_sst.add_argument('integration',**arg_integration_kwargs)
    parser_sst.add_argument('duration',**arg_duration_kwargs)
    parser_sst.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # XST data
    parser_xst = subparsers.add_parser('xst',
                                       help="Make a XST observation")
    parser_xst.set_defaults(func=do_xst)
    parser_xst.add_argument('freqlo', **arg_freqlo_kwargs)
    parser_xst.add_argument('integration',**arg_integration_kwargs)
    parser_xst.add_argument('duration',**arg_duration_kwargs)
    parser_xst.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # TBB data
    parser_tbb = subparsers.add_parser('tbb',
                                       help="Make a TBB observation")
    parser_tbb.set_defaults(func=do_tbb)
    parser_tbb.add_argument('band', **arg_band_kwargs)
    parser_tbb.add_argument('duration', **arg_duration_kwargs)

    args = parser.parse_args()

    myobs = observing.Session()
    # Norminally when the Session is over (i.e. when this script exits), the
    # LCU observing state is stopped. If you want to do further observations
    # you might not want to do this (to avoid the time it takes to boot to
    # an observing state). In this case you should uncomment next line:
    myobs.halt_observingstate_when_finished = False

    args.func(args)
