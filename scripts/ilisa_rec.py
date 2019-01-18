#!/usr/bin/python
"""Start observation on station.
Types of observation are:
* ACC
* BST
* SST
* XST
* TBB

"""

import argparse
import ilisa.observations.session as session


def do_bfs(args):
    myses.do_bfs(args.band, args.duration_tot, args.pointsrc, when=args.begintime,
                 shutdown=args.shutdown)


def do_acc(args):
    myses.do_acc(args.band, args.duration_tot, args.pointsrc)


def do_bst(args):
    """Records BST data in one of the LOFAR bands and creates a header file
    with observational settings.
    """
    myses.do_bsxST('bst', args.freqbnd, args.integration, args.duration_scan, args.pointsrc,
                   when=args.begintime, allsky=args.allsky)


def do_sst(args):
    """Records SST data in one of the LOFAR bands and creates a header file
    with observational settings.
    """
    myses.do_bsxST('sst', args.freqbnd, args.integration, args.duration_scan, args.pointsrc,
                   when=args.begintime, allsky=args.allsky)


def do_xst(args):
    """Records XST data in one of the LOFAR bands and creates a header file
    with observational settings.
    """
    myses.do_bsxST('xst', args.freqbnd, args.integration, args.duration_scan, args.pointsrc,
                   when=args.begintime, allsky=args.allsky)


def do_tbb(args):
    myses.do_tbb(args.duration_tot, args.band)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--allsky', help="Set allsky FoV", action='store_true')
    parser.add_argument('-s', '--shutdown', help="Shutdown observing state when finished",
                        action='store_true')
    parser.add_argument('-b', '--begintime',
                        help="Begin Time (format: YYYY-mm-ddTHH:MM:SS)",
                        type=str, default='NOW')
    parser.add_argument('-p', '--projectid',
                        help="Project ID (format: YYYY-mm-ddTHH:MM:SS)",
                        type=str, default=None)
    subparsers = parser.add_subparsers(title='Observation mode',
                                       description='Select a type of data to record.',
                                       help='Type of datataking:')

    # Specify common parameter args:
    arg_rcuband_kwargs = {'type': str,
                       'help': """(RCU) Band to use: 10_90, 30_90, 110_190, \
                               170_230 or 210_250."""}
    arg_freqband_kwargs = {'type': str,
                           'help': """\
                            Frequency band spec in Hz.\
                            Format: freqct | freqlo:freqhi | freqlo:freqstp:freqhi\
                            """}
    arg_integration_kwargs ={'type': int,
                             'help': "Integration time in s"}
    arg_dur_st_kwargs = {'type': str,
                         'help': "Duration of statistic scan in seconds. \
                                  Can be an arithmetic formula e.g. 24*60*60."}
    arg_dur_tot_kwargs = {'type': str,
                           'help': "Total duration of session in seconds. \
                                    Can be an arithmetic formula e.g. 24*60*60."}
    arg_pointsrc_kwargs = {'type': str, 'nargs': '?', 'default': 'Z',
                           'help': "Pointing [format: RA,DEC,REF with RA,DEC in radians, \
                                    REF=J2000] or a source name. Default 'Z' stands for \
                                    zenith."}

    # ACC data
    parser_acc = subparsers.add_parser('acc',
                                       help="Make an ACC observation.")
    parser_acc.set_defaults(func=do_acc)
    parser_acc.add_argument('band', **arg_rcuband_kwargs)
    parser_acc.add_argument('duration_tot', **arg_dur_tot_kwargs)
    parser_acc.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # BFS data
    parser_bfs = subparsers.add_parser('bfs',
                                       help="Make an BFS observation.")
    parser_bfs.set_defaults(func=do_bfs)
    parser_bfs.add_argument('band', **arg_rcuband_kwargs)
    parser_bfs.add_argument('duration_tot',**arg_dur_tot_kwargs)
    parser_bfs.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # BST data
    parser_bst = subparsers.add_parser('bst',
                                       help="Make a BST observation")
    parser_bst.set_defaults(func=do_bst)
    parser_bst.add_argument('freqbnd', **arg_freqband_kwargs)
    parser_bst.add_argument('integration', **arg_integration_kwargs)
    parser_bst.add_argument('duration_scan',**arg_dur_st_kwargs)
    parser_bst.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # SST data
    parser_sst = subparsers.add_parser('sst',
                                       help="Make a SST observation")
    parser_sst.set_defaults(func=do_sst)
    parser_sst.add_argument('freqbnd', **arg_freqband_kwargs)
    parser_sst.add_argument('integration',**arg_integration_kwargs)
    parser_sst.add_argument('duration_scan',**arg_dur_st_kwargs)
    parser_sst.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # XST data
    parser_xst = subparsers.add_parser('xst',
                                       help="Make a XST observation")
    parser_xst.set_defaults(func=do_xst)
    parser_xst.add_argument('freqbnd', **arg_freqband_kwargs)
    parser_xst.add_argument('integration',**arg_integration_kwargs)
    parser_xst.add_argument('duration_scan',**arg_dur_st_kwargs)
    parser_xst.add_argument('pointsrc', **arg_pointsrc_kwargs)

    # TBB data
    parser_tbb = subparsers.add_parser('tbb',
                                       help="Make a TBB observation")
    parser_tbb.set_defaults(func=do_tbb)
    parser_tbb.add_argument('band', **arg_rcuband_kwargs)
    parser_tbb.add_argument('duration_tot', **arg_dur_tot_kwargs)

    args = parser.parse_args()
    args.duration_scan = eval(args.duration_scan)
    args.duration_tot = eval(args.duration_tot)

    myses = session.Session(projectid=args.projectid,
                            halt_observingstate_when_finished = args.shutdown)

    args.func(args)
