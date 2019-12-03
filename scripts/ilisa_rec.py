#!/usr/bin/python
"""Record local station data with command-line arguments.
"""

import os
import argparse
import ilisa.observations
import ilisa.observations.directions
import ilisa.observations.stationdriver as stationdriver
import ilisa.observations.modeparms as modeparms
import ilisa.observations.dataIO as dataIO
import ilisa.observations.programs as programs


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-a', '--allsky', help="Set allsky FoV", action='store_true')
    parser.add_argument('-s', '--starttime',
                        help="Start Time (format: YYYY-mm-ddTHH:MM:SS)",
                        type=str, default='NOW')
    parser.add_argument('-i', '--integration',
                        help="Integration time [s]",
                        type=float, default=modeparms.MIN_STATS_INTG)
    parser.add_argument('datatype',
                        help="""lofar data type to record.
                        Choose from 'acc', 'bfs', 'bst', 'sst', 'tbb', 'xst', 'nil'.""")
    parser.add_argument('freqspec',
                        help='Frequency spec in Hz.')
    parser.add_argument('duration_tot',
                        help='Duration in seconds. (Can be an arithmetic expression)',
                        type=str)
    parser.add_argument('pointing', nargs='?',
                        help='A direction in az,el,ref (radians) or a source name.')
    args = parser.parse_args()

    accessconf = ilisa.observations.default_access_lclstn_conf()
    stndrv = stationdriver.StationDriver(accessconf['LCU'], accessconf['DRU'], mockrun=False)
    halt_observingstate_when_finished = True
    stndrv.halt_observingstate_when_finished = halt_observingstate_when_finished
    freqbndobj = modeparms.FrequencyBand(args.freqspec)
    try:
        pointsrc = args.pointing
    except AttributeError:
        pointsrc = None
    try:
        pointing = ilisa.observations.directions.std_pointings(pointsrc)
    except KeyError:
        try:
            phi, theta, ref = pointsrc.split(',', 3)
            pointing = pointsrc
        except ValueError:
            raise ValueError(
                "Error: %s invalid pointing syntax".format(args.pointing))
    duration_tot = eval(str(args.duration_tot))

    do_acc = False
    rec_bfs = False
    scanrecs = {}
    scanrec = dataIO.ScanRecInfo()
    sesspath = accessconf['DRU']['LOFARdataArchive']
    if args.datatype == 'nil':
        rec_stat_type = None
    elif args.datatype == 'acc':
        do_acc = True
        rec_stat_type = None
        scanrecs['acc'] = scanrec
        sesspath = os.path.join(sesspath, 'acc')
    elif args.datatype == 'bfs':
        rec_bfs = True
        rec_stat_type = None
        scanrecs['bfs'] = scanrec
        sesspath = os.path.join(sesspath, 'bfs')
    elif args.datatype == 'bst' or args.datatype == 'sst' or args.datatype == 'xst':
        rec_stat_type = args.datatype
        scanrecs['bsx'] = scanrec
        sesspath = os.path.join(sesspath, rec_stat_type)
    elif args.datatype == 'tbb':
        pass
    else:
        raise RuntimeError('Unknown datatype {}'.format(args.datatype))

    if args.datatype != 'tbb':
        bfdsesdumpdir = accessconf['DRU']['BeamFormDataDir']
        scanmeta = stationdriver.ScanMeta(sesspath, bfdsesdumpdir, scanrecs)
        scanpath = programs.record_scan(stndrv, freqbndobj, args.integration, duration_tot, pointing,
                             pointsrc, starttime=args.starttime,
                             rec_stat_type=rec_stat_type, rec_bfs=rec_bfs,
                             duration_frq=None, do_acc=do_acc, allsky=args.allsky,
                             scanmeta=scanmeta)
    else:
        stndrv.do_tbb(duration_tot, freqbndobj.rcubands[0])
    print "Finished"
    import sys
    sys.stdout.flush()
