#!/usr/bin/env python
"""Record local station data with command-line arguments.
"""

import sys
import os
import argparse
import ilisa.observations
import ilisa.observations.directions
import ilisa.observations.stationdriver as stationdriver
import ilisa.observations.modeparms as modeparms
import ilisa.observations.programs as programs


def main():
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
    parser.add_argument('pointing', nargs='?', default='Z',
                        help='A direction in az,el,ref (radians) or a source name.')
    args = parser.parse_args()

    accessconf = ilisa.observations.default_access_lclstn_conf()
    stndrv = stationdriver.StationDriver(accessconf['LCU'], accessconf['DRU'],
                                         mockrun=False)
    halt_observingstate_when_finished = False
    stndrv.halt_observingstate_when_finished = halt_observingstate_when_finished
    freqbndobj = modeparms.FrequencyBand(args.freqspec)
    try:
        pointing = args.pointing
    except AttributeError:
        pointing = None
    dir_bmctl = ilisa.observations.directions.normalizebeamctldir(pointing)
    if not dir_bmctl:
        raise ValueError(
            "Error: %s invalid pointing syntax".format(args.pointing))
    duration_tot = eval(str(args.duration_tot))

    do_acc = False
    rec_bfs = False
    sesspath = accessconf['DRU']['LOFARdataArchive']
    if args.datatype == 'nil':
        bsx_type = None
    elif args.datatype == 'acc':
        do_acc = True
        bsx_type = None
        sesspath = os.path.join(sesspath, 'acc')
    elif args.datatype == 'bfs':
        rec_bfs = True
        bsx_type = None
        sesspath = os.path.join(sesspath, 'bfs')
    elif args.datatype == 'bst' or args.datatype == 'sst' or args.datatype == 'xst':
        bsx_type = args.datatype
        sesspath = os.path.join(sesspath, bsx_type)
    elif args.datatype == 'tbb':
        pass
    else:
        raise RuntimeError('Unknown datatype {}'.format(args.datatype))

    if args.datatype != 'tbb':
        bfdsesdumpdir = accessconf['DRU']['BeamFormDataDir']
        stndrv.scanpath = sesspath
        scanresult = programs.record_scan(
            stndrv, freqbndobj, duration_tot, pointing, starttime=args.starttime,
            rec=[args.datatype], integration=args.integration, allsky=args.allsky,
            duration_frq=None)
        if bsx_type:
            ret_dattyp = 'bsx'
        else:
            ret_dattyp = args.datatype
        print("Saved {} scanrec here: {}".format(args.datatype,
                scanresult[ret_dattyp].get_scanrecpath()))
        scanresult['bsx'].write()
    else:
        stndrv.do_tbb(duration_tot, freqbndobj.rcubands[0])
    print("Finished")
    sys.stdout.flush()


if __name__ == "__main__":
    main()
