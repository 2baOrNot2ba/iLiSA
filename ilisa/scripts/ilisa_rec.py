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
    parser.add_argument(
        '-a', '--allsky', help="Set allsky FoV", action='store_true')
    parser.add_argument('-s', '--starttime',
                        help="Start Time (format: YYYY-mm-ddTHH:MM:SS)",
                        type=str, default='NOW')
    parser.add_argument('-i', '--integration',
                        help="Integration time [s]",
                        type=float, default=modeparms.MIN_STATS_INTG)
    parser.add_argument('datatype',
                        help="""\
lofar data type to record.
Choose from 'acc', 'bfs', 'bst', 'sst', 'tbb', 'xst', 'dmp' or 'None'.""")
    parser.add_argument('freqspec',
                        help='Frequency spec in Hz.')
    parser.add_argument('duration_tot',
                        help='Duration in seconds. '
                             '(Can be an arithmetic expression)',
                        type=str)
    parser.add_argument('pointing', nargs='?', default='Z',
                        help='Direction in az,el,ref (radians) or source name.')
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
        raise ValueError("Invalid pointing syntax: {}".format(args.pointing))
    duration_tot = eval(str(args.duration_tot))

    do_acc = False
    rec_bfs = False
    bsx_type = None
    sesspath = accessconf['DRU']['LOFARdataArchive']
    if args.datatype == 'None':
        # datatype None means running a beam with no recording
        args.datatype = None
    elif args.datatype == 'acc':
        do_acc = True
        sesspath = os.path.join(sesspath, 'acc')
    elif args.datatype == 'bfs':
        rec_bfs = True
        sesspath = os.path.join(sesspath, 'bfs')
    elif (args.datatype == 'bst' or args.datatype == 'sst'
          or args.datatype == 'xst'):
        bsx_type = args.datatype
        sesspath = os.path.join(sesspath, bsx_type)
    elif args.datatype == 'tbb' or args.datatype == 'dmp':
        # 'dmp' is for just recording without setting setting up a beam.
        pass
    else:
        raise RuntimeError('Unknown datatype {}'.format(args.datatype))

    if args.datatype != 'tbb' and args.datatype != 'dmp':
        bfdsesdumpdir = accessconf['DRU']['BeamFormDataDir']
        stndrv.scanpath = sesspath
        use_programs = True
        if use_programs:
            # TODO Remove this block
            scanresult = programs.record_scan(
                stndrv, freqbndobj, duration_tot, pointing,
                starttime=args.starttime, rec=(args.datatype,),
                integration=args.integration, allsky=args.allsky,
                duration_frq=None)
        else:
            dir_bmctl = ilisa.observations.directions.normalizebeamctldir(
                pointing)
            stndrv.streambeams(freqbndobj, dir_bmctl)
            ldatinfo = stndrv.start_scanrec(args.datatype, args.integration,
                                            args.duration_tot)
            scanresult = stndrv.stop_scanrec(ldatinfo, freqbndobj)
        for res in scanresult['rec']:
            print("Saved {} scanrec here: {}".format(
                res, scanresult[res].get_scanrecpath()))
            scanresult[res].write()
        if not scanresult['rec']:
            print("No data recorded ('None' selected)")

    elif args.datatype == 'tbb':
        stndrv.do_tbb(duration_tot, freqbndobj.rcubands[0])
    elif args.datatype == 'dmp':
        stndrv.halt_observingstate_when_finished = False
        stndrv.exit_check = False
        rectime = args.starttime
        lanes = freqbndobj.getlanes()
        band = freqbndobj.rcubands[0]
        scanpath_bfdat = stndrv.bf_data_dir
        stnid = stndrv.get_stnid()
        _datafiles, _logfiles = stndrv.dru_interface.rec_bf_proxy(
            rectime, duration_tot, lanes, band, scanpath_bfdat, stndrv.bf_port0,
            stnid)
    sys.stdout.flush()


if __name__ == "__main__":
    main()
