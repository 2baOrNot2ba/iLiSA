#!/usr/bin/python

import os
import datetime
import argparse
import yaml
import ilisa
import ilisa.observations.stationdriver as stationdriver
from ilisa.observations.scansession import ScanSession, projid2meta


def getswlevel():
    current_swl = stndrv.lcu_interface.get_swlevel()
    print("The current swlevel of {} is {}".format(stndrv.get_stnid(), current_swl))


def halt():
    stndrv.halt_observingstate()


def up():
    stndrv.goto_observingstate()


def adm(args):
    if args.admcmd == 'halt':
        halt()
    elif args.admcmd == 'up':
        up()
    elif args.admcmd == 'swlevel':
        getswlevel()


def obs(args):
    with open(args.file, 'r') as f:
        scansess_in = yaml.load(f)
    scansess_in['start'] = args.time
    scansess_in['mockrun'] = args.mockrun
    scansess_in['projectid'] = args.project

    scnsess = ScanSession(stndrv)
    scnsess.run_scansess(scansess_in)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('-t', '--time', type=str, default='NOW',
                        help="Start Time (format: YYYY-mm-ddTHH:MM:SS)"
                        )
    parser.add_argument('-p', '--project', type=str, default='0', help="Project ID")
    parser.add_argument('-s', '--station', type=str, default=None, help="Station ID")
    parser.add_argument('-m', '--mockrun', action='store_true', help="Mockrun")
    subparsers = parser.add_subparsers(title='LOFAR stand alone commands',
                                       description='Select a command.',
                                       help='LOFAR commands:')

    parser_adm = subparsers.add_parser('adm', help="Admin")
    parser_adm.set_defaults(cmd='adm')
    parser_adm.set_defaults(func=adm)
    parser_adm.add_argument('admcmd', help='Admin command')

    parser_obs = subparsers.add_parser('obs', help="Observe a scansession.")
    parser_obs.set_defaults(cmd='obs')
    parser_obs.set_defaults(func=obs)
    parser_obs.add_argument('file', help='ScanSession file')

    args = parser.parse_args()

    if args.time == 'NOW':
        args.time = datetime.datetime.utcnow()
    else:
        try:
            args.time = datetime.datetime.strptime(args.time, '%Y-%m-%dT%H:%M:%S')
        except:
            raise RuntimeError("Wrong datetime format.")
    projectmeta, accessfiles = projid2meta(args.project)
    if args.station is None:
        try:
            args.station = accessfiles.keys().pop()
        except:
            raise RuntimeError("No stations found for project {}".format(args.project))
    try:
        acf_name = accessfiles[args.station]
    except:
        raise RuntimeError("Station {} not found for project {}".format(args.station))
    print(args.time)
    userilisadir = ilisa.observations.user_conf_dir
    acf_path = os.path.join(userilisadir, acf_name)
    with open(acf_path) as acffp:
        ac = yaml.load(acffp)
    # Initialize stationdriver :
    stndrv = stationdriver.StationDriver(ac['LCU'], ac['DRU'], mockrun=args.mockrun)
    args.func(args)
    print("Performed: {} {} {} {}".format(args.time, args.project, args.station,
                                          args.cmd))
