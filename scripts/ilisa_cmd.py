#!/usr/bin/python

import sys
import os
import datetime
import argparse
import yaml
import ilisa
import ilisa.observations.stationdriver as stationdriver
from ilisa.observations.scansession import ScanSession, projid2meta

LOGFILE = "ilisa_cmds.log"


def down(stndrv):
    """Put station into idle state."""
    stndrv.halt_observingstate()


def up(stndrv):
    """Put station into ready to observe state."""
    stndrv.goto_observingstate()


def handback(stndrv):
    """Handback station to ILT control."""
    down(stndrv)


def adm(stndrv, args):
    """Dispatch admin commands."""
    if args.admcmd == 'up':
        up(stndrv)
    elif args.admcmd == 'down':
        down(stndrv)
    elif args.admcmd == 'handback':
        handback(stndrv)


def obs(stndrv, args):
    """Observe a scansession from ScanSes file."""
    with open(args.file, 'r') as f:
        scansess_in = yaml.load(f)
    scansess_in['start'] = args.time
    scansess_in['mockrun'] = args.mockrun
    scansess_in['projectid'] = args.project

    scnsess = ScanSession(stndrv)
    scnsess.run_scansess(scansess_in)


def parse_cmdline(argv):
    """Parse a schedule commandline."""
    cmdln_prsr = argparse.ArgumentParser()
    cmdln_prsr.add_argument('-t', '--time', type=str, default='NOW',
                            help="Start Time (format: YYYY-mm-ddTHH:MM:SS)"
                            )
    cmdln_prsr.add_argument('-p', '--project', type=str, default='0', help="Project ID")
    cmdln_prsr.add_argument('-s', '--station', type=str, default=None, help="Station ID")
    cmdln_prsr.add_argument('-m', '--mockrun', action='store_true', help="Mockrun")
    cmdln_sbprsr = cmdln_prsr.add_subparsers(title='LOFAR stand alone commands',
                                             description='Select a command.',
                                             help='LOFAR commands:')

    parser_adm = cmdln_sbprsr.add_parser('adm', help="Admin")
    parser_adm.set_defaults(cmd='adm')
    parser_adm.set_defaults(func=adm)
    parser_adm.add_argument('admcmd', help='Admin command')

    parser_obs = cmdln_sbprsr.add_parser('obs', help="Observe a scansession.")
    parser_obs.set_defaults(cmd='obs')
    parser_obs.set_defaults(func=obs)
    parser_obs.add_argument('file', help='ScanSession file')
    args = cmdln_prsr.parse_args(argv)
    if args.time == 'NOW':
        args.time = datetime.datetime.utcnow()
    else:
        try:
            args.time = datetime.datetime.strptime(args.time, '%Y-%m-%dT%H:%M:%S')
        except:
            raise RuntimeError("Wrong datetime format.")
    return args


def exec_cmdline(args):
    """Run a schedule commandline."""
    projectmeta, accessfiles = projid2meta(args.project)
    if args.station is None:
        # Try to get station from accessfiles in project config file
        try:
            args.station = accessfiles.keys().pop()
        except:
            raise RuntimeError("No stations found for project {}".format(args.project))
    try:
        # See if station has an access config file
        acf_name = accessfiles[args.station]
    except:
        raise RuntimeError("Station {} not found for project {}".format(args.station,
                                                                        args.project))
    userilisadir = ilisa.observations.user_conf_dir
    acf_path = os.path.join(userilisadir, acf_name)
    with open(acf_path) as acffp:
        ac = yaml.load(acffp)
    # Initialize stationdriver :
    stndrv = stationdriver.StationDriver(ac['LCU'], ac['DRU'], mockrun=args.mockrun)
    args.func(stndrv, args)
    if args.cmd == 'adm':
        args.state = args.admcmd
    elif args.cmd == 'obs':
        args.state = 'obs:' + args.file
    return args


if __name__ == "__main__":
    args = parse_cmdline(sys.argv[1:])
    args = exec_cmdline(args)
    with open(LOGFILE, 'a') as lgf:
        if args.mockrun:
            mockfld = 'M'
        else:
            mockfld = '1'
        lgf.write("{} {} {} {}\n".format(args.time, mockfld, args.project, args.station,
                                          args.state))
