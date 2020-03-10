#!/usr/bin/python
import subprocess
import time
import datetime
import argparse
import yaml

now_timestamp = time.time()
ut2lt_offset = datetime.datetime.fromtimestamp(now_timestamp) \
               - datetime.datetime.utcfromtimestamp(now_timestamp)
del now_timestamp
DEFAULT_PROJ = 0

def sched2at(schedfile):
    with open(schedfile, 'rb') as file:
        schedtab = yaml.load(file)
    # Setup possible global schedule variables
    mockflag = ''
    try:
        schedtab['mockrun']
    except KeyError:
        pass
    else:
        if schedtab['mockrun']:
            mockflag = '-m'
    try:
        station = schedtab['station']
    except KeyError:
        station = None
    try:
        proj = schedtab['project']
    except KeyError:
        proj = DEFAULT_PROJ

    schedlines = schedtab['schedule']
    for schedline in schedlines:
        # Setup start time
        if type(schedline['start']) == datetime.datetime:
            start_ut = schedline['start']
        elif type(schedline['start']) == str:
            start_ut = datetime.datetime.strptime(schedline['start'], '%Y-%m-%dT%H:%M:%S')
        else:
            raise RuntimeError(type(schedline['start']))
        start_lt =  start_ut + ut2lt_offset
        # Setup `at` pipe
        atcmdv = ['at', '-t']
        atcmdv.append(start_lt.strftime('%Y%m%d%H%M'))
        p = subprocess.Popen(atcmdv, stdin=subprocess.PIPE)
        # Setup project:
        try:
            proj = schedline['project']
        except KeyError:
            pass
        # Setup station:
        try:
            station = schedline['station']
        except KeyError:
            if not station:
                raise RuntimeError('No station specified')
        # Setup commands for ilisa_cmd.py
        schedline_cmd = schedline['cmd']
        try:
            scansesfile = schedline['session']
        except KeyError:
            scansesfile = None
        if scansesfile:
            cmd, arg = 'obs', scansesfile
        else:
            cmd, arg = 'adm', schedline_cmd
        cmdline = 'ilisa_cmd -t {} -p {} -s {} {} {} {}'.format(schedline['start'],
                    proj, station, mockflag, cmd, arg)
        # Send ilisa cmds to pipe
        cmdline_wlog = '{} > err.log\n'.format(cmdline)
        p.communicate(input=cmdline_wlog.encode())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('schedfile', help="Schedule file")
    args = parser.parse_args()
    sched2at(args.schedfile)


if __name__ == "__main__":
    main()
