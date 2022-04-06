#!/usr/bin/python
import logging
import subprocess
import time
import datetime
import argparse
import yaml

import ilisa.operations  # Sets up logging
import ilisa.operations.modeparms as modeparms

__now_timestamp = time.time()
ut2lt_offset = datetime.datetime.fromtimestamp(__now_timestamp) \
               - datetime.datetime.utcfromtimestamp(__now_timestamp)
del __now_timestamp
DEFAULT_PROJ = 0
# Margin of time before starttime 'at' command issued
BEFORE_AT_MARGIN = datetime.timedelta(seconds=8)

def sched2at(schedfile):
    with open(schedfile, 'rb') as file:
        schedtab = yaml.safe_load(file)
    # Setup possible global schedule variables
    mockflag = ''
    mockrun = schedtab.get('mockrun', False)
    if mockrun:
        mockflag = '-m'
    station = schedtab.get('station')
    proj = schedtab.get('project', DEFAULT_PROJ)

    schedlines = schedtab['schedule']
    # Check start times
    for schedline in schedlines:
        if type(schedline['start']) == datetime.datetime:
            start_ut = schedline['start']
        elif type(schedline['start']) == str:
            try:
                start_ut = modeparms.timestr2datetime(schedline['start'])
            except ValueError as verr:
                raise verr
        else:
            raise RuntimeError(type(schedline['start']))
        if start_ut < datetime.datetime.utcnow():
            raise RuntimeError(f'Starttime {start_ut} is in the past.')
        else:
            schedline['start'] = start_ut
    # Now do it for real
    for schedline in schedlines:
        start_ut = schedline['start']
        # Setup start time
        startat_lt = (start_ut + ut2lt_offset - BEFORE_AT_MARGIN)
        sleep_before_at_secs = 0
        if startat_lt.second > 0:
            sleep_before_at_secs = startat_lt.second
            startat_lt = startat_lt.replace(second=0)

        # Setup `at` pipe
        atcmdv = ['at', '-t']
        atcmdv.append(startat_lt.strftime('%Y%m%d%H%M'))
        # Setup project:
        proj = schedline.get('project') if schedline.get('project') else proj
        proj = str(proj)
        # Setup station:
        try:
            station = schedline['station']
        except KeyError:
            if not station:
                raise RuntimeError('No station specified')
        # Setup commands for ilisa_cmd
        cli_name = 'ilisa_cmd'  # Deprecated original monolithic ilisa command
        cli_argv = ['-t',  modeparms.astimestr(schedline['start']),
                    '-s', station, mockflag]
        if schedline['cmd'] == 'obs':
            cli_name = 'ilisa_obs'
            cli_argv += ['-p', proj, schedline.get('session')]
        else:
            cli_name = 'ilisa_adm'
            cli_argv += [schedline['cmd']]
        cmdline = [cli_name] + cli_argv
        cmdline = ' '.join(cmdline)

        # Send shell commands to setup at call
        p = subprocess.Popen(atcmdv, stdin=subprocess.PIPE)
        # First send sleep before ilisa command (at only does down to minutes)
        sleepcmd = ''
        if sleep_before_at_secs:
            sleepcmd = f"sleep {sleep_before_at_secs}\n"
        # Send ilisa_cmd to pipe
        cmdline_wlog = '{} {} >> err_{}.log\n'.format(sleepcmd, cmdline,
                                                      station)
        #print(cmdline_wlog)
        p.communicate(input=cmdline_wlog.encode())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('schedfile', help="Schedule file")
    args = parser.parse_args()
    try:
        sched2at(args.schedfile)
    except (RuntimeError, ValueError) as err:
        logging.error(err)


if __name__ == "__main__":
    main()
