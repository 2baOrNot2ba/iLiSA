#!/usr/bin/python
import logging
import subprocess
from datetime import datetime, timedelta, timezone
import argparse
import warnings

import yaml

import ilisa.operations  # Sets up logging
import ilisa.operations.modeparms as modeparms
import ilisa.operations.scansession as ssss

# Compute time delta between LT to UT using astimezone(None) to convert to LT
ut2lt_offset = datetime.utcoffset(datetime.now(None).astimezone(None))

DEFAULT_PROJ = 0
# Margin of time before starttime 'at' command issued
BEFORE_AT_MARGIN = timedelta(seconds=8)
ATJOBSFILE_TMPLT = 'at_jobs_{}.txt'  # '{}' tobe replaced by station


def sched2at(schedfile, check=False):
    with open(schedfile, 'rb') as file:
        schedtab = yaml.safe_load(file)
    # Setup possible global schedule variables
    mockflag = ''
    mockrun = schedtab.get('mockrun', False)
    if mockrun:
        mockflag = '-m'
    station = schedtab.get('station')
    proj_def = schedtab.get('project', DEFAULT_PROJ)
    bootbefore = schedtab.get('bootbefore', None)
    duration = schedtab.get('duration', None)
    idleafter = schedtab.get('idleafter', None)
    cyclebefore = schedtab.get('cyclebefore', False)

    first_start_ut = None
    last_start_ut = None

    schedlines = schedtab['schedule']
    if not check and len(schedlines):
        atjobid_file = open(ilisa.operations.USER_CACHE_DIR
                            + ATJOBSFILE_TMPLT.format(station), 'w')
        atjobid_file.write('# Schedulefile: '+schedfile+'\n')
    else:
        atjobid_file = None
    # Check start times
    for schedline in schedlines:
        if type(schedline['start']) == datetime:
            start_ut = schedline['start']
        elif type(schedline['start']) == str:
            try:
                start_ut = modeparms.timestr2datetime(schedline['start'])
            except ValueError as verr:
                raise verr
        else:
            raise RuntimeError(type(schedline['start']))
        start_ut.replace(tzinfo=timezone.utc)  # Add UTC timezone to start time
        _utcnow = datetime.now(timezone.utc)
        if start_ut < _utcnow:
            _mess = 'Starttime {} is in the past.'.format(start_ut)
            if not check:
                raise RuntimeError(_mess)
            warnings.warn(_mess)
        schedline['start'] = start_ut
        if not first_start_ut or first_start_ut > start_ut:
            first_start_ut = start_ut
        if not last_start_ut or last_start_ut < start_ut:
            last_start_ut = start_ut
            last_sess_file = schedline.get('session', None)
            if last_sess_file:
                dur_last_sess = ssss.check_scansess(
                    {'file': last_sess_file})['duration_total']
    # Compute start times for boot and idle, and add them to schedule
    if bootbefore:
        boot_start_ut = first_start_ut - modeparms.hmsstr2deltatime(bootbefore)
        bootschedline ={'start': boot_start_ut, 'cmd': 'boot'}
        schedlines.insert(0, bootschedline)
    if cyclebefore:
        # This forces station to swlevel 0 before booting
        cycle_start_ut = boot_start_ut - modeparms.hmsstr2deltatime('1m')
        cycleschedline ={'start': cycle_start_ut, 'cmd': 'idle'}
        schedlines.insert(0, cycleschedline)
    if duration:
        duration_dt = modeparms.hmsstr2deltatime(duration)
        end_ut = first_start_ut + duration_dt
    else:
        end_ut = last_start_ut + dur_last_sess
    if idleafter:
        idle_start_ut = end_ut + modeparms.hmsstr2deltatime(idleafter)
        idleschedline = {'start': idle_start_ut, 'cmd': 'idle'}
        schedlines.append(idleschedline)
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
        proj = proj_def
        if schedline.get('project', None) and schedline.get('project') > -1:
            proj = schedline.get('project')
        # Setup station:
        try:
            station = schedline['station']
        except KeyError:
            if not station:
                raise RuntimeError('No station specified')
        # Setup commands for ilisa_cmd
        cli_name = 'ilisa_cmd'  # Deprecated original monolithic ilisa command
        cli_argv = ['-t', modeparms.astimestr(schedline['start']),
                    '-s', station, mockflag]
        if schedline['cmd'] == 'obs':
            cli_name = 'ilisa_obs'
            postprocess = schedline.get('postprocess')
            if postprocess:
                cli_argv.extend(['-P', "'{}'".format(postprocess)])
            note = schedline.get('note')
            if note:
                cli_argv.extend(['-n', "'{}'".format(note)])
            cli_argv += ['-p', str(proj), schedline.get('session')]
        else:
            cli_name = 'ilisa_adm'
            cli_argv += [schedline['cmd']]

        cmdline = [cli_name] + cli_argv
        cmdline = ' '.join(cmdline)

        # Send shell commands to setup at call
        if check:
            print(' '.join(atcmdv))
        else:
            p = subprocess.Popen(atcmdv, stdin=subprocess.PIPE,
                                 stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        # First send sleep before ilisa command (at only does down to minutes)
        sleepcmd = ''
        if sleep_before_at_secs:
            sleepcmd = "sleep {}\n".format(sleep_before_at_secs)
        # Send ilisa_cmd to pipe
        cmdline_wlog = '{} {} >> err_{}.log\n'.format(sleepcmd, cmdline,
                                                      station)
        if check:
            print(cmdline_wlog)
        else:
            sout, serr = p.communicate(input=cmdline_wlog.encode())
            at_message = serr.decode('utf-8').split('\n')[-2]
            at_job_id = int(at_message.split(' ')[1])
            atjobid_file.write('{} {}\n'.format(at_job_id, cmdline))
    # Provide info on scheduled start and end times in UT:
    print(f"Scheduled observations start {first_start_ut} and end {end_ut}")
    if bootbefore:
        print(f"Booting at {boot_start_ut}...", end='')
    else:
        print(f"Not booting before observations...", end='')
    if idleafter:
        print(f"and will handback at {idle_start_ut}.")
    else:
        print("and will not go idle afterwards.")
    # If open, close atjobid_file (not in `with` to avoid massive indent)
    if atjobid_file is not None:
        atjobid_file.close()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument('-c','--check', action='store_true')
    parser.add_argument('schedfile', help="Schedule file")
    args = parser.parse_args()
    try:
        sched2at(args.schedfile, args.check)
    except (RuntimeError, ValueError) as err:
        logging.error(err)


if __name__ == "__main__":
    main()
