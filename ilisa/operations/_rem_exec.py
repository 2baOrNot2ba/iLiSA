import subprocess
import os
import datetime
try:
    import paramiko
    IMPORTED_PARAMIKO = True
except ImportError:
    IMPORTED_PARAMIKO = False
import logging
import ilisa.operations

_USE_SSH_MASTER_MODE = True  # Whether to use SSH's Master mode option
_LOGGER = logging.getLogger(__name__)


def _exec_rem(remnode, cmdline, nodetype='LCU', background_job=False,
              dryrun=False, accessible=False, quotes="'", stdoutdir='~',
              verbose=True, _slow_conn_time=0):
    return _exec_ssh(remnode, cmdline, nodetype=nodetype,
                     background_job=background_job, dryrun=dryrun,
                     accessible=accessible, quotes=quotes, stdoutdir=stdoutdir,
                     verbose=verbose, _slow_conn_time=_slow_conn_time)


def _exec_ssh(nodeurl, cmdline, nodetype='LCU',
              background_job=False, dryrun=False, accessible=False, quotes="'",
              stdoutdir='~', verbose=True, _slow_conn_time=0):
    """Execute a command on the remnode, either as a background job or in the
    foreground (blocking). Typically access is remote via ssh.
    (To speed things up use the ssh CommandMaster option.)
    """
    _log_exec_exit = False  # Add info line in logging when exec exits
    if nodetype == 'TEST':
        nodetype = nodeurl
    nodeprompt = "On " + nodetype + "> "
    if nodeurl.endswith('localhost'):
        shellinvoc = ''
        quotes = ''
    else:
        shellinvoc = "ssh"
        if _USE_SSH_MASTER_MODE:
            shellinvoc += " -o ControlMaster=auto"  # Set ControlMaster
            shellinvoc += " -o ControlPersist=3600"  # Set
            shellinvoc += (" -S " + ilisa.operations.USER_CACHE_DIR
                           + "/cm-%r@%h:%p")  # Set unique socket in cache dir
        shellinvoc += ' ' + nodeurl
    output = None
    if _slow_conn_time:
        # Simulate a slow connection by executing a sleep before cmdline:
        cmdline = "sleep {};".format(_slow_conn_time) + " " + cmdline
    if background_job:
        # Currently only run_beamctl & run_tbbctl run in background
        # Put stdout & stderr in log in dumpdir
        cmdline = ("(( " + cmdline + " ) > " + stdoutdir + "/"
                   + "lcu_shell_out.log 2>&1) &")
    if dryrun:
        pre_prompt = "(dryrun) "
    else:
        pre_prompt = ""
    if verbose:
        _LOGGER.info(pre_prompt + nodeprompt + cmdline)
    if (not dryrun) and accessible:
        if background_job == 'locally':
            # Runs in background locally rather than in background on LCU
            output = subprocess.run(shellinvoc + " " + cmdline + " &",
                                    shell=True, stdout=subprocess.PIPE).stdout
        else:
            output = subprocess.run(shellinvoc + " "
                                    + quotes + cmdline + quotes,
                                    shell=True, universal_newlines=True,
                                    stdout=subprocess.PIPE).stdout
            if _log_exec_exit:
                _LOGGER.info('End {} {}'.format(nodetype, cmdline))
        if output:
            output = output.rstrip()
    elif not accessible:
        _LOGGER.warning("Not running as " + nodeurl
                        + " since it is not accessible.")
    return output


def __exec_lcu_paramiko(self, cmdline, background_job=False):
    lcuprompt = "LCUp>"
    if self.DryRun:
        preprompt = "(dryrun)"
    else:
        preprompt = ""
    if background_job is True:
        cmdline = "(( " + cmdline + " ) > " + self._home_dir +\
                  "lofarctl.log 2>&1) &"
    if self.verbose:
        print("{} {} {}".format(preprompt, lcuprompt, cmdline))

    client = paramiko.SSHClient()
    client.load_system_host_keys()
    client.set_missing_host_key_policy(paramiko.WarningPolicy())

    ssh_config = paramiko.SSHConfig()
    user_config_file = os.path.expanduser("~/.ssh/config")
    if os.path.exists(user_config_file):
        with open(user_config_file) as f:
            ssh_config.parse(f)
    cfg = {'hostname': self.hostname, 'username': self.user}

    user_config = ssh_config.lookup(cfg['hostname'])
    for k in ('hostname', 'username', 'port'):
        if k in user_config:
            cfg[k] = user_config[k]

    if 'proxycommand' in user_config:
        cfg['sock'] = paramiko.ProxyCommand(user_config['proxycommand'])

    client.connect(**cfg)

    stdin, stdout, stderr = client.exec_command(cmdline)
    print(stdout.read())
    client.close()


def __stdout_ssh(nodeurl, cmdline, nodetype='LCU', dryrun=False,
                 verbose=True):
    """Execute a command on the remnode and return its output."""
    nodeprompt = "On " + nodetype + "> "
    shellinvoc = "ssh " + nodeurl
    if dryrun:
        pre_prompt = "(dryrun) "
    else:
        pre_prompt = ""
    if verbose:
        print(pre_prompt + nodeprompt + cmdline)
    if not dryrun:
        try:
            output = subprocess.check_output(shellinvoc + " '" + cmdline + "'",
                                             shell=True).rstrip()
            output = str(output.decode('UTF8'))
        except subprocess.CalledProcessError as e:
            raise e
    else:
        output = "None"
    return output


def __outfromLCU(self, cmdline, integration, duration):
    """Execute a command on the LCU and monitor progress."""
    LCUprompt = "LCUo> "
    shellinvoc = "ssh " + self.lcu
    if self.DryRun:
        prePrompt = "(dryrun) "
    else:
        prePrompt = ""
    if self.verbose:
        print(prePrompt+LCUprompt+cmdline)
    if self.DryRun is False:
        cmd = subprocess.Popen(shellinvoc+" '"+cmdline+"'",
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, shell=True)
    else:
        return None
    count = 0
    outstrname = 'stderr'
    while cmd.poll() is None:
        if outstrname == 'stdout':
            outstr = cmd.stdout
        elif outstrname == 'stderr':
            outstr = cmd.stderr
        else:
            raise ValueError("Unknown output name {}".format(outstrname))
        try:
            got = cmd.stderr.readline().decode('utf8')
        except IOError:
            raise IOError()
        else:
            # print got
            if "shape(stats)=" in got:
                if count % 2 == 0:
                    print(str(int(round(duration-count/2.0*integration, 0)
                                  )) + "sec left out of " + str(duration))
                count += 1


def ostimenow(url, remunit):
    """\
    Return now time of remunit's operating system
    """
    def _rem_exec(cmdline):
        return _exec_ssh(nodeurl=url, cmdline=cmdline, nodetype=remunit,
                         accessible=True)
    isofmt = ilisa.operations.DATETIMESTRFMT
    datecmdline = "date -u '+{}'".format(isofmt)
    lcunow_str = _rem_exec(cmdline=datecmdline)
    print(lcunow_str)
    lcunow = datetime.datetime.strptime(lcunow_str, isofmt)
    # Since linux date returns nanosecond time but python datetime doesn't,
    # handle it separately
    datecmdline = "date -u '+%N'"
    ns_str = _rem_exec(cmdline=datecmdline)
    # take ns str, convert it to decimal seconds, round to 6 sig dec
    # and convert back to int microseconds
    us = int(round(float('.'+ns_str), 6)*10**6)
    lcunow = lcunow.replace(microsecond=us)
    return lcunow


if __name__ == '__main__':
    import sys
    url = sys.argv[1]
    cmdline = ' '.join(sys.argv[2:])
    print("Executing {} on {}".format(cmdline, url))
    try:
        hn = _exec_ssh(url, cmdline, nodetype='TEST', accessible=True, verbose=True)
        print(hn)
    except KeyboardInterrupt:
        print('INTERRUPTED')
