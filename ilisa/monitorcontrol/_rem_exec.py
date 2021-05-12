import subprocess
import os
try:
    import paramiko
    IMPORTED_PARAMIKO = True
except ImportError:
    IMPORTED_PARAMIKO = False

def _exec_rem(remnode, cmdline, stdoutdir, nodetype='LCU',
              background_job=False, dryrun=False, accessible=False, quotes="'",
              verbose=True):
    return _exec_ssh(remnode, cmdline, stdoutdir, nodetype=nodetype,
                     background_job=background_job, dryrun=dryrun,
                     accessible=accessible, quotes=quotes, verbose=verbose)


def _exec_ssh(nodeurl, cmdline, stdoutdir='~', nodetype='LCU',
              background_job=False, dryrun=False, accessible=False, quotes="'",
              verbose=True):
    """Execute a command on the remnode, either as a background job or in the
    foreground (blocking). Typically access is remote via ssh.
    (To speed things up use the ssh CommandMaster option.)
    """
    nodeprompt = "On " + nodetype + "> "
    if nodeurl.endswith('localhost'):
        shellinvoc = ''
        quotes = ''
    else:
        shellinvoc = "ssh " + nodeurl
    output = None
    if background_job:
        # Currently only run_beamctl & run_tbbctl run in background
        # Put stdout & stderr in log in dumpdir
        cmdline = ("(( " + cmdline + " ) > " + stdoutdir
                   + "lcu_shell_out.log 2>&1) &")
    if dryrun:
        pre_prompt = "(dryrun) "
    else:
        pre_prompt = ""
    if verbose:
        print(pre_prompt + nodeprompt + cmdline)
    if (not dryrun) and accessible:
        if background_job == 'locally':
            # Runs in background locally rather than in background on LCU
            output = subprocess.run(shellinvoc + " " + cmdline + " &",
                                    shell=True, stdout=subprocess.PIPE).stdout
        else:
            output = subprocess.run(shellinvoc + " "
                                    + quotes + cmdline + quotes,
                                    shell=True, universal_newlines = True,
                                    stdout=subprocess.PIPE).stdout
        if output:
            output = output.rstrip()
    elif not accessible:
        print("Warning: not running as " + nodeurl
              + " since it is not accesible.")
    return output


def __exec_lcu_paramiko(self, cmdline, backgroundJOB=False):
    lcuprompt = "LCUp>"
    if self.DryRun:
        preprompt = "(dryrun)"
    else:
        preprompt = ""
    if backgroundJOB is True:
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
        prePrompt = "(dryrun) "
    else:
        prePrompt = ""
    if verbose:
        print(prePrompt + nodeprompt + cmdline)
    if not(dryrun):
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
    shellinvoc = "ssh " + self.lcuURL
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


if __name__ == '__main__':
    hn = _exec_ssh('user6@se607c', 'hostname', accessible=True)
    print(hn)
