import subprocess


def _exec_rem(remnode, cmdline, stdoutdir, nodetype='LCU',
              background_job=False, dryrun=False, accessible=False, quotes="'",
              verbose=True):
    return _exec_ssh(remnode, cmdline, stdoutdir, nodetype=nodetype,
                     background_job=background_job, dryrun=dryrun,
                     accessible=accessible, quotes=quotes, verbose=verbose)


def _exec_ssh(nodeurl, cmdline, stdoutdir, nodetype='LCU',
              background_job=False, dryrun=False, accessible=False, quotes="'",
              verbose=True):
    """Execute a command on the remnode, either as a background job or in the
    foreground (blocking). Typically access is remote via ssh.
    (To speed things up use the ssh CommandMaster option.)
    """
    nodeprompt = "On " + nodetype + "> "
    if nodeurl == 'localhost':
        shellinvoc = ''
        quotes = ''
    else:
        shellinvoc = "ssh " + nodeurl
    output = None
    if background_job is True:
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
    elif not accessible:
        print("Warning: not running as " + nodeurl
              + " since it is not accesible.")
    return output


def _stdout_ssh(remnode, cmdline, nodetype='LCU', dryrun=False,
                verbose=True):
    """Execute a command on the remnode and return its output."""
    nodeprompt = "On " + nodetype + "> "
    shellinvoc = "ssh " + remnode
    if dryrun:
        prePrompt = "(dryrun) "
    else:
        prePrompt = ""
    if verbose:
        print(prePrompt + nodeprompt + cmdline)
    if dryrun is False:
        try:
            output = subprocess.check_output(shellinvoc + " '" + cmdline + "'",
                                             shell=True).rstrip()
            output = str(output.decode('UTF8'))
        except subprocess.CalledProcessError as e:
            raise e
    else:
        output = "None"
    return output
