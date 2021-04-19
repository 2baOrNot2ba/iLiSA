import os
import plumbum
from ilisa.monitorcontrol._rem_exec import _exec_ssh
import ilisa.monitorcontrol.modeparms

# Name of binary executable on DRU to run when capturing
# UDP packets with LOFAR beamformed voltages data.
# (Currently requires manually install putting it in DRU user's PATH env var.
# )
dumpername = 'pl_rec'


class DRUinterface:
    druprompt = "On DRU>"
    verbose = True

    def __init__(self, accessconf_dru, ports=None):
        try:
            hostname = accessconf_dru['hostname']
        except KeyError:
            hostname = 'localhost'
        self.hostname = hostname
        try:
            user = accessconf_dru['user']
        except KeyError:
            user = None
        if hostname != 'localhost':
            dru = plumbum.SshMachine(self.hostname, user=user)
        else:
            dru = plumbum.local
        self.accessible = True
        dru = self._exec_dru_func()
        self.dru = dru
        self.ports = ports

    def _exec_dru_func(self):
        def exec_ssh_inner(cmdline, remnode=self.hostname, stdoutdir=None,
                           nodetype='DRU', background_job=False, dryrun=False,
                           accessible=self.accessible, quotes="'", verbose=self.verbose):
            return _exec_ssh(remnode, cmdline, stdoutdir=stdoutdir,
                             nodetype=nodetype, background_job=background_job,
                             dryrun=dryrun, accessible=accessible,
                             quotes=quotes, verbose=verbose)
        return exec_ssh_inner

    def bfsfilepathslist(self, starttime, band, bf_data_dir, ports, stnid,
                         compress=True):
        port0 = ports[0]
        outdumpdirs = []
        outargs = []
        datafileguesses = []
        dumplognames = []
        for lane in range(len(ports)):
            outdumpdir, outarg, datafileguess, dumplogname = \
                self.bfsfilepaths(lane, starttime, band, bf_data_dir, port0,
                                  stnid, compress=compress)
            outdumpdirs.append(outdumpdir)
            outargs.append(outarg)
            datafileguesses.append(datafileguess)
            dumplognames.append(dumplogname)
        return outdumpdirs, outargs, datafileguesses, dumplognames

    def bfsfilepaths(self, lane, starttime, band, bf_data_dir, port0, stnid,
                     compress=True):
        """Generate paths and name for BFS recording.

        Parameters
        ----------
        lane : int
            Lane number 0,1,2, or 3.
        starttime : str
            The datetime string when the BF stream started.
        band :
            The band name for the BF stream.
        bf_data_dir : str
            Template for BF data lane dump directory. Should have format:
                <pre_bf_dir>?<pst_bf_dir>
            where '?' will be replaced by the lane number.
        port0 :
            The port number of lane 0.
        stnid :
            Station ID.
        compress: bool
            Whether or not compression is used.

        Returns
        -------
        outdumpdir : str
            Directory in which to dump bf data. Has format:
                <outdumpdir>/udp_<stnid>
            where
                <outdumpdir> := <rootdir>/lane<lanenr>/<pst_bf_dir>
        outarg : str
            Argument 'out' passed to dumper CLI.
        datafileguess : str
            Path to data file. Has format:
                _<port>.start.<%Y-%m-%dT%H:%M:%S>.000
        dumplogname :
            Name of dumper's logfile.
        """
        port = port0 + lane
        pre_bf_dir, pst_bf_dir = bf_data_dir.split('?')
        outdumpdir = pre_bf_dir + str(lane) + pst_bf_dir
        outfilepre = "udp_" + stnid
        rcumode = ilisa.monitorcontrol.modeparms.band2rcumode(band)
        outarg = os.path.join(outdumpdir, outfilepre)
        dumplogname = '{}_lane{}_rcu{}.log'.format(dumpername, lane, rcumode)
        # local_hostname = self.dru['hostname']().rstrip()
        local_hostname = self.dru('hostname').rstrip()
        starttime_arg = starttime + '.000'
        datafileguess = outarg + '_' + str(port) + '.' + local_hostname + '.'\
                        + starttime_arg
        if compress:
            datafileguess += '.zst'
        return outdumpdir, outarg, datafileguess, dumplogname

    def _rec_bf_proxy(self, ports, duration, bf_data_dir, starttime='NOW',
                      compress=False, band=None, stnid=None):
        """\
        Record beamformed streams using recording process on DRU

        Note: Blocks until finished recording on DRU
        """
        dumpercmd = dumpername
        startarg = starttime
        if starttime != 'NOW':
            startarg = starttime.strftime("%Y-%m-%dT%H:%M:%S")
        outdumpdirs, outargs, datafiles, logfiles = \
            self.bfsfilepathslist(startarg, band, bf_data_dir, ports, stnid,
                                  compress)
        for outdumpdir in outdumpdirs:
            self.dru('mkdir -p '+outdumpdir)
        portlststr = ','.join([str(p) for p in ports])
        cmdlineargs = ['--ports', portlststr, '--duration', str(duration)]
        if startarg != 'NOW':
            cmdlineargs.extend(['--Start', startarg])
        if compress:
            pass
            # cmdlineargs.append('--compress')
        self.dru(' '.join([dumpercmd] + cmdlineargs))
        return datafiles, logfiles

    def rec_bf_proxy(self, starttime, duration, lanes, band, bf_data_dir,
                     port0, stnid, compress=False):
        """Start recording beamformed streams using an external dumper process.
        """
        recorders = ['ow', 'py']
        _which_recorder = recorders[0]
        # dumpercmd = self.dru.path(self.pipeline_path) / dumpername
        dumpercmd = dumpername
        rec_cmd = self.dru[dumpercmd]
        startarg = starttime
        if starttime != 'NOW':
            startarg = starttime.strftime("%Y-%m-%dT%H:%M:%S")
        reclanes = []
        datafiles = []
        logfiles = []
        for lane in lanes:
            port = port0 + lane
            outdumpdir, outarg, datafileguess, dumplogname = \
                self.bfsfilepaths(lane, startarg, band, bf_data_dir, port0,
                                  stnid, compress)
            self.dru['mkdir']['-p'](outdumpdir)
            cmdlineargs = ['--ports', port, '--check', '--duration', duration,
                           '--timeout', '999', '--out',  outarg]
            if startarg != 'NOW':
                cmdlineargs.extend(['--Start', startarg])
            if compress:
                cmdlineargs.append('--compress')
            if self.verbose:
                cmdlineargstr = ' '.join(map(str, cmdlineargs))
                print("{} {} {}".format(self.druprompt, dumpercmd,
                                        cmdlineargstr))
            reclanes.append((rec_cmd[cmdlineargs] > dumplogname) & plumbum.BG)
            datafiles.append(datafileguess)
            logfiles.append(dumplogname)
        for reclane in reclanes:
            reclane.wait()
        return datafiles, logfiles


if __name__ == "__main__":
    import sys
    from ilisa.monitorcontrol.scansession import get_proj_stn_access_conf
    stnid = sys.argv.pop()
    projid = sys.argv.pop()
    accessconf = get_proj_stn_access_conf(projid, stnid)
    dru_interface = DRUinterface(accessconf['DRU'])
    #dru_interface._rec_bf_proxy([16070, 16071], 5, '/mnt/lane?/BF/SE607/Scans/',
    #                            starttime='NOW', compress=False, band='110_190', stnid='SE607')
