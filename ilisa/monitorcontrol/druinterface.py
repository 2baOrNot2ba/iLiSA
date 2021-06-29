import os
import plumbum
from ilisa.monitorcontrol._rem_exec import _exec_ssh
from ilisa.monitorcontrol.modeparms import band2rcumode, timestr2datetime,\
    normalizetimestr
from ilisa.pipelines.bfbackend import pl_rec_wrapper, dumpername

# dumpername is name of binary executable on DRU which is run by
# pl_rec_wrapper when capturing UDP packets with LOFAR beamformed voltages data.
# (Currently requires manually install putting it in DRU user's PATH env var)


class DRUinterface:
    druprompt = "On DRU>"
    verbose = True

    def __init__(self, accessconf_dru, ports=None):
        self.hostname = accessconf_dru.get('hostname', 'localhost')
        self.user = accessconf_dru.get('user', None)
        if not self.user and self.hostname == 'localhost':
            self.user = os.getlogin()
        if self.hostname != 'localhost':
            dru = plumbum.SshMachine(self.hostname, user=self.user)
        else:
            dru = plumbum.local
        self.accessible = True
        dru = self._exec_dru_func()
        self.dru = dru
        self.ports = ports

    def _exec_dru_func(self):
        nodeurl = '{}@{}'.format(self.user, self.hostname)
        def exec_ssh_inner(cmdline, nodeurl=nodeurl, stdoutdir=None,
                           nodetype='DRU', background_job=False, dryrun=False,
                           accessible=self.accessible, quotes="'", verbose=self.verbose):
            return _exec_ssh(nodeurl, cmdline, stdoutdir=stdoutdir,
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
        rcumode = band2rcumode(band)
        outarg = os.path.join(outdumpdir, outfilepre)
        dumplogname = '{}_lane{}_rcu{}.log'.format(dumpername, lane,
                                                   rcumode)
        dumplogpath = os.path.join(outdumpdir, dumplogname)
        # local_hostname = self.dru['hostname']().rstrip()
        local_hostname = self.dru('hostname').rstrip()
        starttime = normalizetimestr(starttime)
        starttime_arg = starttime + '.000'
        datapathguess = outarg + '_' + str(port) + '.' + local_hostname + '.'\
                        + starttime_arg
        if compress:
            datapathguess += '.zst'
        return outdumpdir, outarg, datapathguess, dumplogpath

    def _rec_bf_proxy(self, ports, duration, bf_data_dir, starttime='NOW',
                      compress=True, band='110_190', stnid=None,
                      mockrun=False):
        """\
        Record beamformed streams using recording process on DRU

        Note: Blocks until finished recording on DRU
        """
        dumpercmd = pl_rec_wrapper
        which_recorder = 'ow'
        outdumpdirs, outargs, datafiles, logfiles = \
            self.bfsfilepathslist(starttime, band, bf_data_dir, ports, stnid,
                                  compress)
        for outdumpdir in outdumpdirs:
            self.dru('mkdir -p '+outdumpdir)
        portlststr = ','.join([str(p) for p in ports])
        rcumode = band2rcumode(band)
        cmdlineargs = ['--which', which_recorder, '--ports', portlststr,
                       '--duration', str(duration),
                       '--bfdatadir', '"'+bf_data_dir+'"', '--rcumode', rcumode,
                       '--stnid', stnid]
        cmdlineargs.extend(['--starttime', starttime])
        if compress:
            cmdlineargs.append('--compress')
        if mockrun:
            cmdlineargs.append('--mockrun')
        self.dru(' '.join([dumpercmd] + cmdlineargs))
        return datafiles, logfiles


if __name__ == "__main__":
    import sys
    from ilisa.monitorcontrol.scansession import get_proj_stn_access_conf
    stnid = sys.argv.pop()
    projid = sys.argv.pop()
    accessconf = get_proj_stn_access_conf(projid, stnid)
    dru_interface = DRUinterface(accessconf['DRU'])
    ports = [0,1,2,3]
    duration_tot = 1.0
    scanpath_bfdat = accessconf['DRU']['BeamFormDataDir']
    rectime = 'ASAP'
    band = '110_190'
    _datafiles, _logfiles = dru_interface._rec_bf_proxy(ports,
            duration_tot, scanpath_bfdat, starttime=rectime, band=band,
            stnid=stnid, mockrun=True)
