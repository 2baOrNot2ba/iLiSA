import os
import datetime

from ilisa.operations._rem_exec import _exec_ssh
from ilisa.operations.modeparms import band2rcumode, normalizetimestr
from ilisa.pipelines.bfbackend import PL_REC_WRAPPER, DUMPERNAME


class DRUinterface:
    druprompt = "On DRU>"
    verbose = True

    def __init__(self, accessconf_dru, ports=None):
        self.hostname = accessconf_dru.get('hostname', 'localhost')
        self.user = accessconf_dru.get('user', None)
        if not self.user and self.hostname == 'localhost':
            self.user = os.getlogin()
        self.url = self.user + "@" + self.hostname
        self.bf_logfile = accessconf_dru.get('BeamFormLogFile')
        self.bf_data_dir = accessconf_dru.get('BeamFormDataDir')
        self.accessible = True
        dru = self._exec_dru_func()
        self.dru = dru
        self.ports = ports

    def _exec_dru_func(self, background_job=False, stdoutdir='~'):
        nodeurl = '{}@{}'.format(self.user, self.hostname)

        def exec_ssh_inner(cmdline, nodeurl=nodeurl, nodetype='DRU',
                           background_job=background_job,  stdoutdir=stdoutdir,
                           dryrun=False, accessible=self.accessible, quotes="'",
                           verbose=self.verbose):
            return _exec_ssh(nodeurl, cmdline, nodetype=nodetype,
                             background_job=background_job, stdoutdir=stdoutdir,
                             dryrun=dryrun, accessible=accessible,
                             quotes=quotes, verbose=verbose)
        return exec_ssh_inner

    def ls(self, path):
        """List files in path on DRU"""
        ls_out = self.dru('ls ' + path)
        if ls_out:
            ls_outs = ls_out.split('\n')
        else:
            ls_outs = []
        return ls_outs

    def bfsfilepathslist(self, starttime, band, bf_data_dir, ports, stnid,
                         compress=True):
        """\
        Generate paths and name for BFS recording

        Parameters
        ----------
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
        def _bfsfilepaths(lane, starttime, band, bf_data_dir, port0, stnid,
                          compress=True):
            """\
            Generate paths and name for a BFS recording lane
            """
            port = port0 + lane
            pre_bf_dir, pst_bf_dir = bf_data_dir.split('?')
            outdumpdir = pre_bf_dir + str(lane) + pst_bf_dir
            outfilepre = "udp_" + stnid
            rcumode = band2rcumode(band)
            outarg = os.path.join(outdumpdir, outfilepre)
            dumplogname = '{}_lane{}_rcu{}.log'.format(DUMPERNAME, lane,
                                                       rcumode)
            dumplogpath = os.path.join(outdumpdir, dumplogname)
            # local_hostname = self.dru['hostname']().rstrip()
            local_hostname = self.dru('hostname').rstrip()
            starttime = normalizetimestr(starttime)
            starttime_arg = starttime + '.000'
            datapathguess = outarg + '_' + str(port) + '.' + local_hostname + '.' \
                            + starttime_arg
            if compress:
                datapathguess += '.zst'
            return outdumpdir, outarg, datapathguess, dumplogpath

        port0 = ports[0]
        outdumpdirs = []
        outargs = []
        datafileguesses = []
        dumplognames = []
        for lane in range(len(ports)):
            outdumpdir, outarg, datafileguess, dumplogname = \
                _bfsfilepaths(lane, starttime, band, bf_data_dir, port0, stnid,
                              compress=compress)
            outdumpdirs.append(outdumpdir)
            outargs.append(outarg)
            datafileguesses.append(datafileguess)
            dumplognames.append(dumplogname)
        return outdumpdirs, outargs, datafileguesses, dumplognames

    def parse_bfs_filename(self, bfs_filepath):
        """\
        Parse name of BFS file to extract basic observation metadata

        Parameters
        ----------
        bfs_filepath : str
            Path to BFS file.

        Returns
        -------
        port: int
            Port number
        hostname: str
            Name of station
        startstr: str
            Datetime of observation start in ISO format.
        ms: str
            Fractional seconds of start
        cmprss_suf: str
            Compression suffix. If equal to 'zst' BFS compresses with Zlib.
        """
        bfs_filename = os.path.basename(bfs_filepath)
        udp, stnid, prt_host_start_ms_compress = bfs_filename.split('_', 2)
        port, hostname, startstr, ms, cmprss_suf =\
            prt_host_start_ms_compress.split('.', 4)
        return port, hostname, startstr, ms, cmprss_suf

    def get_bfs_filenames(self, bf_dat_dir_tmplt):
        """\
        Get BFS data and log file names on DRU

        Parameters
        ----------
        bf_data_dir : str
            Data dir where bfs files were recorded

        Returns
        -------
        paths_data : list
            Paths to the BFS data files.
        paths_logs : list
            Paths to the BFS log files.
        """
        nrlanes = 4
        bf_dir_list = [bf_dat_dir_tmplt]
        # Convert /paths/?/like/this to real paths
        if '?' in bf_dat_dir_tmplt:
            bf_dir_list = [bf_dat_dir_tmplt.replace('?', str(lanenr))
                           for lanenr in range(nrlanes)]
        paths_data = []
        paths_logs = []
        for _bf_data_dir in bf_dir_list:
            ls = self.ls(_bf_data_dir)
            if ls:
                # Find all non .log files
                data_files = sorted(filter(lambda f: not f.endswith('.log'), ls))
                _paths_data = [os.path.join(_bf_data_dir, f) for f in
                               data_files]
                paths_data.extend(_paths_data)
                # Find all .log files
                log_files = sorted(filter(lambda f: f.endswith('.log'), ls))
                _paths_logs = [os.path.join(_bf_data_dir, f) for f in
                               log_files]
                paths_logs.extend(_paths_logs)
        return paths_data, paths_logs

    def start_bf_rec(self, ports, duration, bf_data_dir, starttime,
                     file_dur=None, compress=True, band='110_190', stnid=None,
                     mockrun=False):
        """\
        Start Record beamformed streams using recording process on DRU

        Returns
        -------
        datafiles : list
            The paths to the BFS data files.
        logfiles : list
            The paths to the BFS log files.
        """
        dumpercmd = PL_REC_WRAPPER
        which_recorder = 'ow'
        starttime_str = starttime.strftime('%Y-%m-%dT%H:%M:%S')
        outdumpdirs, outargs, datafiles, logfiles = \
            self.bfsfilepathslist(starttime_str, band, bf_data_dir, ports,
                                  stnid, compress)
        for outdumpdir in outdumpdirs:
            self.dru('mkdir -p '+outdumpdir)
        portlststr = ','.join([str(p) for p in ports])
        rcumode = band2rcumode(band)
        cmdlineargs = ['--which', which_recorder, '--ports', portlststr,
                       '--duration', str(duration),
                       '--bfdatadir', '"'+bf_data_dir+'"', '--rcumode', rcumode,
                       '--stnid', stnid]
        if file_dur and file_dur != duration:
            cmdlineargs += ['--file_duration', str(file_dur)]
        cmdlineargs.extend(['--starttime', starttime_str])
        if compress:
            cmdlineargs.append('--compress')
        if mockrun:
            cmdlineargs.append('--mockrun')
        print("BFS rec process issued to DRU @", datetime.datetime.utcnow())
        self.dru(' '.join([dumpercmd] + cmdlineargs), background_job=True)
        print("BFS_end")
        return datafiles, logfiles

    def _rec_bf_proxy(self, ports, duration, bf_data_dir, starttime,
                      file_dur=None, compress=True, band='110_190', stnid=None,
                      mockrun=False):
        """\
        Record beamformed streams using recording process on DRU

        Note: Blocks until finished recording on DRU
        """
        dumpercmd = PL_REC_WRAPPER
        which_recorder = 'ow'
        starttime_str = starttime.strftime('%Y-%m-%dT%H:%M:%S')
        outdumpdirs, outargs, datafiles, logfiles = \
            self.bfsfilepathslist(starttime_str, band, bf_data_dir, ports,
                                  stnid, compress)
        for outdumpdir in outdumpdirs:
            self.dru('mkdir -p '+outdumpdir)
        portlststr = ','.join([str(p) for p in ports])
        rcumode = band2rcumode(band)
        cmdlineargs = ['--which', which_recorder, '--ports', portlststr,
                       '--duration', str(duration),
                       '--bfdatadir', '"'+bf_data_dir+'"', '--rcumode', rcumode,
                       '--stnid', stnid]
        if file_dur:
            cmdlineargs += ['--file_duration', str(file_dur)]
        cmdlineargs.extend(['--starttime', starttime_str])
        if compress:
            cmdlineargs.append('--compress')
        if mockrun:
            cmdlineargs.append('--mockrun')
        print("BFS rec process issued to DRU @", datetime.datetime.utcnow())
        self.dru(' '.join([dumpercmd] + cmdlineargs))
        print("BFS_end")
        return datafiles, logfiles


if __name__ == "__main__":
    import sys
    from ilisa.operations.scansession import get_proj_stn_access_conf
    stnid = sys.argv.pop()
    projid = sys.argv.pop()
    accessconf = get_proj_stn_access_conf(projid, stnid)
    dru_interface = DRUinterface(accessconf['DRU'])
    ports = [4346, 4347, 4348, 4349]
    duration_tot = 1.0
    scanpath_bfdat = accessconf['DRU']['BeamFormDataDir']
    starttime = datetime.datetime.utcnow()
    band = '110_190'
    _datafiles, _logfiles = dru_interface._rec_bf_proxy(ports,
            duration_tot, scanpath_bfdat, starttime, band=band,
            stnid=stnid, mockrun=False)
