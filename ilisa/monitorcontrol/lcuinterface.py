"""LOFAR module for basic hi-level functionality for observing in LOCAL station
mode.
"""

import time
import os

import ilisa.monitorcontrol
from ilisa.monitorcontrol._rem_exec import _exec_ssh
from ilisa.monitorcontrol.modeparms import parse_lofar_conf_files
# LOFAR constants
from ilisa.monitorcontrol.modeparms import rcumode2band, beamctl_args2cmds,\
    rcusetup_args2cmds, rspctl_stats_args2cmds


class LCUInterface(object):
    """This class provides an interface to the Local Control Unit (LCU) of an
       International LOFAR station."""
    lofarroot = "/opt/lofar_local/"  # "/opt/lofar/"
    lofarbin = "/opt/lofar/bin"
    lofaroperationsbin = "/opt/operations/bin"
    lofarstationtestdir = "/localhome/stationtest/"
    CalServer_conf = lofarroot + "etc/CalServer.conf"
    RSPDriver_conf = lofarroot + "etc/RSPDriver.conf"

    # User dirs & files:
    ## Directory containing user homes dirs
    homespath = "/data/home/"

    def checkLCUenv(self):
        """Check the LCU environment, especially for data taking assumptions."""
        if self.DryRun:
            return True, True
        envpath = self._stdoutLCU("env | grep '^PATH'")
        envpath = envpath.split("=")[1].split(":")
        if self.lofarbin in envpath and self.lofaroperationsbin in envpath:
            path_ok = True
        else:
            path_ok = False
        try:
            self.getdatalist()
        except:
            datadirs_ok = False
        else:
            datadirs_ok = True
        return path_ok, datadirs_ok

    def checkaccess(self):
        """Check that this object has access to LCU.
        Note: It does this by trying to get the MAC version number."""
        mac_version = ""
        try:
            mac_version = self.get_mac_version()
        except Exception:
            print("Error: Station object cannot access station with URL: "
                  + self.lcuURL)
        if mac_version is not "":
            accessible = True
        else:
            accessible = False
        return accessible

    def __init__(self, lcuaccessconf=None):
        if lcuaccessconf is None:
            accessconf = ilisa.monitorcontrol.default_access_lclstn_conf()
            lcuaccessconf = accessconf['LCU']

        # Initialize with user-station configuration:
        self.stnid = lcuaccessconf['stnid']
        self.user = lcuaccessconf['user']
        self.hostname = lcuaccessconf['hostname']
        self.lcuURL = self.user + "@" + self.hostname
        self.DryRun = False  # DryRun means commands to LCU are not executed
        self.verbose = True  # Write out LCU commands
        self.accessible = self.checkaccess()
        if self.accessible and self.verbose:
            print("Established access to LCU on {}.".format(self.stnid))

        # Init some OS paths:
        ## User's home dir
        self._home_dir = self.homespath + self.user + '/'
        ## Cache in user home (equiv. to XDG_CACHE_HOME)
        self.cache_dir = self._home_dir + '.cache/ilisa/'
        # This where the statistics data goes:
        self.lcuDumpDir = self.cache_dir + 'BSX_data/' # lcuaccessconf['dumpdir']
        ## ACCsrcDir is set in CalServer.conf and used when CalServer is running
        self.ACCsrcDir = self.cache_dir  + "ACC_data/"

        # Check LCU OS env:
        path_ok, datadirs_ok = self.checkLCUenv()
        assert path_ok, "Check $PATH on LCU. Needs to have {} and {}." \
            .format(self.lofarbin, self.lofaroperationsbin)
        assert datadirs_ok, "Check that data folders {} and {} on LCU." \
            .format(self.lcuDumpDir, self.ACCsrcDir)

    def __del__(self):
        pass

    def exec_lcu(self, cmdline, backgroundJOB=False, quotes="'"):
        return _exec_ssh(nodeurl=self.lcuURL, cmdline=cmdline,
                         stdoutdir=self.lcuDumpDir, nodetype='LCU',
                         background_job=backgroundJOB, dryrun=self.DryRun,
                         accessible=True, quotes=quotes, verbose=True)

    def _stdoutLCU(self, cmdline):
        return _exec_ssh(nodeurl=self.lcuURL, cmdline=cmdline,
                         stdoutdir=self.lcuDumpDir, nodetype='LCU',
                         background_job=False, dryrun=self.DryRun,
                         accessible=True, quotes="'", verbose=True)

    def _list_dat_files(self, dumpdir):
        ddls = self._stdoutLCU("ls " + dumpdir).split('\n')
        datfiles = []
        for ddfile in ddls:
            if ddfile.endswith('.dat'):
                datfiles.append(ddfile)
        return datfiles

    def getdatalist(self):
        """\
        Get list of data recorded on LCU. Returns a separate list for data
        dumps and ACC data. Note that because data dumps are named by date and
        time (uses YYYYmmdd_HHMMSS) and because default 'ls' sorts names in
        alphabetical order, the list of files will be ordered such that oldest
        data is first.
        """
        dryrun = self.DryRun
        self.DryRun = False  # Override DryRun for this case
        ls_lcuDumpDir = self._list_dat_files(self.lcuDumpDir)
        ls_ACCsrcDir = self._list_dat_files(self.ACCsrcDir)
        self.DryRun = dryrun
        return ls_lcuDumpDir, ls_ACCsrcDir

    def who_servicebroker(self):
        """Check who is running the Service Broker on the LCU. This is an
        indication of who is currently using the station."""
        if self.DryRun:
            return None
        try:
            ps_out = \
              self._stdoutLCU("/bin/ps -CServiceBroker --no-headers -ouser")
        except:
            # Can happen if no ServiceBroker process running. Pretend ps
            # returns blank
            ps_out = ""
        ps_out_lns = ps_out.splitlines()
        if len(ps_out_lns) == 0:
            sb_user = None
        else:
            sb_user = ps_out_lns[0]
        return sb_user

    def getstationswitchmode(self):
        """Get mode of station switch. Can be 'ILT' or 'local' mode,
        or `None` if undefined."""
        if self.DryRun:
            return 'local'
        try:
            getstationmode_out = self._stdoutLCU("getstationmode")
        except:
            # Can happen if no ServiceBroker process running
            print("Warning: ServiceBroker not running")
            stationmode = None
        else:
            stationmode = getstationmode_out.split()[-1]
        return stationmode

    def is_beam_on(self):
        """Check if there is a beamctl running"""
        try:
            ps_beamctl_outs = self._stdoutLCU("/bin/ps -Cbeamctl").splitlines()
        except:
            return False
        if len(ps_beamctl_outs) > 1:
            return True
        return False

    def _get_filecontent(self, filename):
        """Get contents of certain file on LCU by name."""
        filecontents = self._stdoutLCU("cat " + filename)
        return filecontents

    def _rm(self, source, override_mock=True):
        """Remove specified source file(s) on LCU."""
        dryrun = self.DryRun
        if override_mock:
            self.DryRun = False
        self.exec_lcu("rm -f " + source)
        self.DryRun = dryrun

    def cleanup(self):
        """Clean up all local data dumps. It is usually the ObsSession()
        objects responsibility to call this."""
        self._rm(self.lcuDumpDir + "/*.dat")
        self._rm(self.ACCsrcDir + "/*.dat")

    def get_RSPDriver_conf(self):
        """Get RSPDriver configuration settings."""
        filecontents = self._get_filecontent(self.RSPDriver_conf)
        rspdriver = parse_lofar_conf_files(filecontents)
        return rspdriver

    def rcu_disable_list(self, rcumode):
        """Return list of RCUs to be disabled.
         Looks for list in files on LCU (determined by ASTRON).
        """
        filename = (self.lofarstationtestdir + "DISABLED/disabled-mode"
                    + str(rcumode) + ".txt")
        filecontents = self._get_filecontent(filename)
        disabledrcus = filecontents.split(',')
        return disabledrcus

    def get_mac_version(self):
        """Get MAC version of station."""
        macversionstr = self._stdoutLCU("swlevel -V")
        if self.DryRun:
            macversionstr = "Mock-version-2.0.0"
        return macversionstr

    def get_swlevel(self):
        """Get current Software Level of station. Returns a string which
        under normal local mode operations is an integer between 0-3."""
        if not self.DryRun:
            swlevel = int(self._stdoutLCU("swlevel -S"))
        else:
            swlevel = 3
        return swlevel

    def set_swlevel(self, swleveltarget=3, fullreboot=False):
        """Set station's software level. swleveltarget=3 is the swlevel for which
        most monitorcontrol take place."""
        swlevel_changed = False
        if not self.DryRun:
            if not fullreboot:
                swlevel = self.get_swlevel()
                if swlevel != swleveltarget:
                    self.exec_lcu("swlevel {}".format(swleveltarget))
                    swlevel_changed = True
            else:
                # For completeness swlevel 0, but swlevel 1 is faster
                self.exec_lcu("swlevel 0; swlevel {}".format(swleveltarget))
                swlevel_changed = True
        return swlevel_changed

    def stop_beam(self):
        """Stop any running beamctl processes."""
        # Stop any beamctl on lcu.
        if not self.DryRun:
            self.exec_lcu("killall beamctl")
            # Put caltables back to default
            self.selectCalTable('default')
        # print("Beam off at %s"%time.asctime(time.localtime(time.time())))

    def run_rspctl(self, select=None, mode=None, tbbmode=None):
        """Run rspctl command to setup RCUs: rcumode, select.
        """
        # Form commandline argument string ignoring blank arguments:
        rspctl_cmd = "rspctl"
        if select is not None:
            rspctl_cmd += " --select={}".format(select)
        if mode is not None:
            rspctl_cmd += " --mode={}".format(mode)
        if tbbmode is not None:
            rspctl_cmd += " --tbbmode={}".format(tbbmode)
        self.exec_lcu(rspctl_cmd)
        return rspctl_cmd

    def rcusetup(self, bits, attenuation, mode=None):
        """Setup basic RCU setting: bits is 8 or 16, and attenuation is 0 to 31
        (0 means no attenuation & increasing number means more attenutation)"""
        rcusetup_cmds = rcusetup_args2cmds(bits, attenuation, mode=mode)
        # NOTE Looks like bitmode and rcuattenuation have to be set in separate
        #      commands.
        for rcusetup_cmd in rcusetup_cmds:
            self.exec_lcu(rcusetup_cmd)
        if self.DryRun:
            self.bits = bits
        waittime = 1
        print("Waiting {}s for rspctl to settle...".format(waittime))
        time.sleep(waittime)  # Wait for rspctl to settle
        return rcusetup_cmds

    def get_bits(self):
        """Get rcu sample bit-depth."""
        bits = getattr(self, 'bits', None)
        if not bits:
            outrspctl = self._stdoutLCU("rspctl --bitmode").splitlines()
            ans = outrspctl[2:]
            rspbits = []
            for l in ans:
                rsplnnr, head1684,  bitsstr = l.split(':', 2)
                rsp_nr = int(rsplnnr[4:6])
                rspbits.append((rsp_nr, int(bitsstr.lstrip())))
            # Will assume all rsp board have same bit depth:
            bits = rspbits[0][1]
        return bits

    def run_beamctl(self, beamlets, subbands, band, anadigdir, rcus='all',
                    beamdurstr='', backgroundJOB=True):
        """Start a beam using beamctl command. Blocks until ready."""
        beamctl_cmd = beamctl_args2cmds(beamlets, subbands, band, anadigdir,
                                        rcus, beamdurstr)
        self.exec_lcu(beamctl_cmd, backgroundJOB)
        waittime = 11
        print("Waiting {}s for beam to settle...".format(waittime))
        time.sleep(waittime)  # Wait for beam to settle
        return beamctl_cmd

    def run_rspctl_statistics(self, bsxtype, integration, duration, subband=0,
                              directory=None):
        """\
        Run rspctl statistics command
        """
        rspctl_cmds = rspctl_stats_args2cmds(bsxtype, integration, duration,
                                             subband=subband)
        if directory is None:
            directory = self.lcuDumpDir
        rspctl_cmds[-1] += " --directory={}".format(directory)
        for rspctl_cmd in rspctl_cmds:
            self.exec_lcu(rspctl_cmd)
        if self.DryRun:
            self.mockstatistics(bsxtype, integration, duration, directory)
        return rspctl_cmds

    def mockstatistics(self, statistics, integration, duration, directory=None):
        """Make mock statistics data file(s)."""
        if directory is None:
            if statistics != 'acc':
                directory = self.lcuDumpDir
            else:
                directory = self.ACCsrcDir
        dryrun = self.DryRun
        self.DryRun = False
        nrtimsamps = int(duration/integration)
        import ilisa.monitorcontrol.modeparms as modeparms
        import datetime
        dd_cmdbase = 'dd if=/dev/zero'
        now = datetime.datetime.utcnow()
        dtstamp = now.strftime('%Y%m%d_%H%M%S')
        if statistics == 'bst':
            bits = self.get_bits()
            # Write mock bst files. (mock files contain just zeros)
            nrbls = modeparms.NRBEAMLETSBYBITS[bits]
            dd_bs = 8
            dd_count = nrbls * nrtimsamps
            dd_cmdbase += ' bs={} count={}'.format(dd_bs, dd_count)
            for pol in ['X', 'Y']:
                bstfilename = "{}_bst_00{}.dat".format(dtstamp, pol)
                fpath = os.path.join(directory, bstfilename)
                dd_cmd = dd_cmdbase + ' of={}'.format(fpath)
                self.exec_lcu(dd_cmd)
        elif statistics == 'sst':
            # Write mock sst files
            nrsbs = modeparms.TotNrOfsb
            dd_bs = 8
            dd_count = nrsbs * nrtimsamps
            dd_cmdbase += ' bs={} count={}'.format(dd_bs, dd_count)
            for rcunr in range(modeparms.nrofrcus):
                sstfilename = "{}_sst_rcu{:03}.dat".format(dtstamp, rcunr)
                fpath = os.path.join(directory, sstfilename)
                dd_cmd = dd_cmdbase + ' of={}'.format(fpath)
                self.exec_lcu(dd_cmd)
        elif statistics == 'xst':
            # Write mock sst files
            nrrcus = modeparms.nrofrcus
            dd_bs = 2*8
            dd_count = nrrcus*nrrcus * nrtimsamps
            dd_cmdbase += ' bs={} count={}'.format(dd_bs, dd_count)
            xstfilename = "{}_xst.dat".format(dtstamp)
            fpath = os.path.join(directory, xstfilename)
            dd_cmd = dd_cmdbase + ' of={}'.format(fpath)
            self.exec_lcu(dd_cmd)
        elif statistics == 'acc':
            # Write mock acc files
            integration = 1
            nrrcus = modeparms.nrofrcus
            dd_bs = 2*8
            dd_count = nrrcus*nrrcus * 512
            dd_cmdbase += ' bs={} count={}'.format(dd_bs, dd_count)
            nraccfiles, restsamp = divmod(duration+7, 519)
            filetime = now
            for accnr in range(nraccfiles):
                dtstamp = filetime.strftime('%Y%m%d_%H%M%S')
                accfilename = "{}_acc_512x196x196.dat".format(dtstamp)
                fpath = os.path.join(directory, accfilename)
                dd_cmd = dd_cmdbase + ' of={}'.format(fpath)
                self.exec_lcu(dd_cmd)
                filetime += datetime.timedelta(seconds=519)
        # Wait duration seconds (disregard time code above takes)
        time.sleep(duration)
        self.DryRun = dryrun

    def run_tbbctl(self, select=None, alloc=False, free=False, record=False,
                   stop=False, mode=None, storage=None, readall=None,
                   cepdelay=None, backgroundJOB=False):
        """Run the tbbctl command on the LCU with arguments given."""
        tbbctl_args = ""

        if alloc:
            tbbctl_args += " --alloc"
        if free:
            tbbctl_args += " --free"
        if record:
            tbbctl_args += " --record"
        if stop:
            tbbctl_args += " --stop"
        if mode is not None:
            tbbctl_args += " --mode={}".format(mode)
        if storage is not None:
            tbbctl_args += " --storage={}".format(storage)
        if readall is not None:
            tbbctl_args += " --readall={}".format(readall)
        if cepdelay is not None:
            tbbctl_args += " --cepdelay={}".format(cepdelay)
        if select is not None:
            # `select` needs to come last
            tbbctl_args += " --select={}".format(select)
        if tbbctl_args != "":
            tbbctl_cmd = "tbbctl"+tbbctl_args
            self.exec_lcu(tbbctl_cmd, backgroundJOB)
        else:
            tbbctl_cmd = None
        return tbbctl_cmd

    def acc_mode(self, enable=True):
        """Enable or disable ACC mode.
        If enableacc=True, ACC files will be written to file when CalServer is running
        (swlevel>=2). If enableacc=False, ACC files will not be written.
        """
        if enable:
            self.exec_lcu(
            r"""sed -i.orig 's/^CalServer.DisableACMProxy=1/CalServer.DisableACMProxy=0/;\
            s/^CalServer.WriteACCToFile=0/CalServer.WriteACCToFile=1/;\
            s,^CalServer.DataDirectory=.*,CalServer.DataDirectory={}, ' {}"""
                .format(self.ACCsrcDir, self.CalServer_conf), quotes='"')
        else:
            self.exec_lcu(
            r"""sed -i 's/^CalServer.DisableACMProxy=0/CalServer.DisableACMProxy=1/;\
            s/^CalServer.WriteACCToFile=1/CalServer.WriteACCToFile=0/;\
            s,^CalServer.DataDirectory=.*,CalServer.DataDirectory=/localhome/data,' {}"""
                .format(self.CalServer_conf), quotes='"')

    def getCalTableInfo(self, rcumode):
        """Fetch and return the caltable info from the LCU."""
        calinfo = None
        if self.DryRun:
            return ""
        if int(rcumode) == 4:
            # Band 30_90 not correctly implemented in "beamctl --calinfo".
            # It uses the 10_90 caltab anyways so:
            rcumode = 3
        calinfoout =self._stdoutLCU("beamctl --calinfo")
        if self.DryRun:
            return ""
        # Convert output into a list of dict per antset
        calinfolist = []
        # Strip off first initial lines and split on blank lines
        calinfooutlist = (''.join(calinfoout.splitlines(True)[2:])).split('\n\n')
        if calinfooutlist[0] == '':
            # No calinfo. Return blank
            return ""
        for calinfooutlistitem in calinfooutlist:
            calinfolistitem = {}
            for calinfolistitemline in calinfooutlistitem.split('\n'):
                calkey, calval = calinfolistitemline.split(':', 1)
                calkey = calkey.rstrip()
                calval = calval.lstrip()
                calinfolistitem[calkey] = calval
            calinfolist.append(calinfolistitem)
        band = rcumode2band(rcumode)
        for calinfolistitem in calinfolist:
            if calinfolistitem['Band'] == band:
                calinfo = calinfolistitem
                break
        return calinfo

    def selectCalTable(self, which='default'):
        """This is specific to an lcu which has the script SelectCalTable.sh
        with which a user can switch between different caltables.
        """
        pass

    def turnoffLBA_LNAs(self, select="0:191"):
        """Turn-off the LNAs on LBA.
        (Used as an indication of system temperature)"""
        rspctl_cmd = "rspctl --rcu=0x00034880 --select=" + select
        self.exec_lcu(rspctl_cmd)
        if self.verbose:
            print("Turning OFF LBA LNAs...")
        time.sleep(30)
        return rspctl_cmd

    # **SEPTON
    setElem_ON = 128
    setElem_OFF = 2
    elementsInTile = 16
    nrTiles = 96

    def _tiles2rcus(self, tiles):
        rcus = []
        for tile in tiles:
            rcus.extend([2 * tile, 2 * tile + 1])  # Set same delay for both X&Y pol rcu
        return rcus

    def turnoffElinTile_byTile(self, elemsOn):
        """"Turn off all elements per tile except the one specificied in list.
        Execution is done by tile, which is more intuitive but slower."""
        self.run_rspctl(mode='5')
        for tileNr in range(self.nrTiles):
            # Start with all elements in tile off
            # (2 is OFF)
            tileMap = [self.setElem_OFF for elemNr in range(self.elementsInTile)]
            # Turn on the appropriate element
            tileMap[elemsOn[tileNr]] = self.setElem_ON  # 128 is ON
            lcucmd = "rspctl --hbadelay="\
                     + str(tileMap).strip('[]').replace(" ", "")\
                     + " --select="+str(2*tileNr)+","+str(2*tileNr+1)
            self.exec_lcu(lcucmd)

    def turnoffElinTile_byEl(self, elems_on):
        """"Turn off all elements per tile except the one specificied in list.
        Execution is done by element, which is less intuitive but faster."""
        self.run_rspctl(mode='5')
        for elNr in range(self.elementsInTile):
            tiles = [ind for ind in range(self.nrTiles) if elNr == elems_on[ind]]
            if len(tiles) == 0:
                continue
            tile_map = [self.setElem_OFF for elemNr in range(self.elementsInTile)]
            tile_map[elNr] = self.setElem_ON
            rcus = self._tiles2rcus(tiles)
            lcucmd = "rspctl --hbadelay="\
                     + str(tile_map).strip('[]').replace(" ", "")\
                     + " --select="+str(rcus).strip('[]').replace(" ", "")
            self.exec_lcu(lcucmd)
