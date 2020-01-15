"""LOFAR module for basic hi-level functionality for observing in LOCAL station
mode.
"""

import time
import subprocess
import os
try:
    import paramiko
    IMPORTED_PARAMIKO = True
except ImportError:
    IMPORTED_PARAMIKO = False

import ilisa.observations
from ilisa.observations.modeparms import parse_lofar_conf_files
# LOFAR constants
from ilisa.observations.modeparms import nrofrcus, band2antset, rcumode2band


class LCUInterface(object):
    """This class provides an interface to the Local Control Unit (LCU) of an
       International LOFAR station."""
    lofarroot = "/opt/lofar_local/"
    lofarbin = "/opt/lofar/bin"
    lofaroperationsbin = "/opt/operations/bin"
    lofarstationtestdir = "/localhome/stationtest/"
    # ACCsrcDir is set in CalServer.conf and used when CalServer is running.
    ACCsrcDir = "/localhome/data/ACCdata/"
    CalServer_conf = lofarroot + "/etc/CalServer.conf"
    RSPDriver_conf = lofarroot + "/etc/RSPDriver.conf"

    def setupaccess(self, accessconf):
        """Initialize with user-station configuration."""
        self.stnid = accessconf['stnid']
        self.user = accessconf['user']
        self.hostname = accessconf['hostname']
        self.lcuURL = self.user+"@"+self.hostname
        self.lcuHome = accessconf['home']
        # This where the statistics data goes:
        self.lcuDumpDir = accessconf['dumpdir']
        # Should lcu scripts be used?:
        self.usescriptonlcu = accessconf['usescriptonlcu']
        # TODO Implement condition: if self.usescriptonlcu:
        # This is where the scripts are:
        # TODO Remove dependency on scriptsDir:
        # (scripts should run on system-wide PATH)
        self.scriptsDir = "/data/home/"+accessconf['user']+"/scripts/"
        self.DryRun = accessconf['DryRun']
        self.verbose = True  # Write out LCU commands
        if self.checkaccess() and self.verbose:
            print("Established access to LCU.")
        pathOK, datadirsOK = self.checkLCUenv()
        assert pathOK, "Check $PATH on LCU. Needs to have {} and {}."\
            .format(self.lofarbin, self.lofaroperationsbin)
        assert datadirsOK, "Check that data folders {} and {} on LCU."\
            .format(self.lcuDumpDir, self.ACCsrcDir)

    def checkLCUenv(self):
        """Check the LCU environment, especially for datataking assumptions."""
        if self.DryRun:
            return True, True
        envpath = self._stdoutLCU("env | grep '^PATH'")
        envpath = envpath.split("=")[1].split(":")
        if self.lofarbin in envpath and self.lofaroperationsbin in envpath:
            pathOK = True
        else:
            pathOK = False
        try:
            ls_lcuDumpDir, ls_ACCsrcDir = self.getdatalist()
        except OSError as err:
            datadirsOK = False
        else:
            datadirsOK = True
        return pathOK, datadirsOK

    def checkaccess(self):
        """Check that this object has access to LCU.
        Note: It does this by trying to get the MAC version number."""
        self.MACversion = ""
        try:
            self.MACversion = self.get_mac_version()
        except Exception:
            print("Error: Station object cannot access station with URL: "
                  + self.lcuURL)
        if self.MACversion is not "":
            self.accessible = True
        else:
            self.accessible = False
        return self.accessible

    def __init__(self, lcuaccessconf=None):
        if lcuaccessconf is None:
            accessconf = ilisa.observations.default_access_lclstn_conf()
            lcuaccessconf = accessconf['LCU']
        self.setupaccess(lcuaccessconf)

    def __del__(self):
        pass

    def exec_lcu(self, cmdline, backgroundJOB=False, quotes="'", method='ssh'):

        if method == 'ssh':
            self._exec_lcu_ssh(cmdline, backgroundJOB=backgroundJOB, quotes=quotes)
        elif method == 'paramiko':
            self._exec_lcu_paramiko(cmdline, backgroundJOB=backgroundJOB)

    def _exec_lcu_paramiko(self, cmdline, backgroundJOB=False):
        lcuprompt = "LCUp>"
        if self.DryRun:
            preprompt = "(dryrun)"
        else:
            preprompt = ""
        if backgroundJOB is True:
            cmdline = "(( "+cmdline+" ) > "+self.lcuHome+"lofarctl.log 2>&1) &"
        if self.verbose:
            print("{} {} {}".format(preprompt,lcuprompt,cmdline))

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

    def _exec_lcu_ssh(self, cmdline, backgroundJOB=False, quotes="'"):
        """Execute a command on the LCU, either as a background job or in the
        foreground (blocking). Typically access is remote via ssh.
        (To speed things up use the ssh CommandMaster option.)
        """
        LCUprompt = "On LCU> "
        shellinvoc = "ssh "+self.lcuURL
        if backgroundJOB is True:
            cmdline = "(( "+cmdline+" ) > "+self.lcuHome+"lofarctl.log 2>&1) &"
        if self.DryRun:
            prePrompt = "(dryrun) "
        else:
            prePrompt = ""
        if self.verbose:
            print(prePrompt+LCUprompt+cmdline)
        if self.DryRun is False and self.accessible:
            if backgroundJOB == 'locally':
                # Runs in background locally rather than in background on LCU
                lcuproc = subprocess.call(shellinvoc+" "+cmdline+" &",
                                          shell=True)
            else:
                if quotes == "'":
                    lcuproc = subprocess.call(shellinvoc+" "+"'"+cmdline+"'",
                                              shell=True)
                elif quotes == '"':
                    lcuproc = subprocess.call(shellinvoc+" "+'"'+cmdline+'"',
                                              shell=True)
                else:
                    lcuproc = subprocess.call(shellinvoc+" "+cmdline,
                                              shell=True)
        elif not self.accessible:
            print("Warning: not running as "+self.lcuURL
                  + " since it is not accesible.")

    def _stdoutLCU(self, cmdline):
        """Execute a command on the LCU and check its output."""
        LCUprompt = "On LCU> "
        shellinvoc = "ssh "+self.lcuURL
        if self.DryRun:
            prePrompt = "(dryrun) "
        else:
            prePrompt = ""
        if self.verbose:
            print(prePrompt+LCUprompt+cmdline)
        if self.DryRun is False:
            try:
                output = subprocess.check_output(shellinvoc+" '"+cmdline+"'",
                                                 shell=True).rstrip()
                output = str(output.decode('UTF8'))
            except subprocess.CalledProcessError as e:
                raise e
        else:
            output = "None"
        return output

    def outfromLCU(self, cmdline, integration, duration):
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
                got = cmd.stderr.readline()
            except IOError:
                raise IOError()
            else:
                # print got
                if "shape(stats)=" in got:
                    if count % 4 == 0:
                        print(str(int(round(duration-count/4.0*integration, 0)
                                      )) + "sec left out of " + str(duration))
                    count += 1

    def _list_dat_files(self, dumpdir):
        return self._stdoutLCU("ls "+dumpdir).split('\n')

    def getdatalist(self):
        """Get list of data recorded on LCU. Returns a separate list for data dumps
        and ACC data. Note that because data dumps are named by date and time
        (uses YYYYmmdd_HHMMSS) and because default 'ls' sorts names in alphabetical
        order, the list of files will be ordered such that oldest data is first.
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
        """Get mode of station switch. Can be 'ILT' or 'local' mode."""
        if self.DryRun:
            return 'local'
        try:
            getstationmode_out = self._stdoutLCU("getstationmode")
        except:
            # Can happen if no ServiceBroker process running
            print("Warning: ServiceBroker not running")
            getstationmode_out = ""
        stationmode = getstationmode_out.split()[-1]
        return stationmode

    def rm(self, source):
        """Remove specified source file(s) on LCU."""
        self.exec_lcu("rm -fr " + source)

    def cleanup(self):
        """Clean up all local data dumps. It is usually the ObsSession()
        objects responsibility to call this."""
        self.exec_lcu("rm -fr " + self.lcuDumpDir + "/*")
        self.exec_lcu("rm " + self.ACCsrcDir + "/*.dat")

    def get_RSPDriver_conf(self):
        """Get RSPDriver configuration settings."""
        filecontents = self._stdoutLCU("cat "+self.RSPDriver_conf)
        rspdriver = parse_lofar_conf_files(filecontents)
        return rspdriver

    def getDISABLEDRCUs(self, rcumode):
        """Return list of RCUs to be disabled (as determined by ASTRON)."""
        filename = (self.lofarstationtestdir + "DISABLED/disabled-mode"
                    + str(rcumode) + ".txt")
        filecontents = self._stdoutLCU("cat "+filename)
        disabledrcus = filecontents.split(',')
        return disabledrcus

    def selectrcustr(self, rcumode):
        disabledrcustr = self.getDISABLEDRCUs(rcumode)
        if disabledrcustr[0] != '':
            disabledrcus = [int(rcustr) for rcustr in disabledrcustr]
            allrcus = range(nrofrcus)
            enabledrcus = [rcu for rcu in allrcus if rcu not in disabledrcus]
            # difenabrcu = numpy.diff(enabledrcus)
            enabledrcuflagstr = str(enabledrcus[0])+":"
            for rcuidx in range(1, len(enabledrcus)-1):
                if enabledrcus[rcuidx]-enabledrcus[rcuidx-1] != 1:
                    enabledrcuflagstr += \
                      str(enabledrcus[rcuidx-1])+","+str(enabledrcus[rcuidx])\
                      + ":"
            enabledrcuflagstr += str(enabledrcus[-1])
        else:
            enabledrcuflagstr = "0:{}".format(nrofrcus-1)
        return enabledrcuflagstr

    def get_mac_version(self):
        """Get MAC version of station."""
        macversionstr = self._stdoutLCU("swlevel -V")
        if self.DryRun:
            macversionstr = "Mock-version-2.0.0"
        return macversionstr

    def get_swlevel(self):
        """Get current Software Level of station. Returns a string which
        under normal local mode operations is an integer between 0-3."""
        swlevel = self._stdoutLCU("swlevel -S")
        if self.DryRun:
            swlevel = 'Mock 3'
        return swlevel

    def set_swlevel(self, swleveltarget=3, FullReboot=False):
        """Set station's software level. swleveltarget=3 is the swlevel for which
        most observations take place."""
        swlevel_changed = False
        if not FullReboot:
            print("Checking swlevel (prior to running observations)")
            if not self.DryRun:
                swlevel = self.get_swlevel()
                if swlevel != str(swleveltarget):
                    self.exec_lcu("swlevel " + str(swleveltarget))
                    swlevel_changed = True
        else:
            # For completeness swlevel 0, but swlevel 1 is faster
            self.exec_lcu("swlevel 0; swlevel " + str(swleveltarget))
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
        argsdict = {'select': select, 'mode': mode, 'tbbmode': tbbmode}
        # Form commandline argument string ignoring blank arguments:
        argsstr = ''
        for arg in argsdict.keys():
            if argsdict[arg] is not None:
                argsstr += " --{}={}".format(arg, argsdict[arg])
        if argsstr != '':
            rspctl_cmd = "rspctl"+argsstr
            self.exec_lcu(rspctl_cmd)
        else:
            rspctl_cmd = None
        return rspctl_cmd

    def rcusetup(self, bits, attenuation):
        """Setup basic RCU setting: bits is 8 or 16, and attenuation is 0 to 31
        (0 means no attenuation & increasing number means more attenutation)"""
        rcu_setup_cmds = ""
        rcu_setup_cmds += "rspctl --bitmode="+str(bits)+" ; "
        # NOTE Looks like bitmode and rcuattenuation have to be set in separate
        #      commands.
        if attenuation:
            # NOTE attenuation only set when beamctl is runnning.
            rcu_setup_cmds += "rspctl --rcuattenuation="+str(attenuation)+" ; "
        #rcu_setup_CMD = self.rspctl_cmd(str(bits), attenuation)
        self.exec_lcu(rcu_setup_cmds)
        if self.DryRun:
            self.bits = bits
        waittime = 1
        print("Waiting {}s for rspctl to settle...".format(waittime))
        time.sleep(waittime)  # Wait for rspctl to settle
        return rcu_setup_cmds

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
                rspbits.append( (rsp_nr, int(bitsstr.lstrip())) )
            # Will assume all rsp board have same bit depth:
            bits = rspbits[0][1]
        self.bits = bits
        return self.bits

    def _setup_beamctl(self, beamlets, subbands, band, anadigdir, rcus,
                       beamdurstr=''):
        """Create a beamctl command string from the given arguments."""
        if beamdurstr != '':
            beamdurstr = ',' + beamdurstr
        anadir = anadigdir
        digdir = anadigdir

        try:
            # See if band is actually old rcumode 3,5,7 etc
            band = rcumode2band(band)
        except ValueError:
            pass    # It's not an rcumode. Assume it's a proper band descriptor
        antset = band2antset(band)
        beamctl_cmd = ("beamctl --antennaset=" + antset +" --rcus=" + rcus
                       + " --band=" + band +" --beamlets=" + beamlets
                       + " --subbands=" + subbands
                       + " --anadir=" + anadir + beamdurstr
                       + " --digdir=" + digdir + beamdurstr)
        return beamctl_cmd

    def run_beamctl(self, beamlets, subbands, band, anadigdir, rcus='0:191',
                    beamdurstr='', backgroundJOB=True):
        """Start a beam using beamctl command. Blocks until ready."""
        beamctl_cmd = self._setup_beamctl(beamlets, subbands, band, anadigdir, rcus,
                                          beamdurstr)
        self.exec_lcu(beamctl_cmd, backgroundJOB)
        waittime = 11
        print("Waiting {}s for beam to settle...".format(waittime))
        time.sleep(waittime)  # Wait for beam to settle
        return beamctl_cmd

    def rec_bst(self, integration, duration, directory=None):
        """Convenience function to record BST data on LCU."""
        if directory is None:
            directory = self.lcuDumpDir
        rspctl_cmd = ("rspctl --statistics=beamlet"
                      + " --integration="+str(integration)
                      + " --duration="+str(duration)
                      + " --directory="+directory)
        self.outfromLCU(rspctl_cmd, integration, duration)
        if self.DryRun:
            self.mockstatistics('bst', integration, duration, directory)
        return rspctl_cmd

    def rec_sst(self, integration, duration, directory=None):
        """Convenience function to record SST data on LCU."""
        if directory is None:
            directory = self.lcuDumpDir
        rspctl_cmd = ("rspctl --statistics=subband"
                      + " --integration="+str(integration)
                      + " --duration="+str(duration)
                      + " --directory="+directory)
        self.exec_lcu(rspctl_cmd)
        if self.DryRun:
            self.mockstatistics('sst', integration, duration, directory)
        return rspctl_cmd

    def rec_xst(self, sb, integration, duration, directory=None):
        """Convenience function to record BST data on LCU."""
        if directory is None:
            directory = self.lcuDumpDir
        rspctl_cmds = ""
        # NOTE:  Seems like this has to be sent before xstats
        rspctl_cmd = ("rspctl --xcsubband="+str(sb))
        self.exec_lcu(rspctl_cmd)
        rspctl_cmds += rspctl_cmd + "\n"
        rspctl_cmd = ("rspctl --xcstatistics"
                      + " --integration="+str(integration)
                      + " --duration="+str(duration)
                      + " --directory="+directory)
        self.exec_lcu(rspctl_cmd)
        if self.DryRun:
            self.mockstatistics('xst', integration, duration, directory)
        rspctl_cmds += rspctl_cmd
        return rspctl_cmds

    def mockstatistics(self, statistics, integration, duration, directory=None):
        """Make mock statistics data file(s)."""
        if directory is None:
            directory = self.lcuDumpDir
        self.DryRun = False
        nrtimsamps = int(duration/integration)
        import ilisa.observations.modeparms as modeparms
        import datetime
        dd_cmdbase = 'dd if=/dev/zero'
        bits = self.get_bits()
        dtstamp = datetime.datetime.utcnow().strftime('%Y%m%d_%H%M%S')
        if statistics == 'bst':
            # Write mock bst files. (mock files contain just zeros)
            nrbls = modeparms.FrequencyBand().get_maxbeamletsbybits(bits)
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
        self.DryRun = True

    def run_tbbctl(self, select=None, alloc=False, free=False, record=False, stop=False,
                   mode=None, storage=None, readall=None, cepdelay=None,
                   backgroundJOB = False):
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
        r"sed -i.orig 's/^CalServer.DisableACMProxy=1/CalServer.DisableACMProxy=0/ ; s/^CalServer.WriteACCToFile=0/CalServer.WriteACCToFile=1/ ; s,^CalServer.DataDirectory=.*,CalServer.DataDirectory={}, ' {}"\
            .format(self.ACCsrcDir, self.CalServer_conf), quotes='"')
        else:
            self.exec_lcu(
        r"sed -i 's/^CalServer.DisableACMProxy=0/CalServer.DisableACMProxy=1/; s/^CalServer.WriteACCToFile=1/CalServer.WriteACCToFile=0/; s,^CalServer.DataDirectory=.*,CalServer.DataDirectory=/localhome/data,' {}"\
            .format(self.CalServer_conf), quotes='"')

    def getCalTableInfo(self, rcumode):
        """Fetch and return the caltable info from the LCU."""
        if int(rcumode) == 4:
            # Band 30_90 not correctly implemented in "beamctl --calinfo".
            # It uses the 10_90 caltab anyways so:
            rcumode = 3
        if self.usescriptonlcu:
            calinfo = subprocess.check_output("ssh "
                                                   + self.lcuURL
                                                   + " "+"infoCalTable.sh"+" "
                                                   + str(rcumode),
                                                   shell=True)
        else:
            calinfoout =self._stdoutLCU("beamctl --calinfo")
            if self.DryRun:
                return ""
            # Convert output into a list of dict per antset
            calinfolist = []
            # Strip off first initial lines and split on blank lines
            calinfooutlist = (''.join(calinfoout.splitlines(1)[2:])).split('\n\n')
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

    def selectCalTable(self, which):
        """This is specific to an lcu which has the script SelectCalTable.sh
        with which a user can switch between different caltables. se607c has
        this at /opt/lofar_local/bin/
        """
        if self.usescriptonlcu:
            if which == 'default':
                SelCalTabArg = "0"
            elif which == 'local':
                SelCalTabArg = "1"
            self.exec_lcu("SelectCalTable.sh" + " " + SelCalTabArg)
        else:
            # TODO Implement select CalTable without using script on lcu
            pass

    def turnoffLBA_LNAs(self,):
        """Turn-off the LNAs on LBA. (Used as an indication of system
        temperature."""
        # TODO allow selection of rcus rather than always all.
        rspctl_CMD = "rspctl --rcu=0x00034880 --sel=0:191"
        time.sleep(30)
        self.exec_lcu(rspctl_CMD)
        print("Warning: Turning OFF LBA LNAs.")
        time.sleep(30)
        return rspctl_CMD

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

    def turnoffElinTile_byEl(self, elemsOn):
        """"Turn off all elements per tile except the one specificied in list.
        Execution is done by element, which is less intuitive but faster."""
        self.run_rspctl(mode='5')
        for elNr in range(self.elementsInTile):
            tiles = [ind for ind in range(self.nrTiles) if elNr == elemsOn[ind]]
            if len(tiles) == 0:
                continue
            tileMap = [self.setElem_OFF for elemNr in range(self.elementsInTile)]
            tileMap[elNr] = self.setElem_ON
            rcus = self._tiles2rcus(tiles)
            lcucmd = "rspctl --hbadelay="\
                     + str(tileMap).strip('[]').replace(" ", "")\
                     + " --select="+str(rcus).strip('[]').replace(" ", "")
            self.exec_lcu(lcucmd)
