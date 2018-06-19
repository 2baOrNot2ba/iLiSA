"""LOFAR module for basic hi-level functionality for observing in LOCAL station
mode.
"""

import time
import subprocess
import argparse
import numpy

# LOFAR convience
ALLRCUs = "0:191"

# LOFAR constants
Nqfreq = 100.0e6  # Nyquist frequency in Hz
TotNrOfsb = 512  # Total number of subbands. (Subbands numbered 0:511)
nrofrcus = 192  # Number of RCUs
NrBeamletsPerLane = 61
baseSB = 0  # Nominally 0.
beamlets_lane0 = range(0,                   1*NrBeamletsPerLane)
beamlets_lane1 = range(1*NrBeamletsPerLane, 2*NrBeamletsPerLane)
beamlets_lane2 = range(2*NrBeamletsPerLane, 3*NrBeamletsPerLane)
beamlets_lane3 = range(3*NrBeamletsPerLane, 4*NrBeamletsPerLane)
maxNrSubbands = 4*NrBeamletsPerLane
elementsInTile = 16
nrTiles = 96
beamboottime = 12.0  # Time it takes for beams to settle, depends on LBA or HBA


def band2antset(band):
    """Map band to antennaset, which is used in beamctl arguments.
    Assumption is that one wants to use as many of antennas in field as
    possible.
    """
    if band == "10_90" or band == "30_90":
        antset = "LBA_INNER"
    elif band == "110_190" or band == "170_230" or band == "210_250":
        antset = "HBA_JOINED"
    else:
        raise ValueError("Undefined band: {}.".format(band))
    return antset


def band2rcumode(band):
    """Map band to rcumode string (Inverse of rcumode2band())."""
    if band == "10_90":
        rcumode = "3"
    elif band == "30_90":
        rcumode = "4"
    elif band == "110_190":
        rcumode = "5"
    elif band == "170_230":
        rcumode = "6"
    elif band == "210_250":
        rcumode = "7"
    else:
        raise ValueError('Undefined band %{}'.format(band))
    return rcumode


def rcumode2band(rcumode):
    """Map rcumode to band string as used in beamctl arguments."""
    rcumode = str(rcumode)
    if rcumode == "3":
        band = "10_90"
    elif rcumode == "4":
        band = "30_90"
    elif rcumode == "5":
        band = "110_190"
    elif rcumode == "6":
        band = "170_230"
    elif rcumode == "7":
        band = "210_250"
    else:
        raise ValueError('Undefined rcumode %{}'.format(rcumode))
    return band


def rcumode2freqrange(rcumode):
    return rcumode2band(rcumode).split('_')


def band2freqrange(band):
    return tuple(1e6*f for f in map(int, band.split('_')))


def rcumode2antset(rcumode):
    """Map rcumode to antennaset, which is used in beamctl arguments.
    Assumption is that one wants to use as many of antennas in field as
    possible. (This function may soon be deprecated.)
    """
    # NOTE new/more antennasets are now available.
    if rcumode == "3":
        antset = "LBA_INNER"
    elif rcumode == "5" or rcumode == "7":
        antset = "HBA_JOINED"
    else:
        raise ValueError("Undefined rcumode: {}.".format(rcumode))
    return antset


def rcumode2sbfreqs(rcumode):
    """Get the frequencies (in Hz) of the subbands for the given rcumode.
    Returns an array of frequencies where index is subband number."""
    NZ = (int(rcumode)-3)/2
    # Note the endpoint=False here. Before it 2018-03-22 it was missing.
    freqs = numpy.linspace(NZ*Nqfreq, (NZ+1)*Nqfreq, TotNrOfsb, endpoint=False)
    return freqs


def rcumode2NyquistZone(rcumode):
    NZ = int((int(rcumode)-3)/2)
    return NZ


def NyquistZone2rcumode(NZ):
    rcumode = NZ*2+3
    return str(rcumode)


def freq2sb(freq):
    """Convert frequency in Hz to subband number and Nyquist zone."""
    absSB = int(round(freq/Nqfreq*TotNrOfsb))
    sb = absSB % TotNrOfsb
    NqZone = absSB / TotNrOfsb
    return sb, NqZone


def sb2freq(sb, NqZone):
    """Convert subband in a given Nyquist zone to a frequency."""
    freq = Nqfreq*(int(sb)/float(TotNrOfsb)+int(NqZone))
    return freq


def freq2beamParam(freq):
    """Get beam parameters antset, rcumode and subband for with frequency (Hz).
    """
    sb, sbNZ = freq2sb(freq)
    rcumode = NyquistZone2rcumode(sbNZ)
    antset = rcumode2antset(rcumode)
    return (antset, rcumode, str(sb))


def freqBand2sb(freqBand, wordsize_bits=16):
    """Convert a frequency band to subbands."""
    maxsbs = maxNrOfBeamlets(wordsize_bits)
    freqLo = freqBand[0]
    sbLo, sbLoNZ = freq2sb(freqLo)
    if len(freqBand) == 2:
        freqHi = freqBand[1]
        sbHi, sbHiNZ = freq2sb(freqHi)
    else:
        sbHi = sbLo+min(maxsbs, TotNrOfsb-sbLo)-1
        sbHiNZ = sbLoNZ
    NrSBs = sbHi-sbLo+1
    if NrSBs > maxsbs:
        print "Frequency range not permitted: too many subbands"
        raise ValueError
    if sbLoNZ != sbHiNZ:
        print "Frequency range not permitted: different Nyquist zones"
        raise ValueError
    rcumode = NyquistZone2rcumode(sbLoNZ)
    antset = rcumode2antset(rcumode)
    beamlets = "0:"+str(NrSBs-1)
    subbands = str(sbLo)+":"+str(sbHi)
    return (antset, rcumode, beamlets, subbands)


def maxNrOfBeamlets(wordsize_bits):
    """Return maximum number of subbands, depending on word size of ADC
    samples: 16-bit, 8-bit modes."""
    maxNrOfSBs_16 = 244
    if wordsize_bits == 16:
        maxNrOfSBs = maxNrOfSBs_16
    elif wordsize_bits == 8:
        maxNrOfSBs = 2*maxNrOfSBs_16
    else:
        print "Unknown wordsize: ", str(wordsize_bits)
        raise('ValueError')
    return maxNrOfSBs


def tiles2rcus(tiles):
    rcus = []
    for tile in tiles:
        rcus.extend([2*tile, 2*tile+1])  # Set same delay for both X&Y pol rcu
    return rcus


def parse_multibeamctl_args(beamctl_strs):
    antennaset_list = []
    rcus_list = []
    rcumode_list = []
    beamlets_list = []
    subbands_list = []
    anadir_list = []
    digdir_list = []
    beamctl_str_lst = beamctl_strs.split("; ")
    for beamctl_str in beamctl_str_lst:
        (antennaset, rcus, rcumode, beamlets, subbands, anadir, digdir
         ) = parse_beamctl_args(beamctl_str)
        antennaset_list.append(antennaset)
        rcus_list.append(rcus)
        rcumode_list.append(rcumode)
        beamlets_list.append(beamlets)
        subbands_list.append(subbands)
        anadir_list.append(anadir)
        digdir_list.append(digdir)
    return (antennaset_list, rcus_list, rcumode_list, beamlets_list,
            subbands_list, anadir_list, digdir_list)


def parse_beamctl_args(beamctl_str):
    """Parse beamctl command arguments"""
    beamctl_str_normalized = beamctl_str.replace('=', ' ')
    beamctl_parser = argparse.ArgumentParser()
    beamctl_parser.add_argument('--antennaset')
    beamctl_parser.add_argument('--rcus')
    beamctl_parser.add_argument('--rcumode')  # Obsolete
    beamctl_parser.add_argument('--band')
    beamctl_parser.add_argument('--beamlets')
    beamctl_parser.add_argument('--subbands')
    beamctl_parser.add_argument('--integration')
    beamctl_parser.add_argument('--duration')
    beamctl_parser.add_argument('--anadir')
    beamctl_parser.add_argument('--digdir')
    args = beamctl_parser.parse_args(beamctl_str_normalized.split()[1:])
    rcumode = band2rcumode(args.band)
    return (args.antennaset, args.rcus, rcumode, args.beamlets,
            args.subbands, args.anadir, args.digdir)


#######################################
# Basic station control START


class Station(object):
    """This class manages an International LOFAR station."""
    lofarroot = "/opt/lofar_local/"
    lofarstationtestdir = "/localhome/stationtest/"
    # ACCsrcDir is set in CalServer.conf and used when CalServer is running.
    ACCsrcDir = "/localhome/data/ACCdata/"
    CalServer_conf = lofarroot + "/etc/CalServer.conf"

    def setupaccess(self, accessconf):
        """Initialize with user-station configuration."""
        self.stnid = accessconf['stnid']
        self.user = accessconf['user']
        self.lcuURL = self.user+"@"+accessconf['hostname']
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
            print "Established access to LCU."


    def checkaccess(self):
        """Check that this object has access to LCU.
        Note: It does this by trying to get the MAC version number."""
        self.MACversion = ""
        try:
            self.MACversion = self.getMACversion()
        except Exception:
            print("Error: Station object cannot access station with URL: "
                  + self.lcuURL)
        if self.MACversion is not "":
            self.accessible = True
        else:
            self.accessible = False
        return self.accessible

    def __init__(self, lcuaccessconf=None):
        if lcuaccessconf is not None:
            self.setupaccess(lcuaccessconf)

    def __del__(self):
        pass

    def execOnLCU(self, cmdline, backgroundJOB=False, quotes="'"):
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
            print prePrompt+LCUprompt+cmdline
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
            print prePrompt+LCUprompt+cmdline
        if self.DryRun is False:
            try:
                output = subprocess.check_output(shellinvoc+" '"+cmdline+"'",
                                                 shell=True).rstrip()
            except subprocess.CalledProcessError, e:
                raise Exception('Access LCU error: {}'.format(e))
        else:
            output = "None"
        return output

    def outfromLCU(self, cmdline, integration, duration):
        print "LCUo>", cmdline
        cmd = subprocess.Popen("ssh "+self.lcuURL+" "+"'"+cmdline+"'",
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, shell=True)
        count = 0
        outstrname = 'stderr'
        while cmd.poll() is None:
            if outstrname == 'stdout':
                outstr = cmd.stdout
            elif outstrname == 'stderr':
                outstr = cmd.stderr
            else:
                print "Unknown output name"
                exit(1)
            try:
                got = cmd.stderr.readline()
            except IOError:
                print("ioerror")
                exit(1)
            else:
                # print got
                if "shape(stats)=" in got:
                    if count % 4 == 0:
                        print(str(int(round(duration-count/4.0*integration, 0)
                                      )) + "sec left out of " + str(duration))
                    count += 1

    def getdatalist(self):
        ls_lcuDumpDir = self._stdoutLCU("ls "+self.lcuDumpDir).split('\n')
        ls_ACCsrcDir = self._stdoutLCU("ls "+self.ACCsrcDir).split('\n')
        return ls_lcuDumpDir, ls_ACCsrcDir

    def whoServiceBroker(self):
        """Check who is running the Service Broker on the LCU. This is an
        indication of who is currently using the station."""
        try:
            ps_out = \
              self._stdoutLCU("/bin/ps -CServiceBroker --no-headers -ouser")
        except:
            # Can happen if no ServiceBroker process running. Pretend ps
            # returns blank
            ps_out = ""
        ps_out_lns = ps_out.splitlines()
        if len(ps_out_lns) == 0:
            sb_user = 'None'
        else:
            sb_user = ps_out_lns[0]
        return sb_user

    def getstationswitchmode(self):
        """Get mode of station switch. Can be 'ILT' or 'local' mode."""
        try:
            getstationmode_out = self._stdoutLCU("getstationmode")
        except:
            # Can happen if no ServiceBroker process running
            print "Caught"
            getstationmode_out = ""
        stationmode = getstationmode_out.split()[-1]
        return stationmode

    def rm(self, source):
        """Remove specified source file(s) on LCU."""
        self.execOnLCU("rm -fr "+source)

    def cleanup(self):
        """Clean up all local data dumps. It is usually the ObsSession()
        objects responsibility to call this."""
        self.execOnLCU("rm -fr "+self.lcuDumpDir+"/*")
        self.execOnLCU("rm "+self.ACCsrcDir+"/*.dat")

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
            enabledrcuflagstr = ALLRCUs
        return enabledrcuflagstr

    def getMACversion(self):
        """Get MAC version of station."""
        macversionstr = self._stdoutLCU("swlevel -V")
        if self.DryRun:
            macversionstr = "Mock-version-2.0.0"
        return macversionstr

    def getswlevel(self):
        """Get current Software Level of station. Returns a string which
        under normal local mode operations is an integer between 0-3."""
        swlevel = self._stdoutLCU("swlevel -S")
        if self.DryRun:
            swlevel = 'Mock 3'
        return swlevel

    def bootToObservationState(self, swleveltarget=3, FullReboot=False):
        """Get station to observation state."""
        # Stop any beam and remove any data in lcu datadump
        self.execOnLCU("killall beamctl")

        if FullReboot is not True:
            print "Checking swlevel (prior to running observations)"
            if not self.DryRun:
                swlevel = self.getswlevel()
            else:
                swlevel = "undefined"
            print "Found swlevel="+swlevel
            if swlevel != str(swleveltarget):
                # FullReboot = True
                self.execOnLCU("swlevel "+str(swleveltarget))
        if FullReboot is True:
            # May need to be swlevel 0, but swlevel 1 is faster
            self.execOnLCU("swlevel 0; swlevel "+str(swleveltarget))
        # TODO check if we own the swlevel

    def stopBeam(self):
        """Stop any running beamctl processes."""
        # Stop any beamctl and clean up datadump dir on lcu.
        self.execOnLCU("killall beamctl")
        # Put caltables back to default
        self.selectCalTable('default')
        # print("Beam off at %s"%time.asctime(time.localtime(time.time())))

    def shutdownObservationState(self):
        # Go to swlevel 0. No LOFAR services will be running after this.
        self.execOnLCU('swlevel 0')

##Basic station data taking commands BEGIN

    def rspctl_cmd(self, bits=16, attenuation=None):
        """Return rspctl command to setup RCUs: bits is 8 or 16, and
        attenuation is 0 to 31.
        (Nothing is executed on LCU)"""
        # TODO should maybe check if bit mode is already correct...
        # NOTE Looks like bitmode and rcuattenuation have to be set in separate
        #      commands.
        rspctl_CMDs = ""
        rspctl_CMDs += "rspctl --bitmode="+str(bits)+" ; "
        if attenuation:
            # NOTE attenuation only set when beamctl is runnning.
            rspctl_CMDs += "rspctl --rcuattenuation="+str(attenuation)+" ; "
        return rspctl_CMDs

    def rcusetup(self, bits, attenuation):
        """Setup basic RCU setting: bits is 8 or 16, and attenuation is 0 to 31
        (0 means no attenuation & increasing number means more attenutation)"""
        rcu_setup_CMD = self.rspctl_cmd(str(bits), attenuation)
        self.execOnLCU(rcu_setup_CMD)
        return rcu_setup_CMD

    def setupbeamlets(self, beamletIDs, subbandNrs, band, pointing,
                      RCUSflag=ALLRCUs, beamletDurStr=""):
        """ """
        if beamletDurStr != "":
            beamletDurStr = ","+beamletDurStr
        anadir = pointing
        digdir = pointing

        try:
            # See if band is actually old rcumode 3,5,7 etc
            band = rcumode2band(band)
        except ValueError:
            pass    # It's not an rcumode. Assume it's a proper band descriptor
        antset = band2antset(band)
        beamctl_CMD = ("beamctl --antennaset="+antset+" --rcus="+RCUSflag
                       + " --band="+band+" --beamlets="+beamletIDs
                       + " --subbands="+subbandNrs
                       + " --anadir="+anadir+beamletDurStr
                       + " --digdir="+digdir+beamletDurStr)
        return beamctl_CMD

    def runbeamctl(self, beamletIDs, subbandNrs, rcumode, pointing,
                   RCUSflag=ALLRCUs, beamletDurStr="", backgroundJOB=True):
        """Start a beam using beamctl command. Blocks until ready."""
        beamctl_CMD = self.setupbeamlets(beamletIDs, subbandNrs, rcumode,
                                         pointing, RCUSflag, beamletDurStr)
        self.execOnLCU(beamctl_CMD, backgroundJOB)
        waittime = 10
        print "Waiting "+str(waittime)+" seconds"
        time.sleep(waittime)  # Wait for beam to settle
        return beamctl_CMD

    def streambeam(self, freqBand, pointings, recDuration=float('inf'),
                   bits=16, attenuation=0, DUMMYWARMUP=False):
        # Multiple sbs with same pointing
        # or one sb with pointing
        beamctl_CMDs = ""
        MultiBeamctl = False

        if type(freqBand) is tuple:
            antset, rcumode, beamletIDs, subbands = freqBand2sb(freqBand, bits)
        elif type(freqBand) is float:
            antset, rcumode, subbands = freq2beamParam(freqBand)
            if type(pointings) is tuple:
                MultiBeamctl = True
            else:
                # Nominally:
                # (beamctl does not allow allocation of only one beamlet,
                # so we add one more above)
                beamletIDs = '0:1' 
                subbands += ':'+ str(int(subbands)+1) 
                # (special test used):
                # nrBLs = 61*4*16/bits
                # beamletIDs=','.join([str(b) for b in range(nrBLs)])
                # #subbands = ((subbands+',')*nrBLs).rstrip(',')

        # Select good rcus
        enabledrcus = ALLRCUs
        # enabledrcus = self.selectrcustr(rcumode)

        if MultiBeamctl:
            for beamletNr in range(0, len(pointings)):
                beamctl_CMDs = (beamctl_CMDs
                                + self.setupbeamlets(str(beamletNr), subbands,
                                                     rcumode,
                                                     pointings[beamletNr])
                                + " & \\"+"\n"
                                )
            # beamctl_CMDs=beamctl_CMDs[:-2] #Trim off trailing "; "
        else:
            beamctl_main = self.setupbeamlets(beamletIDs, subbands, rcumode,
                                              pointings, RCUSflag=enabledrcus)
            backgroundJOB = True
            beamctl_CMDs = beamctl_main
            # beamctl_CMDs="(("+beamctl_CMDs+") & )"
        # print "(About to stream beams)"
        # recStop_CMD=":" #bash no-op
        if backgroundJOB is True:
            beamctl_CMDs += " & "
        # Setup rspctl settings
        rspctl_SET = self.rspctl_cmd(bits, attenuation)
        # beamctl_CMDs = rspctl_SET + beamctl_CMDs # rspctl cmd before beamctl
        beamctl_CMDs = beamctl_CMDs + rspctl_SET  # rspctl cmd *after* beamctl
        if recDuration != float('inf'):
            recDuration_CMD = " sleep "+str(recDuration+beamboottime)+"; "
            recStop_CMD = "killall beamctl"
            # Start beam streaming and stop after recDuration seconds
            # beamctl_CMD_BG="(("+beamctl_CMDs+") &)" #Put beamctl in bg
            beamctl_CMDs += recDuration_CMD
            beamctl_CMDs += recStop_CMD
            backgroundJOB = True  # FIX this
        if DUMMYWARMUP:
            print "Warm-up with dummy beam"
            self.execOnLCU(beamctl_CMDs+" sleep "+str(1)
                           + "; killall beamctl", False)
        self.execOnLCU(beamctl_CMDs, backgroundJOB)
        if backgroundJOB:
            beamSettleTime = 6.0    # Time in seconds before beam settles.
                                    # Can be only 5s if bitmode unchanged
            print "Waiting "+str(beamSettleTime)+"s for beam to settle"
            time.sleep(beamSettleTime)
        # FIX separation between beamctl_CMDs & _main
        return beamctl_CMDs, rspctl_SET, beamctl_main

### Basic station "statistic" datataking BEGIN
    def rec_bst(self, integration, duration):
        rspctl_CMD = ("rspctl --statistics=beamlet"
                      + " --integration="+str(integration)
                      + " --duration="+str(duration)
                      + " --directory="+self.lcuDumpDir)
        self.outfromLCU(rspctl_CMD, integration, duration)
        return rspctl_CMD

    def rec_sst(self, integration, duration):
        rspctl_CMD = ("rspctl --statistics=subband"
                      + " --integration="+str(integration)
                      + " --duration="+str(duration)
                      + " --directory="+self.lcuDumpDir)
        self.execOnLCU(rspctl_CMD)
        return rspctl_CMD

    def rec_xst(self, sb, integration, duration):
        rspctl_CMDs = ""
        # NOTE:  Seems like this has to be sent before xstats
        rspctl_CMD = ("rspctl --xcsubband="+sb)
        self.execOnLCU(rspctl_CMD)
        rspctl_CMDs += rspctl_CMD + "\n"
        rspctl_CMD = ("rspctl --xcstatistics"
                      + " --integration="+str(integration)
                      + " --duration="+str(duration)
                      + " --directory="+self.lcuDumpDir)
        self.execOnLCU(rspctl_CMD)
        rspctl_CMDs += rspctl_CMD
        return rspctl_CMDs
### Basic station "statistic" datataking END

### TBB control BEGIN
    def setupTBBs(self):
        """Setup transient buffer boards for recording."""
        if self.usescriptonlcu:
            self.execOnLCU("scripts/tbb_setup.sh")
        else:
            print "Freeing TBBs"
            self.execOnLCU("tbbctl --free")
            print "Setting TBB transient mode on rspctl"
            self.execOnLCU("rspctl --tbbmode=transient")
            time.sleep(1)
            print "Allocating TBBs"
            self.execOnLCU("tbbctl --alloc")
            time.sleep(1)
            print "Setting TBB transient mode on tbbctl"
            self.execOnLCU("tbbctl --mode=transient")
            time.sleep(1)
            print "Start TBB recording"
            self.execOnLCU("tbbctl --record")
            print "Finished setting up TBBs & started recording"

    def freezeTBBdata(self):
        if self.usescriptonlcu:
            self.execOnLCU("scripts/tbb_stop.sh")
        else:
            print "Stopping TBB recording"
            self.execOnLCU("tbbctl --stop")
            print "Stopping any dummy beam"
            self.stopBeam()

    def startTBBdataStream(self, duration):
        """Stream duration seconds of TBB data out of the LCU to
        datataking node."""
        # Set delay between subsequent frames. One delay unit is 5us.
        udpdelay = 500  # (Previously 100)
        nrpages = str(int(duration*2*Nqfreq/1024))  # One page is 1024 samples.
                                                    # Normal sampling frequency
                                                    # is 200MHz.
        if self.usescriptonlcu:
            self.execOnLCU("scripts/tbb_dumpall.sh"+" "+nrpages)
        else:
            print "Streaming TBB data"
            self.execOnLCU("tbbctl --storage=lofarA1 --select=0:15,16:31,32:47"
                           )
            self.execOnLCU("tbbctl --storage=lofarA2 --select=48:63,64:79,80:95"
                           )
            self.execOnLCU("tbbctl --storage=lofarA3 --select=96:111,112:127,128:143")
            self.execOnLCU("tbbctl --storage=lofarA4 --select=144:159,160:175,176:191")

            self.execOnLCU("tbbctl --cepdelay="+str(udpdelay))

            self.execOnLCU("tbbctl --readall="+nrpages+" --select=0:15",
                           backgroundJOB='locally')
            self.execOnLCU("tbbctl --readall="+nrpages+"  --select=48:63",
                           backgroundJOB='locally')
            self.execOnLCU("tbbctl --readall="+nrpages+"  --select=96:111",
                           backgroundJOB='locally')
            self.execOnLCU("tbbctl --readall="+nrpages+"  --select=144:159",
                           backgroundJOB='locally' )

            self.execOnLCU("tbbctl --readall="+nrpages+"  --select=16:31",
                           backgroundJOB='locally')
            self.execOnLCU("tbbctl --readall="+nrpages+"  --select=64:79",
                           backgroundJOB='locally')
            self.execOnLCU("tbbctl --readall="+nrpages+"  --select=112:127",
                           backgroundJOB='locally')
            self.execOnLCU("tbbctl --readall="+nrpages+"  --select=160:175",
                           backgroundJOB='locally')

            self.execOnLCU("tbbctl --readall="+nrpages+"  --select=32:47",
                           backgroundJOB='locally')
            self.execOnLCU("tbbctl --readall="+nrpages+"  --select=80:95",
                           backgroundJOB='locally')
            self.execOnLCU("tbbctl --readall="+nrpages+"  --select=128:143",
                           backgroundJOB='locally')
            # Last one is not put in background so the parent process blocks
            # until finished.
            self.execOnLCU("tbbctl --readall="+nrpages+"  --select=176:191",
                           backgroundJOB='locally')

### TBB control END
### ACC control BEGIN
    def runACC(self, rcumode, duration, pointing, sst_integration=600):
        """Perform ACC calibration observation on station.

        ACC files are autocovariance-cubes: the covariance of all array
        elements with each as a function of subband. These files are generated
        by the MAC service called CalServer. It run at swlevel 3 and is
        configured in the file lofar/etc/CalServer.conf. Note subband
        integration is always 1s, so ACC file is dumped after 512 seconds.
        """

        # Make dump data directory for concurrent SST data
        self.execOnLCU(
          'if [ ! -d "'+self.lcuDumpDir+'" ]; then mkdir '
          + self.lcuDumpDir+'; fi', quotes="'")

        # Make sure swlevel=<2
        self.bootToObservationState(2)

        # Set CalServ.conf to dump ACCs:
        if self.usescriptonlcu:
            self.execOnLCU(self.scriptsDir+'/CalServDump.sh 1')
        else:
            self.execOnLCU(
r"sed -i.orig 's/^CalServer.DisableACMProxy=1/CalServer.DisableACMProxy=0/ ; s/^CalServer.WriteACCToFile=0/CalServer.WriteACCToFile=1/ ; s,^CalServer.DataDirectory=.*,CalServer.DataDirectory={}, ' {}".format(self.ACCsrcDir, self.CalServer_conf)
            , quotes='"')

        # Boot to swlevel 3 so the calserver service starts (
        self.bootToObservationState()

        freqs = rcumode2sbfreqs(rcumode)

        midfreq = float(freqs[len(freqs)/2])
        beamctl_CMD, rspctl_SET, beamctl_main = self.streambeam(midfreq,
                                                                pointing)

        # Run for $duration seconds
        rspctl_CMD = self.rec_sst(sst_integration, duration)
        self.stopBeam()

        # Switch back to normal state i.e. turn-off ACC dumping:
        if self.usescriptonlcu:
            self.execOnLCU(self.scriptsDir+'/CalServDump.sh 0')
        else:
            self.execOnLCU(
r"sed -i 's/^CalServer.DisableACMProxy=0/CalServer.DisableACMProxy=1/; s/^CalServer.WriteACCToFile=1/CalServer.WriteACCToFile=0/; s,^CalServer.DataDirectory=.*,CalServer.DataDirectory=/localhome/data,' {}".format(self.CalServer_conf)
            , quotes='"')

        return beamctl_CMD, rspctl_CMD

### ACC control END
##Basic station data taking commands END

## Special commands START
    def getCalTableInfo(self, rcumode):
        """Fetch and return the caltable info from the LCU."""
        # TODO Consider use of "beamctl --calinfo"
        if self.usescriptonlcu:
            caltableInfo = subprocess.check_output("ssh "
                                                   + self.lcuURL
                                                   + " "+"infoCalTable.sh"+" "
                                                   + str(rcumode),
                                                   shell=True)
        else:
            caltableInfo = ""
            for line in subprocess.check_output("ssh "+self.lcuURL
                                                + " "
                                                + "cat /opt/lofar/etc/CalTable_mode"
                                                + str(rcumode)+".dat",
                                                shell=True).split('\n'):
                if line == "HeaderStop": break
                if line == "HeaderStart": continue
                caltableInfo += line + '\n'
        return caltableInfo

    def selectCalTable(self, which):
        """This is specific to an lcu which has the script SelectCalTable.sh
        with which a user can switch between different caltables. se607c has
        this at /opt/lofar_local/bin/
        """
        if which == 'default':
            SelCalTabArg = "0"
        elif which == 'local':
            SelCalTabArg = "1"
        self.execOnLCU("SelectCalTable.sh"+" "+SelCalTabArg)

    def turnoffLBA_LNAs(self,):
        """Turn-off the LNAs on LBA. (Used as an indication of system
        temperature."""
        # TODO allow selection of rcus rather than always all.
        rspctl_CMD = "rspctl --rcu=0x00034880 --sel=0:191"
        time.sleep(30)
        self.execOnLCU(rspctl_CMD)
        print "Warning: Turning OFF LBA LNAs."
        time.sleep(30)
        return rspctl_CMD

    def setrcumode(self, rcumode):
        """Set the rcumode."""
        self.execOnLCU("rspctl --mode={}".format(rcumode))
        time.sleep(2.0)


# SEPTON
    setElem_ON = 128
    setElem_OFF = 2

    def turnoffElinTile_byTile(self, elemsOn):
        """"Turn off all elements per tile except the one specificied in list.
        Execution is done by tile, which is more intuitive but slower."""
        self.setrcumode(5)
        for tileNr in range(nrTiles):
            # Start with all elements in tile off
            # (2 is OFF)
            tileMap = [self.setElem_OFF for elemNr in range(elementsInTile)]
            # Turn on the appropriate element
            tileMap[elemsOn[tileNr]] = self.setElem_ON  # 128 is ON
            lcucmd = "rspctl --hbadelay="\
                     + str(tileMap).strip('[]').replace(" ", "")\
                     + " --select="+str(2*tileNr)+","+str(2*tileNr+1)
            self.execOnLCU(lcucmd)

    def turnoffElinTile_byEl(self, elemsOn):
        """"Turn off all elements per tile except the one specificied in list.
        Execution is done by element, which is less intuitive but faster."""
        self.setrcumode(5)
        for elNr in range(elementsInTile):
            tiles = [ind for ind in range(nrTiles) if elNr == elemsOn[ind]]
            if len(tiles) == 0:
                continue
            tileMap = [self.setElem_OFF for elemNr in range(elementsInTile)]
            tileMap[elNr] = self.setElem_ON
            rcus = tiles2rcus(tiles)
            lcucmd = "rspctl --hbadelay="\
                     + str(tileMap).strip('[]').replace(" ", "")\
                     + " --select="+str(rcus).strip('[]').replace(" ", "")
            self.execOnLCU(lcucmd)

## Special commands END

# Basic station control END
#######################################
