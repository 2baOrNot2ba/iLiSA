#!/usr/bin/python
"""Package that provides functions for setting up and running observations
via station controller objects. This package knows about the data archive and
should not run anything directly on LCU."""


import math
import time
import subprocess
import os
import yaml
import multiprocessing
import copy
import ilisa.observations.stationcontrol as stationcontrol
import ilisa.observations.dataIO as dataIO


# SEPTON configurations:
#        1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24
from ilisa.observations.modeparms import stdPointings, elementMap2str

elOn_step = [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15, 0, 1, 2, 3, 4, 5, 6, 7,
             8, 9,10,11,12,13,14,15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,
             0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15, 0, 1, 2, 3, 4, 5, 6, 7,
             8, 9,10,11,12,13,14,15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15]
elOn_gILT = [15, 0,15, 3, 9,15,14, 2, 0, 3, 4,14,10, 8, 5,15,12, 0, 2,11, 3,12,12, 1,
              5, 4, 4, 8, 6, 3, 0, 5, 3,11, 3, 2, 8,15,13, 8, 3, 2, 9, 1,14, 8, 8, 0,
             12,13, 0,11,15, 3,12, 3,13, 3,10, 5, 0,10, 1, 6, 4,10, 3,15, 3,14, 0,12,
              0, 7, 0,12, 7, 3,13, 0, 7, 3,15, 4,14, 4, 3, 8, 4, 9,12, 0,14, 9, 3,11]
elOn_Generic_Int_201512 = \
            [ 0, 5, 3, 1, 8, 3,12,15,10,13,11, 5,12,12, 5, 2,10, 8, 0, 3, 5, 1, 4, 0,
             11, 6, 2, 4, 9,14,15, 3, 7, 5,13,15, 5, 6, 5,12,15, 7, 1, 1,14, 9, 4, 9,
              3, 9, 3,13, 7,14, 7,14, 2, 8, 8, 0, 1, 4, 2, 2,12,15, 5, 7, 6,10,12, 3,
              3,12, 7, 4, 6, 0, 5, 9, 1,10,10,11, 5,11, 7, 9, 7, 6, 4, 4,15, 4, 1,15]
elOn_same_el = 0
elOn_same = [elOn_same_el for elemNr in range(stationcontrol.nrTiles)]
elemsOn = elOn_Generic_Int_201512  # elOn_same or elOn_step or elOn_gILT or ...


class Session(object):
    """observing.Session class is a client type class typically running on the
    compute node."""

    def cleanup(self):
        pass

    def checkobservingallowed(self):
        """Check whether observations are allowed. This occurs is someone else
        is using the station."""
        serviceuser = self.stationcontroller.whoServiceBroker()

        if serviceuser == 'None' or serviceuser == self.stationcontroller.user:
            stationswitchmode = self.stationcontroller.getstationswitchmode()
            if stationswitchmode == 'local':
                return True
            else:
                print "Warning: Station is not in stand-alone mode."
                return False
        else:
            print "Warning: Someone else ({}) is using LCU".format(serviceuser)
            print "         (You are running as {})".format(self.stationcontroller.user)
            return False

    def __init__(self, accessconffile=None,
                 goto_observingstate_when_starting=True):
        """Initialize a Session object, which has access to a station via
        a Station object configured with setting given by accessconfile.
        When goto_observingstate_when_starting is True, boot the
        station to swlevel 3, but only if observations are allowed on the
        station.
        """
        if accessconffile is None:
            accessconffile = os.path.expanduser('~/.iLiSA/access_config.yml')
        with open(accessconffile) as cfigfilep:
            accessconf = yaml.load(cfigfilep)
        self.observer = accessconf['OBSERVING']['observer']
        self.project = accessconf['OBSERVING']['project']
        lcuaccessconf = accessconf['LCU']
        dpuaccessconf = accessconf['DPU']
        self.stationcontroller = stationcontrol.Station(lcuaccessconf)

        self.LOFARdataArchive = dpuaccessconf['LOFARdataArchive']
        self.bf_data_dir =      dpuaccessconf['BeamFormDataDir']
        self.bf_port0 =     int(dpuaccessconf['BeamFormPort0'])
        self.bf_logfile =       dpuaccessconf['BeamFormLogFile']
        self.tbbraw2h5cmd =     dpuaccessconf['TBBraw2h5Cmd']
        self.tbbh5dumpdir =     dpuaccessconf['TBBh5dumpDir']

        self.DryRun = lcuaccessconf['DryRun']
        self.exit_check = True
        self.halt_observingstate_when_finished = True
        self.cleanup()
        # Check whether the station is being used by someone else:
        if self.checkobservingallowed() and goto_observingstate_when_starting:
            self.stationcontroller.bootToObservationState()
            # TODO: Figure out what should happen when goto_observingstate_when_starting==False
            # but later user wants to observe (but station might not be in observation state)

    def haltobservingstate(self):
        """Halt observing state on station."""
        if self.checkobservingallowed():
            self.stationcontroller.shutdownObservationState()
            self.cleanup()
            # Cleanup any data left on LCU.
            self.stationcontroller.cleanup()
        # Close down stationcontroller:
        del self.stationcontroller

    def __del__(self):
        """Normally, shutdown observation mode on station."""
        if self.halt_observingstate_when_finished:
            self.haltobservingstate()
        elif self.exit_check:
            if self.checkobservingallowed():
                swlevel = self.stationcontroller.getswlevel()
                if int(swlevel) != 0:
                    print "Warning: You are leaving station in swlevel {} != 0".format(swlevel)

    def movefromlcu(self, source, dest, recursive=False):
        """Move file(s) off LCU to DPU."""
        if not os.path.isdir(dest):
            os.mkdir(dest)
        movecmd = "scp"
        if recursive:
            movecmd += " -r"
        fullcmd = movecmd+" "+self.stationcontroller.lcuURL+":"+source+" "+dest
        cmdprompt = "on DPU>"
        if self.stationcontroller.DryRun:
            cmdprompt = "(dryrun) "+cmdprompt
        if self.stationcontroller.verbose:
            print "{} {}".format(cmdprompt, fullcmd)
        if not self.stationcontroller.DryRun:
            subprocess.call(fullcmd, shell=True)
            self.stationcontroller.rm(source)

    def get_data_timestamp(self, order=0):
        """Get timestamp of datafiles on LCU. order is the temporal order of the data
        file, order=0 is oldest, order=-1 is newest."""
        dd_dir, acc_dir = self.stationcontroller.getdatalist()
        # Assumes only one file in datadump dir with
        # format YYYYmmdd_HHMMSS_[bsx]st.dat
        # FIXME handle ACC datafiles too
        obsdate, obstime, obssuff = dd_dir[order].split('_', 2)
        obsdatetime_stamp = obsdate+'_'+obstime
        return obsdatetime_stamp

#######################################
# Begin: Basic obs modes

    def streambeams(self, freqbndobj, pointing, recDuration=float('inf'),
                    attenuation=0, DUMMYWARMUP=False):
        """Form beams with station."""
        bits = freqbndobj.bits
        if DUMMYWARMUP:
            print("Warning warnup not currently implemented")
        beamctl_cmds = []
        for bandbeamidx in range(len(freqbndobj.rcumodes)):
            antset = freqbndobj.antsets[bandbeamidx]
            rcumode = freqbndobj.rcumodes[bandbeamidx]
            beamletIDs = freqbndobj.beamlets[bandbeamidx]
            subbands =  freqbndobj.sb_range[bandbeamidx]
            rcusel = freqbndobj.rcusel[bandbeamidx]
            beamctl_main = self.stationcontroller.runbeamctl(beamletIDs, subbands,
                                                             rcumode, pointing, rcusel)
            beamctl_cmds.append(beamctl_main)
        rcu_setup_CMD = self.stationcontroller.rcusetup(bits, attenuation)
        return rcu_setup_CMD, beamctl_cmds


    def do_bst(self, freqbndobj, integration, duration, pointing, obsinfo):
        """Do a Beamlet STatistic (bst) recording on station. frequency in Hz.
        """

        CALTABLESRC = 'default'   # FIXME put this in args

        # Setup Calibration:
        ## (Only BST uses calibration tables)
        # Choose between 'default' or 'local'
        self.stationcontroller.selectCalTable(CALTABLESRC)
        # Start beamforming
        (rcu_setup_CMD, beamctl_cmds) = self.streambeams(freqbndobj, pointing)

        # Get some metadata about operational settings:
        if len(freqbndobj.rcumodes) == 1:
            print freqbndobj.rcumodes[0]
            caltabinfo = self.stationcontroller.getCalTableInfo(freqbndobj.rcumodes[0])
        else:
            # TODO implement storing of multiband caltab
            caltabinfo = ''

        # Record data
        # waittime = 0
        # print "Waiting extra", str(waittime) +" seconds" #Seems necessary
        # time.sleep(waittime)
        rspctl_CMD = self.stationcontroller.rec_bst(integration, duration)
        # Move data to archive
        obsdatetime_stamp = self.get_data_timestamp()
        obsinfo.setobsinfo_fromparams('bst', obsdatetime_stamp, beamctl_cmds, rspctl_CMD,
                                      caltabinfo)
        bsxSTobsEpoch, datapath = obsinfo.getobsdatapath(self.LOFARdataArchive)
        self.movefromlcu(self.stationcontroller.lcuDumpDir+"/*00[XY].dat",
                         datapath)
        # beamlet statistics also generate empty *01[XY].dat so remove:
        self.stationcontroller.rm(
                              self.stationcontroller.lcuDumpDir+"/*01[XY].dat")
        obsinfo.create_LOFARst_header(datapath)
        return datapath

    def do_sst(self, freqbndobj, integration, duration, pointing, obsinfo):
        """Run an sst static."""

        # Use this to do a system temperature measurement.
        SYS_TEMP_MEAS = False

        # Start beamforming
        rcu_setup_CMD, beamctl_cmds = \
            self.streambeams(freqbndobj, pointing)
        if SYS_TEMP_MEAS:
            lbbalnaoff_CMD = self.stationcontroller.turnoffLBA_LNAs()
            #beamctl_CMD += "\n"+lbbalnaoff_CMD

        # Get some metadata about operational settings:
        caltabinfo = ""    # No need for caltab info

        # Record data
        rspctl_CMD = self.stationcontroller.rec_sst(integration, duration)

        # Move data to archive
        obsdatetime_stamp = self.get_data_timestamp()
        obsinfo.setobsinfo_fromparams('sst', obsdatetime_stamp, beamctl_cmds, rspctl_CMD,
                                      caltabinfo)
        bsxSTobsEpoch, datapath = obsinfo.getobsdatapath(self.LOFARdataArchive)
        self.movefromlcu(self.stationcontroller.lcuDumpDir+"/*.dat", datapath,
                         recursive=True)
        obsinfo.create_LOFARst_header(datapath)
        return datapath

    def do_xst(self, freqbndobj, integration, duration, pointing, obsinfo):
        """Run an xst statistic towards the given pointing. This corresponds to
        a crosscorrelation of all elements at the given frequency and
        integration repeated for a duration of seconds."""

        caltabinfo = ""  # No need for caltab info for xst data
        obsinfolist = []
        # Get subbands to do
        for sb_rcumode in freqbndobj.sb_range:
            # Start beamforming
            rcu_setup_CMD, beamctl_cmds = self.streambeams(freqbndobj, pointing)
            if ':' in sb_rcumode:
                sblo, sbhi = sb_rcumode.split(':')
                subbands = range(int(sblo),int(sbhi)+1)
            else:
                subbands = [int(sb) for sb in sb_rcumode.split(',')]
            for subband in subbands:
                # Record data
                rspctl_CMD = self.stationcontroller.rec_xst(subband, integration,
                                                            duration)
                obsdatetime_stamp = self.get_data_timestamp(-1)
                curr_obsinfo = copy.copy(obsinfo)
                curr_obsinfo.setobsinfo_fromparams('xst', obsdatetime_stamp, beamctl_cmds,
                                                   rspctl_CMD, caltabinfo)
                obsinfolist.append(curr_obsinfo)
            self.stationcontroller.stopBeam()

        obsinfo = copy.copy(obsinfolist[0])
        obsinfo.sb = freqbndobj.sb_range[0]
        # Move data to archive
        bsxSTobsEpoch, datapath = obsinfo.getobsdatapath(self.LOFARdataArchive)

        self.movefromlcu(self.stationcontroller.lcuDumpDir+"/*.dat", datapath)
        for curr_obsinfo in obsinfolist:
            curr_obsinfo.create_LOFARst_header(datapath)
        return datapath

    def bsxST(self, statistic, freqbndobj, integration, duration, pointSrc):
        """Run a statisics observation.

        Parameters
        ----------
        statistic : str
            Statistic can be either 'BST', 'SST' or 'XST'
            (corresponding to spectral, beamlet or crosscorrelation).
        frequency :
            frequency in Hz.
        integration : int
            integration time in seconds.
        duration : int
            duration of observation in seconds.
        pointSrc : str
            point direction as a beamctl triplet.
        """
        if not self.checkobservingallowed():
            raise RuntimeError, "Station is being used by someone else."

        try:
            pointing = stdPointings(pointSrc)
        except KeyError:
            try:
                phi, theta, ref = pointSrc.split(',', 3)
                # FIXME:  (not always going to be correct)
                pointing = pointSrc
            except ValueError:
                raise ValueError, "Error: %s invalid pointing syntax".format(pointSrc)
        if integration > duration:
            raise ValueError, "Integration {} is longer than duration {}.\
                               ".format(integration, duration)

        obsinfo = dataIO.ObsInfo(self.stationcontroller.stnid, self.project,
                                 self.observer)
        datapath = None
        if statistic == 'bst':
             datapath = self.do_bst(freqbndobj, integration, duration, pointing,
                                    obsinfo)
        elif statistic == 'sst':
             datapath = self.do_sst(freqbndobj, integration, duration, pointing,
                                    obsinfo)
        elif statistic == 'xst':
             datapath = self.do_xst(freqbndobj, integration, duration, pointing,
                                    obsinfo)

        self.stationcontroller.stopBeam()
        return datapath

# End: Basic obs modes
#######################################

    def do_acc(self, band, duration_req, pointSrc='Z', exit_obsstate=True):
        """Perform calibration observation mode on station. Also known as ACC
        mode. The duration may be longer than requested so as to fit within the
        cadence of whole ACC aquisitions (512+7=519 seconds). swlevel needs to
        cycle down to 2 (or less) and then to 3. If swlevel is kept at 3
        (i.e. exit_obsstate=False), then ACC will continue to be produced, until
        swlevel goes below 2.

        ACC files are autocovariance-cubes: the covariance of all array
        elements with each as a function of subband. These files are generated
        by the MAC service called CalServer. It run at swlevel 3 and is
        configured in the file lofar/etc/CalServer.conf. Note subband
        integration is always 1s, so ACC file is dumped after 512 seconds.

        Parameters
        ----------
            band: str
            duration_req: int
            pointSrc: str
        """
        if not self.checkobservingallowed():
            raise RuntimeError
        try:
            rcumode = stationcontrol.band2rcumode(band)
        except ValueError:
            raise
        try:
            pointing = stdPointings(pointSrc)
        except KeyError:
            try:
                phi, theta, ref = pointSrc.split(',', 3)
            except ValueError:
                print("Error: %s invalid pointing syntax".format(pointSrc))
            pointing = pointSrc
        # Get timings
        # Also duration of ACC sweep since each sb is 1 second.
        nrACCsbs = stationcontrol.TotNrOfsb
        # Time between end of one ACC sweep and beginning of next one.
        timebetweenACCs = 7
        ACCcadence = float(nrACCsbs+timebetweenACCs) # =519s time between two ACCs
        duration = int(math.ceil((duration_req-nrACCsbs)/ACCcadence)
                       * (ACCcadence)+nrACCsbs+timebetweenACCs)
        if duration != duration_req:
            print("""Warning: using longer duration {}s to fit with ACC
                  cadence.""".format(duration))
        obsStartDate = time.strftime("%Y%m%d_%H%M%S", time.gmtime())
        sst_integration = int(ACCcadence)

        # Run ACC mode
        #beamctl_CMD, rspctl_CMD = \
        #    self.stationcontroller.runACC(rcumode, duration, pointing)

        # Make sure swlevel=<2
        self.stationcontroller.bootToObservationState(2)

        # Set CalServ.conf to dump ACCs:
        self.stationcontroller.acc_mode(enable=True)

        # Boot to swlevel 3 so the calserver service starts (
        self.stationcontroller.bootToObservationState()

        # Beamlet & Subband allocation does not matter here
        # since niether ACC or SST cares
        beamctl_CMD = self.stationcontroller.runbeamctl('0', '255', rcumode, pointing)

        # Run for $duration seconds
        rspctl_CMD = self.stationcontroller.rec_sst(sst_integration, duration)
        self.stationcontroller.stopBeam()

        # Switch back to normal state i.e. turn-off ACC dumping:
        self.stationcontroller.acc_mode(enable=False)

        if exit_obsstate:
            self.stationcontroller.bootToObservationState(0)

        # Transfer data from LCU to DAU
        ACCsrcFiles = self.stationcontroller.ACCsrcDir+"/*.dat"
        ACCdestDir = \
            os.path.join(self.LOFARdataArchive, 'acc',
                       '{}_{}_rcu{}_dur{}'.format(self.stationcontroller.stnid,
                         obsStartDate, rcumode, duration))
        if int(rcumode) > 3:
            ACCdestDir += "_"+pointSrc
        ACCdestDir += "_acc"
        if os.path.isdir(ACCdestDir):
            print("Appropriate directory exists already (will put data here)")
        else:
            print("Creating directory "+ACCdestDir+" for ACC "+str(duration)\
                  + " s rcumode="+rcumode+" calibration")
            os.mkdir(ACCdestDir)

        # Move ACC dumps to storage
        obsdatetime_stamp = self.get_data_timestamp()
        self.movefromlcu(ACCsrcFiles, ACCdestDir)

        # Move concurrent data to storage
        obsinfo = dataIO.ObsInfo(self.stationcontroller.stnid, self.project, self.observer)
        obsinfo.setobsinfo_fromparams('sst', obsdatetime_stamp, beamctl_CMD, rspctl_CMD)
        bsxSTobsEpoch, datapath = obsinfo.getobsdatapath(self.LOFARdataArchive)
        self.movefromlcu(self.stationcontroller.lcuDumpDir+"/*", datapath,
                         recursive=True)
        obsinfo.create_LOFARst_header(datapath)
        self.stationcontroller.cleanup()

        # Postprocess?
        if False:
            if not self.stationcontroller.DryRun and not rcumode == '3':
                # Translate the incoming data into hdf or ms:
                import acc2bst
                acc2bst.main(ACCdestDir, self.stationcontroller.stnid,
                             pointSrc)
            else:
                print("Not reducing ACCs on DPU in {} (stnid={}) into BSTs...\
                      ".format(ACCdestDir, self.stationcontroller.stnid))
        return ACCdestDir

    def do_SEPTON(self, statistic,  frqbndobj, integration, duration, elemsOn=elOn_gILT):
        """Record xst or sst data in SEPTON mode.
        Setup Single Element per Tile ON mode. This only valid for HBA and
        currently only rcumode=5."""
        #subband, NqZone = stationcontrol.freq2sb(frequency)
        #rcumode = stationcontrol.NyquistZone2rcumode(NqZone)
        subband = frqbndobj.sb_range[0]
        rcumode = frqbndobj.rcumodes[0]
        pointing = ""
        rcusetup_CMD =""
        rspctl_SET = ""
        beamctl_CMD = ""
        # NOTE: LCU must be in swlevel=2 to run SEPTON!
        self.stationcontroller.bootToObservationState(2)
        # self.stationcontroller.turnoffElinTile_byTile(elemsOn) # Alternative
        self.stationcontroller.turnoffElinTile_byEl(elemsOn)
        caltabinfo = ""  # No need for caltab info
        # Record data
        if statistic == 'xst':
            rspctl_CMD = self.stationcontroller.rec_xst(subband, integration,
                                                    duration)
        elif statistic == 'sst':
            rspctl_CMD = self.stationcontroller.rec_sst(integration, duration)
        else:
            raise ValueError("Only xst or sst data allowed.")
        LOFARdatTYPE = statistic + "-SEPTON"
        # Collect observational metadata
        obsdatetime_stamp = self.get_data_timestamp()

        obsinfo = dataIO.ObsInfo(self.stationcontroller.stnid, self.project, self.observer)
        obsinfo.setobsinfo_fromparams(LOFARdatTYPE, obsdatetime_stamp, beamctl_CMD, rspctl_CMD,
                                      caltabinfo, septonconf=elementMap2str(elemsOn))
        # Move data to archive
        bsxSTobsEpoch, datapath = obsinfo.getobsdatapath(self.LOFARdataArchive)
        self.movefromlcu(self.stationcontroller.lcuDumpDir+"/*.dat", datapath)
        obsinfo.create_LOFARst_header(datapath)
        return (bsxSTobsEpoch, rspctl_SET, beamctl_CMD, rspctl_CMD, caltabinfo, datapath)

    #####################
    # BEGIN: TBB services

    def do_tbb(self, duration, band):
        """Record duration seconds of TBB data from rcumode."""

        observer = self.observer
        project = self.project
        observationID = "Null"
        # Start a dummy beam
        pointing = stdPointings('Z')
        freqBand = stationcontrol.band2freqrange(band)
        freqmid = (freqBand[0]+freqBand[1])/2.0
        antset = stationcontrol.band2antset(band)
        self.streambeams(freqmid, pointing)

        print "Set up TBBs"
        self.stationcontroller.setupTBBs()

        print "Wait a while for TBBs to fill up"
        time.sleep(20)
        # Start data capture process locally
        dalcap = \
            multiprocessing.Process(target=capture_data_DAL1,
                                    args=(self.tbbraw2h5cmd, self.tbbh5dumpdir,
                                          observer, antset, project,
                                          observationID,False))
        dalcap.start()
        time.sleep(20)  # Arbitrary time to trigger
        print "Send trigger to TBBs"
        self.stationcontroller.freezeTBBdata()
        print "Start streaming "+str(duration)+" s of TBB data out of LCU"
        self.stationcontroller.startTBBdataStream(float(duration))
        dalcap.join()

    # END: TBB services
    ###################

# END: Session
##############

# TBBh5dumpDir="/home/tobia/lofar/data/tbb/h5/"
# "/mnt/lane0/TBB/" #Should end with "/" character.
# h5fileprefix="L"

def capture_data_DAL1(tbbraw2h5cmd, TBBh5dumpDir, observer, antennaSet,
                      project, observationID, background=False):
    """Start process on DPU to capture streamed TBB data on LCU."""
    if background:
        ground = 'bg'
        grdcmd = '&'
    else:
        ground = 'fg'
        grdcmd = ''
    print("Starting TBBraw2h5 ({})".format(ground))
    subprocess.call(("cd "+TBBh5dumpDir+"; "
                     + " export LD_LIBRARY_PATH=/mnt/old/usr/lib/ ; "
                     + tbbraw2h5cmd
                     + " "
                     + " --observer="+observer
                     + " --antennaSet="+antennaSet
                     + " --project="+project
                     + " --observationID="+observationID
                     + " -P 31664 -P 31665 -P 31666 "
                     + " -P 31667 -P 31668 -P 31669 "
                     + " -P 31670 -P 31671 -P 31672 "
                     + " -P 31673 -P 31674 -P 31675 "
                     + grdcmd), shell=True)
