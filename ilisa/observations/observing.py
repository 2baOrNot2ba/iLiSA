#!/usr/bin/python
"""Package that provides functions for setting up and running observations
via station controller objects. This package knows about the data archive and
should not run anything directly on LCU."""


import math
import time
import subprocess
import os
import numpy
import yaml
import multiprocessing
import ilisa.observations.stationcontrol as stationcontrol
import ilisa.observations.dataIO as dataIO


# SEPTON configurations:
#        1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24
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


# Tile parameter functions
def elementMap2str(elmap):
    elmapStr = ""
    for el in elmap:
        elmapStr = elmapStr+hex(el).lstrip('0').lstrip('x')
    return elmapStr


#######################################
# Begin Direction related code

def pointingGrid(NrAzDirs=8, NrElDirs=7):
    """Returns a tuple of LOFAR beamctl directions strings of a spherical
    grid of pointings around zenith."""
    # Nr of pointings NrAzDirs*NrElDirs+1 (where 1 is zenith)
    # should be less than 61 if they are to fit in one lane.
    az = numpy.linspace(0, 2*math.pi, NrAzDirs+1)
    numpy.delete(az, NrAzDirs)
    az = az-0*math.pi/4
    el = numpy.linspace(0, math.pi/2, NrElDirs+1)
    numpy.delete(el, NrElDirs)
    pntGridStrs = []
    # El Major:
    for ElDirNr in range(NrElDirs):
        for AzDirNr in range(NrAzDirs):
            nextPnting = str(az[AzDirNr])+","+str(el[ElDirNr])+",AZELGEO"
            pntGridStrs.append(nextPnting)
    zenith = '0.0,'+str(math.pi/2)+',AZELGEO'
    pntGridStrs.append(zenith)
    return tuple(pntGridStrs)


def parsebeamctldir(beamctldirarg):
    """Parse a beamctl direction string into direction tuple.

    Parameters
    ----------
    beamctldirarg : str
        String with format 'angle1,angle2,refsys'

    Returns
    -------
    dirtuple : tuple or None
        Direction tuple defined as (angle1: float, angle2: float, refsys: str)
        or None if beamctldirarg was not correct format.
    """
    try:
        angle1str, angle2str, refsys = beamctldirarg.split(',')
        angle1, angle2 = float(angle1str), float(angle2str)
        # TODO should check that refsys is one of the valid refsys strings.
        dirtuple = (angle1, angle2, refsys)
        return dirtuple
    except ValueError:
        return None


def stdPointings(directionterm):
    """Find beamctl direction string based on direction term.

    Parameters
    ----------
    directionterm : str or None
        Source name or term for a direction. E.g. 'Z' for Zenith
        or 'CasA' for Cassiopeia A. If set to None, the function returns the
        direction terms it knows about.

    Returns
    -------
    beamctldir : str
        Argument suitable for beamctl direction arguments --anadir and
        --digdir. If input is None it will return all the direction terms it
        knows.

    Raises
    ------
    KeyError
        Thrown if directionterm is not understood.
    """
    term2beamstr = {  # 1e-6 rad < 1arcsec
          'N':    str(0*math.pi/2)+",0.,AZELGEO",
          'E':    str(1*math.pi/2)+",0.,AZELGEO",
          'S':    str(2*math.pi/2)+",0.,AZELGEO",
          'W':    str(3*math.pi/2)+",0.,AZELGEO",
          'Z':    '0.,'+str(math.pi/2)+',AZELGEO',
          'CasA': '6.123487,1.026515,J2000',
          'CygA': '5.233660,0.710940,J2000',
          'TauA': '1.459672,0.384225,J2000',
          'VirA': '3.276086,0.216265,J2000',
          'Sun':  '0.,0.,SUN',
          'Jupiter': '0.,0.,JUPITER',
          'Moon': '0.,0.,MOON',
          'PSR_LGM': '5.0691,0.3819,J2000',
          'NCP':  '0.,'+str(math.pi/2)+',ITRF'
    }
    if directionterm is None:
        return term2beamstr.keys()
    if directionterm in term2beamstr:
        return term2beamstr[directionterm]
    else:
        raise KeyError('Requested source {} unknown.'.format(directionterm))


def normalizebeamctldir(gendirstr):
    """Parse a general direction string.

    Parameters
    ----------
    gendirstr : str
        This could be one of the direction terms or it could a beamctl
        direction string.

    Returns
    -------
    beamctldirstr : str
        The beamctl direction string corresponding to input gendirstr.
    """
    beamctldir = parsebeamctldir(gendirstr)
    if beamctldir is None:
        try:
            beamctldirstr = stdPointings(gendirstr)
        except:
            raise ValueError('General direction term {} unknown.'.format(gendirstr))
    else:
        beamctldirstr = gendirstr
    return beamctldirstr

# End Direction related code
#######################################


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
        self.bits = 16  # Default to 16
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

    def get_data_timestamp(self):
        """Get timestamp of datafiles on LCU."""
        dd_dir, acc_dir = self.stationcontroller.getdatalist()
        # Assumes only one file in datadump dir with
        # format YYYYmmdd_HHMMSS_[bsx]st.dat
        obsdate, obstime, obssuff = dd_dir[0].split('_', 2)
        obsdatetime_stamp = obsdate+'_'+obstime
        return obsdatetime_stamp

#######################################
# Begin: Basic obs modes
    def getbsxstdatapath(self, LOFARdatTYPE, datetime_stamp, rcumode, sb,
                         integration, duration, pointing):
            """Create name and destination path for folders (on the DPU) in
            which to save the various LOFAR data products.
            """
            stDataArchive = os.path.join(self.LOFARdataArchive, LOFARdatTYPE)
            stObsEpoch = datetime_stamp
            st_extName = stObsEpoch+"_rcu"+str(rcumode)
            if str(sb) != "":
                st_extName += "_sb"+str(sb)
            st_extName += "_int"+str(integration)+"_dur"+str(duration)
            if str(pointing) != "":
                st_extName += "_dir"+str(pointing)
            st_extName += "_"+LOFARdatTYPE
            datapath = os.path.join(stDataArchive, st_extName)
            return stObsEpoch, datapath

    def do_bst(self, frequency, integration, duration, pointing):
        """Do a Beamlet STatistic (bst) recording on station. frequency in Hz.
        """
        if not self.checkobservingallowed():
            raise RuntimeError

        CALTABLESRC = 'default'   # FIXME put this in args
        freqLo = frequency
        # Setup Calibration:
        ## (Only BST uses calibration tables)
        # Choose between 'default' or 'local'
        self.stationcontroller.selectCalTable(CALTABLESRC)
        # Start beamforming
        (beamctl_CMD, rspctl_SET, beamctl_main
         ) = self.stationcontroller.streambeam((freqLo,),
                                               pointing, bits=self.bits,
                                               DUMMYWARMUP=True)
        rcusetup_CMD = ""
        # rcusetup_CMD = self.stationcontroller.rcusetup(bits, attenuation)
        # Get some metadata about operational settings:
        (antset, rcus, rcumode, beamlets, subband, anadir, digdir
         ) = stationcontrol.parse_beamctl_args(beamctl_main)
        caltabinfo = self.stationcontroller.getCalTableInfo(rcumode)
        # Record data
        # waittime = 0
        # print "Waiting extra", str(waittime) +" seconds" #Seems necessary
        # time.sleep(waittime)
        rspctl_CMD = self.stationcontroller.rec_bst(integration, duration)
        # Move data to archive
        obsdatetime_stamp = self.get_data_timestamp()
        bsxSTobsEpoch, datapath = \
            self.getbsxstdatapath('bst', obsdatetime_stamp, rcumode, subband,
                                  integration, duration, pointing)
        self.movefromlcu(self.stationcontroller.lcuDumpDir+"/*00[XY].dat",
                         datapath)
        # beamlet statistics also generate empty *01[XY].dat so remove:
        self.stationcontroller.rm(
                              self.stationcontroller.lcuDumpDir+"/*01[XY].dat")
        return (bsxSTobsEpoch, rcusetup_CMD, beamctl_CMD, rspctl_CMD,
                caltabinfo, datapath)

    def do_sst(self, frequency, integration, duration, pointing):
        """Run an sst static."""

        if not self.checkobservingallowed():
            raise RuntimeError
        # Use this to do a system temperature measurement.
        SYS_TEMP_MEAS = False

        # Specify freq range
        freqLo = frequency
        # Start beamforming
        (beamctl_CMD, rspctl_SET, beamctl_main) = \
            self.stationcontroller.streambeam((freqLo,), pointing)
        if SYS_TEMP_MEAS:
            lbbalnaoff_CMD = self.stationcontroller.turnoffLBA_LNAs()
            beamctl_CMD += "\n"+lbbalnaoff_CMD
        # rcusetup_CMD=""
        # rcusetup_CMD=self.stationcontroller.rcusetup(bits, attenuation)
        # Get some metadata about operational settings:
        (antset, rcus, rcumode, beamlets, subband, anadir, digdir
         ) = stationcontrol.parse_beamctl_args(beamctl_main)
        caltabinfo = ""    # No need for caltab info
        # Record data
        rspctl_CMD = self.stationcontroller.rec_sst(integration, duration)
        # Move data to archive
        obsdatetime_stamp = self.get_data_timestamp()
        bsxSTobsEpoch, datapath = \
            self.getbsxstdatapath('sst', obsdatetime_stamp, rcumode, "",
                                  integration, duration, "")
        self.movefromlcu(self.stationcontroller.lcuDumpDir+"/*.dat", datapath,
                         recursive=True)
        return (bsxSTobsEpoch, rspctl_SET, beamctl_CMD, rspctl_CMD, caltabinfo,
                datapath)

    def do_xst(self, frequency, integration, duration, pointing):
        """Run an xst statistic towards the given pointing. This corresponds to
        a crosscorrelation of all elements at the given frequency and
        integration repeated for a duration of seconds."""

        if not self.checkobservingallowed():
            raise RuntimeError

        # Start beamforming
        beamctl_CMDs, rspctl_SET, beamctl_main = \
            self.stationcontroller.streambeam(frequency, pointing)
        # FIX: Include rcusetup into streambeam call.
        # rcusetup_CMD = self.stationcontroller.rcusetup(bits, attenuation)
        # Get some metadata about operational settings:
        (antset, rcus, rcumode, beamlets, subband, anadir, digdir
         ) = stationcontrol.parse_beamctl_args(beamctl_main)
        caltabinfo = ""  # No need for caltab info
        # Record data
        rspctl_CMD = self.stationcontroller.rec_xst(subband, integration,
                                                    duration)
        # Move data to archive
        obsdatetime_stamp = self.get_data_timestamp()
        bsxSTobsEpoch, datapath = \
            self.getbsxstdatapath('xst', obsdatetime_stamp, rcumode, subband,
                                  integration, duration, pointing)

        self.movefromlcu(self.stationcontroller.lcuDumpDir+"/*.dat", datapath)
        return (bsxSTobsEpoch, rspctl_SET, beamctl_CMDs, rspctl_CMD,
                caltabinfo, datapath)

    def bsxST(self, statistic, frequency, integration, duration, pointSrc):
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

        if statistic == 'bst':
            (bsxSTobsEpoch, rcusetup_CMD, beamctl_CMD, rspctl_CMD, caltabinfo,
             datapath) = self.do_bst(frequency, integration, duration, pointing
                                     )
        elif statistic == 'sst':
            (bsxSTobsEpoch, rcusetup_CMD, beamctl_CMD, rspctl_CMD, caltabinfo,
             datapath) = self.do_sst(frequency, integration, duration, pointing
                                     )
        elif statistic == 'xst':
            (bsxSTobsEpoch, rcusetup_CMD, beamctl_CMD, rspctl_CMD, caltabinfo,
             datapath) = self.do_xst(frequency, integration, duration, pointing
                                     )
        self.create_LOFARst_header(statistic, datapath, bsxSTobsEpoch,
                                   rcusetup_CMD, beamctl_CMD, rspctl_CMD,
                                   caltabinfo)
        self.stationcontroller.stopBeam()
        return datapath

# End: Basic obs modes
#######################################

    def do_acc(self, band, duration_req, pointSrc='Z', exit_obsstate=False):
        """Perform calibration observation mode on station. Also known as ACC
        mode. The duration may be longer than requested so as to fit within the
        cadence of whole ACC aquisitions (512+7=519 seconds). swlevel needs to
        cycle down to 2 (or less) and then to 3.

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
        timeintervalbetweenACCs = 7
        ACCcadence = float(nrACCsbs+timeintervalbetweenACCs)
        endbuftime = timeintervalbetweenACCs
        duration = int(math.ceil((duration_req-nrACCsbs)/ACCcadence)
                       * (ACCcadence)+nrACCsbs+endbuftime)
        if duration != duration_req:
            print("Warning: will use longer duration {}s to fit with ACC\
                  cadence.".format(duration))
        obsStartDate = time.strftime("%Y%m%d_%H%M%S", time.gmtime())

        # Run ACC mode
        sst_integration = 600
        beamctl_CMD, rspctl_CMD = \
            self.stationcontroller.runACC(rcumode, duration, pointing,
                                          sst_integration)
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
            print "Appropriate directory exists already (will put data here)"
        else:
            print "Creating directory "+ACCdestDir+" for ACC "+str(duration)\
                  + " s rcumode="+rcumode+" calibration"
            os.mkdir(ACCdestDir)

        # Move ACC dumps to storage
        obsdatetime_stamp = self.get_data_timestamp()
        self.movefromlcu(ACCsrcFiles, ACCdestDir)

        # Move concurrent data to storage
        bsxSTobsEpoch, datapath = \
            self.getbsxstdatapath('sst', obsdatetime_stamp,
                                  rcumode, "", sst_integration, duration, "")
        self.movefromlcu(self.stationcontroller.lcuDumpDir+"/*", datapath,
                         recursive=True)
        self.create_LOFARst_header('sst', datapath, bsxSTobsEpoch, "",
                                   beamctl_CMD, rspctl_CMD, "")
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

    def setupSEPTON(self, elemsOn=elOn_gILT):
        """Setup Single Element per Tile ON mode. This only valid for HBA and
        currently only rcumode=5."""
        # NOTE: LCU must be in swlevel=2 to run SEPTON!
        self.stationcontroller.bootToObservationState(2)
        # self.stationcontroller.turnoffElinTile_byTile(elemsOn) # Alternative
        self.stationcontroller.turnoffElinTile_byEl(elemsOn)

    def do_SEPTON(self, subband, integration, duration, elemsOn=elOn_gILT, statistic='xst'):
        """Record xst or sst data in SEPTON mode."""
        rcumode = 5
        pointing = ""
        rcusetup_CMD =""
        rspctl_SET = ""
        beamctl_CMD = ""
        self.setupSEPTON(elemsOn)
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
        obsinfo.setobsinfo(LOFARdatTYPE, obsdatetime_stamp, rcumode, subband, integration,
                           duration, pointing
                           )
        # Move data to archive
        bsxSTobsEpoch, datapath = obsinfo.getobsdatapath(self.LOFARdataArchive)
        self.movefromlcu(self.stationcontroller.lcuDumpDir+"/*.dat", datapath)
        obsinfo.create_LOFARst_header(statistic, datapath, bsxSTobsEpoch, rcusetup_CMD,
                                      beamctl_CMD, rspctl_CMD, caltabinfo,
                                      septonconfig=elementMap2str(elemsOn))
        return (bsxSTobsEpoch, rspctl_SET, beamctl_CMD, rspctl_CMD, caltabinfo, datapath)

    #####################
    # BEGIN: TBB services

    def do_tbb(self, duration, band):
        """Record duration seconds of TBB data from rcumode."""

        observer = "TobiaC"
        project = "LOCAL"
        observationID = "Null"
        # Start a dummy beam
        pointing = stdPointings('Z')
        freqBand = stationcontrol.band2freqrange(band)
        freqmid = (freqBand[0]+freqBand[1])/2.0
        antset = stationcontrol.band2antset(band)
        self.stationcontroller.streambeam(freqmid, pointing)

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

    def create_LOFARst_header(self, LOFARstTYPE, datapath, LOFARstObsEpoch,
                              rspsetup_CMD, beamctl_CMD, rspctl_CMD,
                              caltableInfo=""):
        """Create a header file for LOFAR standalone observation."""
        def indenttext(txt):
            indentstr = "  "
            return indentstr+txt.replace("\n","\n"+indentstr)
        headerversion = "2"
        if (LOFARstTYPE != 'bst' and LOFARstTYPE != 'sst'
                and LOFARstTYPE != 'xst' and LOFARstTYPE != 'bf'):
            raise ValueError, "Unknown LOFAR statistic type {}.\
                              ".format(LOFARstTYPE)
        LOFARstHeaderFile = LOFARstObsEpoch+"_"+LOFARstTYPE+".h"
        f = open(os.path.join(datapath, LOFARstHeaderFile), "w")
        f.write("# HeaderType: bsxSTdata (YAML)\n")
        f.write("# Header version {}\n".format(headerversion))
        f.write("Observer: {}\n".format(self.observer))
        f.write("Project: {}\n".format(self.project))
        f.write("DataType: {}\n".format(LOFARstTYPE))
        f.write("StationID: "+self.stationcontroller.stnid+"\n")
        starttime = LOFARstObsEpoch[0:4]+'-'+LOFARstObsEpoch[4:6]+'-'\
                        + LOFARstObsEpoch[6:8]+'T'+LOFARstObsEpoch[9:11]+':'\
                        + LOFARstObsEpoch[11:13]+':'+LOFARstObsEpoch[13:15]
        f.write("StartTime: "+starttime+"\n")
        f.write("BeamctlCmds: |-\n")
        f.write(indenttext(beamctl_CMD)+"\n")
        # f.write(rspsetup_CMD+"\n")
        # FIX separation of beamctl and rspsetup
        # (Currently rspsetup is in beamctl)
        f.write("RspctlCmds: |-\n")
        f.write(indenttext(rspctl_CMD)+"\n")
        if LOFARstTYPE == 'bst':
            f.write("CalTabInfo: |-\n")
            #f.write(indenttext(caltableInfo))
            f.write(str(caltableInfo))
        f.close()


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
