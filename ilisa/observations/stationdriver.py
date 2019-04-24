#!/usr/bin/python
"""Package that provides functions for setting up and running observations
via station controller objects. This package knows about the data archive and
should not run anything directly on LCU."""


import time
import datetime
import subprocess
import os
import multiprocessing
import copy

import ilisa.observations.lcuinterface as stationcontrol
import ilisa.observations.modeparms as modeparms
import ilisa.observations.beamformedstreams.bfbackend as bfbackend
import ilisa.observations.programs as programs


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
# elOn_same = [elOn_same_el for elemNr in range(stationcontrol.nrTiles)]
elemsOn = elOn_Generic_Int_201512  # elOn_same or elOn_step or elOn_gILT or ...


class StationDriver(object):
    """observing.Session class is a client type class typically running on the
    compute node."""

    def cleanup(self):
        pass

    def checkobservingallowed(self):
        """Check whether observations are allowed. This occurs is someone else
        is using the station."""
        serviceuser = self.lcu_interface.who_servicebroker()

        if serviceuser == 'None' or serviceuser == self.lcu_interface.user:
            stationswitchmode = self.lcu_interface.getstationswitchmode()
            if stationswitchmode == 'local':
                return True
            else:
                print("Warning: Station is not in stand-alone mode.")
                return False
        else:
            print("Warning: Someone else ({}) is using LCU".format(serviceuser))
            print("         (You are running as {})".format(self.lcu_interface.user))
            return False

    def __init__(self, accessconf, stnsesinfo,
                 goto_observingstate_when_starting=True):
        """Initialize a StationDriver object, which has access to a station via
        a LCUInterface object configured with setting given by accessconf dict.
        When goto_observingstate_when_starting is True, boot the
        station to swlevel 3, but only if observations are allowed on the
        station.
        """
        lcuaccessconf = accessconf['LCU']
        dpuaccessconf = accessconf['DPU']
        self.lcu_interface = stationcontrol.LCUInterface(lcuaccessconf)
        self.stnsesinfo = stnsesinfo
        self.stnsesinfo.set_stnid(self.lcu_interface.stnid)

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
        if goto_observingstate_when_starting:
            try:
                self.goto_observingstate()
            except RuntimeError:
                raise RuntimeError('Observations not allowed on this station')

    def goto_observingstate(self):
        """Put station into the main observing state."""
        if not self.checkobservingallowed():
            raise RuntimeError('Observations not allowed')
        self.lcu_interface.set_swlevel(3)

    def haltobservingstate(self):
        """Halt observing state on station."""
        if self.checkobservingallowed():
            self.lcu_interface.set_swlevel(0)
            self.cleanup()
            # Cleanup any data left on LCU.
            self.lcu_interface.cleanup()
        # Close down stationcontroller:
        del self.lcu_interface

    def __del__(self):
        """Normally, shutdown observation mode on station."""
        if self.halt_observingstate_when_finished:
            self.haltobservingstate()
        elif self.exit_check:
            if self.checkobservingallowed():
                swlevel = self.lcu_interface.get_swlevel()
                if int(swlevel) != 0:
                    print("Warning: You are leaving station in swlevel {} != 0"
                          .format(swlevel))

    def movefromlcu(self, source, dest, recursive=False):
        """Move file(s) off LCU to DPU."""
        if not os.path.isdir(dest):
            os.mkdir(dest)
        movecmd = "scp"
        if recursive:
            movecmd += " -r"
        fullcmd = movecmd +" " + self.lcu_interface.lcuURL + ":" + source + " " + dest
        cmdprompt = "on DPU>"
        if self.lcu_interface.DryRun:
            cmdprompt = "(dryrun) "+cmdprompt
        if self.lcu_interface.verbose:
            print("{} {}".format(cmdprompt, fullcmd))
        if not self.lcu_interface.DryRun:
            subprocess.call(fullcmd, shell=True)
            self.lcu_interface.rm(source)

    def get_data_timestamp(self, order=0, ACC=False):
        """Get timestamp of datafiles on LCU. order is the temporal order of the data
        file, order=0 is oldest, order=-1 is newest."""
        dd_dir, acc_dir = self.lcu_interface.getdatalist()
        # Assumes only one file in datadump dir with
        # format YYYYmmdd_HHMMSS_[bsx]st.dat
        if not ACC:
             the_dir = dd_dir
        else:
            the_dir = acc_dir
        obsdate, obstime, obssuff = the_dir[order].split('_', 2)
        obsdatetime_stamp = obsdate+'_'+obstime
        return obsdatetime_stamp

    def get_stnid(self):
        """Return the station id that this StationDriver is managing."""
        return self.lcu_interface.stnid

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
            beamctl_main = self.lcu_interface.run_beamctl(beamletIDs, subbands,
                                                          rcumode, pointing, rcusel)
            beamctl_cmds.append(beamctl_main)
        rcu_setup_CMD = self.lcu_interface.rcusetup(bits, attenuation)
        return rcu_setup_CMD, beamctl_cmds

    def _waittoboot(self, starttimestr, pause):
        """Before booting, wait until time given by starttimestr. This includes
         a dummy beam warmup."""
        nw = datetime.datetime.utcnow()
        st = datetime.datetime.strptime(starttimestr, "%Y-%m-%dT%H:%M:%S")

        maxminsetuptime = datetime.timedelta(seconds=105 + pause)  # Longest minimal time
        # before observation
        # start to set up
        d = (st - maxminsetuptime) - nw
        timeuntilboot = d.total_seconds()
        if timeuntilboot < 0.:
            timeuntilboot = 0
        print("Will boot to observe state after " + str(timeuntilboot) + " seconds...")
        time.sleep(timeuntilboot)
        return st

    def do_SEPTON(self, statistic,  frqbndobj, integration, duration_scan,
                  elemsOn=elOn_gILT):
        """Record xst or sst data in SEPTON mode.
        Setup Single Element per Tile ON mode. This only valid for HBA and
        currently only rcumode=5."""
        try:
            self.goto_observingstate()
        except RuntimeError as e:
            raise RuntimeError(e)
        subband = frqbndobj.sb_range[0]
        rcumode = frqbndobj.rcumodes[0]
        rspctl_SET = ""
        beamctl_CMD = ""
        # NOTE: LCU must be in swlevel=2 to run SEPTON!
        self.lcu_interface.set_swlevel(2)
        # self.stationcontroller.turnoffElinTile_byTile(elemsOn) # Alternative
        self.lcu_interface.turnoffElinTile_byEl(elemsOn)
        caltabinfo = ""  # No need for caltab info
        # Record data
        if statistic == 'xst':
            rspctl_CMD = self.lcu_interface.rec_xst(subband, integration, duration_scan)
        elif statistic == 'sst':
            rspctl_CMD = self.lcu_interface.rec_sst(integration, duration_scan)
        else:
            raise ValueError("Only xst or sst data allowed.")
        LOFARdatTYPE = statistic + "-SEPTON"
        # Collect observational metadata
        obsdatetime_stamp = self.get_data_timestamp()

        stnsesinfo = copy.deepcopy(self.stnsesinfo)
        stnsesinfo.new_obsinfo()
        stnsesinfo.obsinfos[-1].setobsinfo_fromparams(
            LOFARdatTYPE, obsdatetime_stamp, beamctl_CMD, rspctl_CMD, caltabinfo,
            septonconf = modeparms.elementMap2str(elemsOn))

        # Move data to archive
        bsxSTobsEpoch, datapath = stnsesinfo.obsinfos[-1].getobsdatapath(
            self.LOFARdataArchive)
        self.movefromlcu(self.lcu_interface.lcuDumpDir + "/*.dat", datapath)
        stnsesinfo.obsinfos[-1].create_LOFARst_header(datapath)
        stnsesinfo.set_obsfolderinfo(LOFARdatTYPE, bsxSTobsEpoch,
                                     modeparms.rcumode2band(rcumode), integration,
                                     duration_scan)
        stnsesinfo.write_session_header(datapath)
        data_url = "{}:{}".format(self.get_stnid(), datapath)
        return data_url

    #####################
    # BEGIN: TBB services

    def _setupTBBs(self):
        """Setup transient buffer boards and start recording."""
        print("In setupTBBs: Freeing TBBs")
        self.lcu_interface.run_tbbctl(free=True)
        print("In setupTBBs: Setting TBB transient mode on rspctl")
        self.lcu_interface.run_rspctl(tbbmode='transient')
        time.sleep(1)
        print("In setupTBBs: Allocating TBBs")
        self.lcu_interface.run_tbbctl(alloc=True)
        time.sleep(1)
        print("In setupTBBs: Setting TBB transient mode on tbbctl")
        self.lcu_interface.run_tbbctl(mode='transient')
        time.sleep(1)
        print("In setupTBBs: Start TBB recording")
        self.lcu_interface.run_tbbctl(record=True)
        print("In setupTBBs: Finished setting up TBBs & started recording")

    def _freezeTBBdata(self):
        """Stop TBB recording."""
        print("In freezeTBBdata: Stopping TBB recording")
        self.lcu_interface.run_tbbctl(stop=True)
        print("In freezeTBBdata: Stopping any dummy beam")
        self.lcu_interface.stop_beam()

    def _startTBBdataStream(self, duration_scan):
        """Stream duration_scan seconds of TBB data out of the LCU to
        datataking node."""
        # Set delay between subsequent frames. One delay unit is 5us.
        udpdelay = 500  # (Previously 100)
        Nqfreq = 100e6
        nrpages = str(int(duration_scan*2*Nqfreq/1024))  # One page is 1024 samples.
                                                    # Normal sampling frequency
                                                    # is 200MHz.
        self.lcu_interface.run_tbbctl(select='0:15,16:31,32:47', storage='lofarA1')
        self.lcu_interface.run_tbbctl(select='48:63,64:79,80:95', storage='lofarA2')
        self.lcu_interface.run_tbbctl(select='96:111,112:127,128:143', storage='lofarA3')
        self.lcu_interface.run_tbbctl(select='144:159,160:175,176:191', storage='lofarA4')

        self.lcu_interface.run_tbbctl(cepdelay=str(udpdelay))

        self.lcu_interface.run_tbbctl(select='0:15', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='48:63', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='96:111', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='144:159', readall=nrpages, backgroundJOB='locally')

        self.lcu_interface.run_tbbctl(select='16:31', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='64:79', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='112:127', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='160:175', readall=nrpages, backgroundJOB='locally')

        self.lcu_interface.run_tbbctl(select='32:47', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='80:95', readall=nrpages, backgroundJOB='locally')
        self.lcu_interface.run_tbbctl(select='128:143', readall=nrpages, backgroundJOB='locally')
        # Last one is not put in background so the parent process blocks
        # until finished.
        self.lcu_interface.run_tbbctl(select='176:191', readall=nrpages) # backgroundJOB='locally')

    def do_tbb(self, duration_scan, band, start_after=0):
        """Record duration_scan seconds of TBB data from rcumode."""

        observer = self.stnsesinfo.projectmeta['observer']
        project = self.stnsesinfo.projectmeta['projectname']
        observationID = "Null"

        # Start a beam
        pointing = modeparms.stdPointings('Z')
        freqband = modeparms.FrequencyBand(band)
        # FrequencyBand obtained from band spec sets 8 bit mode,
        # so create a new FrequencyBand object with only center frequency
        freqlo, freqhi = freqband.edgefreqs()
        freq0 = (freqlo+freqhi)/2.0
        actualfb = modeparms.FrequencyBand(freq0)
        antset = actualfb.antsets[0]
        self.streambeams(actualfb, pointing)

        print("Setting up TBBs")
        self._setupTBBs()

        print("Will start TBB recording in {}s".format(start_after))
        time.sleep(start_after)

        print("Will freeze TBB recording in {}s".format(duration_scan))
        time.sleep(duration_scan)  # Arbitrary time to trigger
        print("Sending trigger to TBBs")
        self._freezeTBBdata()

        # Start data capture process locally
        dalcap = \
            multiprocessing.Process(target=capture_data_DAL1,
                                    args=(self.tbbraw2h5cmd, self.tbbh5dumpdir,
                                          observer, antset, project,
                                          observationID,False))
        dalcap.start()

        print("Streaming {}s of TBB data out of LCU".format(duration_scan))
        self._startTBBdataStream(float(duration_scan))
        dalcap.join()

    # END: TBB services
    ###################

    def executeblock(self, scans):
        for scan in scans:
            prg = programs.BasicObsPrograms(self)

            if scan['obsprog'] != 'None':
                obsfun, obsargs_sig = prg.getprogram(scan['obsprog'])

                # Prepare observation arguments:
                freqbndoobj = modeparms.FrequencyBand(scan['beam']['freqspec'])
                pointsrc = scan['beam']['pointing']
                try:
                    pointing = modeparms.stdPointings(pointsrc)
                except KeyError:
                    try:
                        phi, theta, ref = pointsrc.split(',', 3)
                        # FIXME:  (not always going to be correct)
                        pointing = pointsrc
                    except ValueError:
                        raise ValueError(
                            "Error: %s invalid pointing syntax".format(pointsrc))
                if scan['rec_stat'] is not None:
                    integration = scan['rec_stat']['integration']
                else:
                    integration = 1
                duration_tot = scan['duration_tot']
                if integration > duration_tot:
                    raise (ValueError, "integration {} is longer than duration_scan {}."
                           .format(integration, duration_tot))
                obsargs_in = {'starttime': scan['starttime'],
                              'freqbndobj': freqbndoobj,
                              'pointsrc': pointsrc,
                              'pointing': pointing,
                              'duration_tot': duration_tot,
                              'integration': integration}
                obsargs = {k: obsargs_in[k] for k in obsargs_sig}
                self.do_obsprog(scan['starttime'], obsfun, obsargs)
            else:
                pass

    def do_obsprog(self, starttime, obsfun, obsargs):
        """At starttime execute the observation program specified by the obsfun method
        pointer and run with arguments specified by obsargs dict.
        """

        # Setup Calibration tables on LCU:
        CALTABLESRC = 'default'   # FIXME put this in args
        ## (Only BST uses calibration tables)
        # Choose between 'default' or 'local'
        self.lcu_interface.selectCalTable(CALTABLESRC)

        # Prepare for obs program.
        try:
            self.goto_observingstate()
        except RuntimeError as e:
            raise RuntimeError(e)

        # Run the observation program:
        obsinfolist = obsfun(**obsargs)
        # Stop program beam
        self.lcu_interface.stop_beam()

        if obsinfolist is not None:
            # Get some metadata about operational settings:
            # e.g. caltables used
            caltabinfos = []
            freqbndobj = obsargs['freqbndobj']
            for rcumode in freqbndobj.rcumodes:
                caltabinfo = self.lcu_interface.getCalTableInfo(rcumode)
                caltabinfos.append(caltabinfo)
            for i in range(len(obsinfolist)):
                obsinfolist[i].caltabinfos = caltabinfos
            obsinfo = copy.copy(obsinfolist[0])
            obsinfo.sb = freqbndobj.sb_range[0]

            # Move data to archive
            bsxSTobsEpoch, datapath = obsinfo.getobsdatapath(self.LOFARdataArchive)

            self.movefromlcu(self.lcu_interface.lcuDumpDir + "/*.dat", datapath)
            for curr_obsinfo in obsinfolist:
                curr_obsinfo.create_LOFARst_header(datapath)
            # Prepare metadata for session on this station.
            stnsesinfo = copy.deepcopy(self.stnsesinfo)
            stnsesinfo.set_obsfolderinfo(obsinfo.LOFARdatTYPE, bsxSTobsEpoch,
                                         freqbndobj.arg, obsinfo.integration,
                                         obsinfo.duration_scan, obsinfo.pointing)
            stnsesinfo.write_session_header(datapath)
            data_url = "{}:{}".format(self.get_stnid(), datapath)
        return None

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
