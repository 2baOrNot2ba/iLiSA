#!/usr/bin/python
"""Package that provides functions for setting up and running observations
via station controller objects. This package knows about the data archive and
should not run anything directly on LCU."""


import time
import datetime
import subprocess
import os
import sys
import multiprocessing
import copy

import ilisa.observations.lcuinterface as stationcontrol
import ilisa.observations.modeparms as modeparms
import ilisa.observations.dataIO as dataIO
import ilisa.observations.beamformedstreams.bfbackend as bfbackend


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
        if not os.path.exists(dest):
            os.makedirs(dest)
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
        rcu_setup_cmd = self.lcu_interface.rcusetup(bits, attenuation)
        return rcu_setup_cmd, beamctl_cmds

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
        # From swlevel 0 it takes about 1:30min? to reach swlevel 3
        print("Booting @ {}".format(datetime.datetime.utcnow()))
        return st

    def _setup_tof(self, elemsOn=modeparms.elOn_gILT):
        """Setup (HBA) tiling off mode."""
        # NOTE: LCU must be in swlevel=2 to run SEPTON!
        self.lcu_interface.set_swlevel(2)
        # self.stationcontroller.turnoffElinTile_byTile(elemsOn) # Alternative
        self.lcu_interface.turnoffElinTile_byEl(elemsOn)
        septonconf = modeparms.elementMap2str(elemsOn)
        return septonconf

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

    def main_scan(self, freqbndobj, integration, duration_tot, pointing, pointsrc,
                  starttime='NOW', rec_stat_type=None, rec_bfs=False, duration_scan=None,
                  do_acc=False, allsky=False, warmup=False):
        """Run a generic scan.

        Parameters
        ----------
            freqbndobj: FrequencyBand
                Frequency specification of scan.
            integration: float
                Integration time in seconds.
            duration_tot: float
                Total duration of scan in seconds.
            pointing:
                Pointing direction of scan.
            pointsrc: str
                The field name.
            starttime: str
                The time at which scan should start.
            rec_stat_type: str
                The type of LOFAR statistic to be recorded. Can be either: 'bst', 'sst',
                or 'xst'.
            rec_bfs: bool
                Record the beam-formed stream.
            duration_scan: float
                The duration of one frequency sweep within the scan.
            do_acc: bool
                Perform calibration observation mode on station. Also known as ACC
                mode. The actual total duration will at most be duration_tot_req.
                (Usually it will be shorter so it fits within the cadence of whole ACC
                aquisition, which is 512+7=519 seconds).

                Some technical details: swlevel needs to cycle down to 2 (or less) and then to 3.
                If swlevel is kept at 3 (i.e. exit_obsstate=False), then ACC will continue to be
                produced, until swlevel goes below 2.

                ACC files are autocovariance-cubes: the covariance of all array
                elements with each as a function of subband. These files are generated
                by the MAC service called CalServer. It run at swlevel 3 and is
                configured in the file lofar/etc/CalServer.conf. Note subband
                integration is always 1s, so ACC file is dumped after 512 seconds.
            allsky: bool, optional
                HBA allsky mode.
            warmup: bool, optional
                Do a beam warmup prior to doing scan.
        """
        req_allsky = allsky; del allsky
        arraytype = freqbndobj.antsets[-1][:3]
        # Mode logic
        todo_tof = False  # This is the case when arraytype=='LBA'
        if req_allsky and arraytype == 'HBA':
            if pointing is not None:
                raise ValueError('hba-allsky and beam cannot run simultaneously.')
            if rec_stat_type == 'bst':
                raise ValueError('bst and hba-allsky cannot be combined.')
            if do_acc:
                raise ValueError('acc and hba-allsky cannot be combined.')
            if rec_bfs:
                raise ValueError('bfs and hba-allsky cannot be combined.')
            todo_tof = True
        if not req_allsky and pointing is None:
            pointing = 'Z'

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

        duration_tot_req = duration_tot
        band = freqbndobj.rcubands[0]
        rcumode = freqbndobj.rcumodes[0]

        if todo_tof:
            septonconf = self._setup_tof()
        else:
            septonconf = None

        if do_acc:
            # Also duration of ACC sweep since each sb is 1 second.
            dur1acc = modeparms.TotNrOfsb  # Duration of one ACC
            interv2accs = 7  # time between end of one ACC and start of next one
            acc_cadence = dur1acc+interv2accs  # =519s time between start of two ACCs
            (nraccs, timrest) = divmod(duration_tot_req, acc_cadence)
            if timrest > dur1acc:
                nraccs += 1
            duration_tot = nraccs*acc_cadence-interv2accs
            if duration_tot != duration_tot_req:
                print("""Note: will use total duration {}s to fit with ACC
                      cadence.""".format(duration_tot))
            #sst_integration = int(acc_cadence)

            # Make sure swlevel=<2
            self.lcu_interface.set_swlevel(2)

            # Set CalServ.conf to dump ACCs:
            self.lcu_interface.acc_mode(enable=True)

            # Boot to swlevel 3 so the calserver service starts
            self.lcu_interface.set_swlevel(3)

        # Wait until it is time to start
        pause = 5  # Sufficient?
        if starttime != "NOW":
            starttimestr = starttime.strftime("%Y-%m-%dT%H:%M:%S")
        else:
            starttimestr = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
        st = self._waittoboot(starttimestr, pause)

        # Necessary since fork creates multiple instances of myobs and each one
        # will call it's __del__ on completion and __del__ shutdown...
        shutdown = self.halt_observingstate_when_finished
        self.halt_observingstate_when_finished = False
        self.exit_check = False

        caltabinfos = []
        if pointing is not None:
            # Get metadata about caltables to be used
            for rcumode in freqbndobj.rcumodes:
                caltabinfo = self.lcu_interface.getCalTableInfo(rcumode)
                caltabinfos.append(caltabinfo)
            if warmup:
                # Dummy or hot beam start: (takes about 10sec)
                # TODO: This seems necessary, otherwise beamctl will not start up next time,
                #       although it should not have to necessary.)
                print("Running warmup beam... @ {}".format(datetime.datetime.utcnow()))
                self.streambeams(freqbndobj, pointing)
                self.lcu_interface.stop_beam()

            print("Pause {}s after boot.".format(pause))
            time.sleep(pause)

            # Real beam start:
            print("Now running real beam... @ {}".format(datetime.datetime.utcnow()))
            rcu_setup_cmd, beamctl_cmds = self.streambeams(freqbndobj, pointing)
            beamstarted = datetime.datetime.utcnow()
            timeleft = st - beamstarted
            if timeleft.total_seconds() < 0.:
                starttimestr = beamstarted.strftime("%Y-%m-%dT%H:%M:%S")
            print("(Beam started) Time left before recording: {}".format(
                timeleft.total_seconds()))
        else:
            print("No pointing...")
            rcu_setup_cmd = ""
            beamctl_cmds = ""

        if rec_bfs:
            bf_data_dir = os.path.join(self.bf_data_dir,'proj{}'.format(
                self.stnsesinfo.projectmeta['projectID']),'')
            port0 = self.bf_port0
            stnid = self.lcu_interface.stnid
            lanes = tuple(freqbndobj.getlanes().keys())
            bfbackend.rec_bf_streams(starttimestr, duration_tot, lanes, band, bf_data_dir,
                                     port0, stnid)
            bfspaths = []
            for lane in lanes:
                outdumpdir, outarg, datafileguess, dumplogname =\
                    bfbackend.bfsfilepaths(lane, starttimestr, band, bf_data_dir, port0,
                                           stnid)
                bfspaths.append(datafileguess)
        else:
            print("Not recording")
        sys.stdout.flush()

        if rec_stat_type is not None:
            # Record statistic for duration_tot seconds
            obsinfolist = []
            if rec_stat_type == 'bst':
                rspctl_cmd = self.lcu_interface.rec_bst(integration, duration_tot)
                # beamlet statistics also generate empty *01[XY].dat so remove:
                self.lcu_interface.rm(self.lcu_interface.lcuDumpDir + "/*01[XY].dat")
                obsdatetime_stamp = self.get_data_timestamp(-1)
                curr_obsinfo = dataIO.ObsInfo()
                curr_obsinfo.setobsinfo_fromparams('bst', obsdatetime_stamp, beamctl_cmds,
                                                   rspctl_cmd, caltabinfos)
                obsinfolist.append(curr_obsinfo)
            elif rec_stat_type == 'sst':
                caltabinfo = ""
                rspctl_cmd = self.lcu_interface.rec_sst(integration, duration_tot)
                obsdatetime_stamp = self.get_data_timestamp(-1)
                curr_obsinfo = dataIO.ObsInfo()
                curr_obsinfo.setobsinfo_fromparams('sst', obsdatetime_stamp, beamctl_cmds,
                                                   rspctl_cmd, caltabinfo)
                obsinfolist.append(curr_obsinfo)
            elif rec_stat_type == 'xst':
                caltabinfo = ""  # No need for caltab info for xst data
                nrsubbands = freqbndobj.nrsubbands()
                if duration_scan is None:
                    if nrsubbands > 1:
                        duration_scan = integration
                    else:
                        duration_scan = duration_tot
                # FIXME When duration_tot is too small for 1 rep this will fail badly.
                (rep, rst) = divmod(duration_tot, duration_scan * nrsubbands)
                rep = int(rep)
                # Repeat rep times (the freq sweep)
                for itr in range(rep):
                    # Start freq sweep
                    for sb_rcumode in freqbndobj.sb_range:
                        if ':' in sb_rcumode:
                            sblo, sbhi = sb_rcumode.split(':')
                            subbands = range(int(sblo), int(sbhi) + 1)
                        else:
                            subbands = [int(sb) for sb in sb_rcumode.split(',')]
                        for subband in subbands:
                            # Record data
                            rspctl_cmd = self.lcu_interface.rec_xst(subband, integration,
                                                                    duration_scan)
                            obsdatetime_stamp = self.get_data_timestamp(-1)
                            curr_obsinfo = dataIO.ObsInfo()
                            curr_obsinfo.setobsinfo_fromparams('xst', obsdatetime_stamp,
                                                               beamctl_cmds, rspctl_cmd,
                                                               caltabinfos=caltabinfo,
                                                               septonconf=septonconf)
                            obsinfolist.append(curr_obsinfo)
            else:
                raise Exception('LOFAR statistic datatype "{}" unknown.\
                                (Known are bst, sst, xst)'.format(rec_stat_type))
        else:
            obsinfolist = None

        if not do_acc and not rec_bfs and rec_stat_type is None:
            # Since we're not recording anything, just do nothing for the duration_tot.
            time.sleep(duration_tot)

        # Finished recording
        self.lcu_interface.stop_beam()

        # Determine path where this scan should be stored
        projpath = os.path.join(self.LOFARdataArchive,
                                "proj{}".format(self.stnsesinfo.projectmeta['projectID']))
        # - Make a scan ID. (Currently based on starttime of scan)
        try:
            scanid = self.get_data_timestamp(-1)
        except Exception:
            scanid = beamstarted.strftime("%Y%d%m_%H%M%S")
        scanpath = os.path.join(projpath, "scan_{}".format(scanid))

        if todo_tof:
            self.lcu_interface.set_swlevel(3)

        acc_url = None
        if do_acc:
            # Switch back to normal state i.e. turn-off ACC dumping:
            self.lcu_interface.set_swlevel(2)
            self.lcu_interface.acc_mode(enable=False)
            self.lcu_interface.set_swlevel(3)

            # Transfer data from LCU to DAU
            obsdatetime_stamp = self.get_data_timestamp(ACC=True)
            accsrcfiles = self.lcu_interface.ACCsrcDir + "/*.dat"
            acc_destfolder = \
                os.path.join(scanpath,
                             '{}_{}_rcu{}_dur{}'.format(self.lcu_interface.stnid,
                                                        obsdatetime_stamp, rcumode,
                                                        duration_tot))
            if int(rcumode) > 3:
                acc_destfolder += "_"+pointsrc
            acc_destfolder += "_acc"
            if os.path.exists(acc_destfolder):
                print("Appropriate directory exists already (will put data here)")
            else:
                print("Creating directory "+acc_destfolder+" for ACC "+str(duration_tot)
                      + " s rcumode="+str(rcumode)+" calibration")
                os.makedirs(acc_destfolder)

            # Move ACC dumps to storage
            self.movefromlcu(accsrcfiles, acc_destfolder)
            accdestfiles = os.listdir(acc_destfolder)

            # - Create project header

            sesinfo_acc = copy.deepcopy(self.stnsesinfo)
            acc_integration = 1.0
            sesinfo_acc.set_obsfolderinfo('acc', obsdatetime_stamp, band, acc_integration,
                                          duration_tot, pointing)
            sesinfo_acc.write_scan_rec(acc_destfolder)

            # - Create header for each ACC file
            for destfile in accdestfiles:
                filedatestr, filetimestr, _ = destfile.split('_', 2)
                filedtstr = '_'.join([filedatestr, filetimestr])
                sesinfo_acc.new_obsinfo()
                sesinfo_acc.obsinfos[-1].setobsinfo_fromparams('acc', filedtstr,
                                                               beamctl_cmds, '')
                sesinfo_acc.obsinfos[-1].create_LOFARst_header(acc_destfolder)

            acc_url = "{}:{}".format(self.get_stnid(), acc_destfolder)

        stat_url = None
        if rec_stat_type is not None and obsinfolist is not None:
            obsinfo = copy.copy(obsinfolist[0])
            obsinfo.sb = freqbndobj.sb_range[0]

            # Move data to archive
            bsxSTobsEpoch, datapath = obsinfo.getobsdatapath(scanpath)
            self.movefromlcu(self.lcu_interface.lcuDumpDir + "/*.dat", datapath)
            for curr_obsinfo in obsinfolist:
                curr_obsinfo.create_LOFARst_header(datapath)
            # Prepare metadata for session on this station.
            stnsesinfo = copy.deepcopy(self.stnsesinfo)
            stnsesinfo.set_obsfolderinfo(obsinfo.LOFARdatTYPE, bsxSTobsEpoch,
                                         freqbndobj.arg, obsinfo.integration,
                                         obsinfo.duration_scan, obsinfo.pointing)
            stnsesinfo.write_scan_rec(datapath)
            stat_url = "{}:{}".format(self.get_stnid(), datapath)

        bfs_url = None
        if rec_bfs:
            headertime = datetime.datetime.strptime(starttimestr, "%Y-%m-%dT%H:%M:%S"
                                                    ).strftime("%Y%m%d_%H%M%S")
            stnsesinfo_bfs = copy.deepcopy(self.stnsesinfo)
            stnsesinfo_bfs.new_obsinfo()
            stnsesinfo_bfs.obsinfos[-1].setobsinfo_fromparams('bfs', headertime,
                                                              beamctl_cmds, rcu_setup_cmd,
                                                              caltabinfos)
            bsxSTobsEpoch, datapath = stnsesinfo_bfs.obsinfos[-1].getobsdatapath(scanpath)
            print("Creating BFS destination folder on DPU:\n{}".format(datapath))
            os.makedirs(datapath)
            # Make soft links to actual BFS files
            for lane in lanes:
                os.symlink(bfspaths[lane], os.path.join(datapath,
                                                        os.path.basename(bfspaths[lane])))
            # TODO: Move dump log files
            stnsesinfo_bfs.obsinfos[-1].create_LOFARst_header(datapath)
            integration = None
            stnsesinfo_bfs.set_obsfolderinfo('bfs', headertime, band, integration,
                                             duration_tot, pointing)
            stnsesinfo_bfs.write_scan_rec(datapath)
            bfs_url = "{}:{}".format(self.get_stnid(), datapath)

        self.lcu_interface.cleanup()
        # Necessary due to possible forking
        self.halt_observingstate_when_finished = shutdown

        return bfs_url, stat_url, acc_url

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
            stnsesinfo.write_scan_rec(datapath)
            data_url = "{}:{}".format(self.get_stnid(), datapath)
        return None


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
