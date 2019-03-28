import ilisa.observations.modeparms as modeparms
import ilisa.observations.dataIO as dataIO

class BasicObsPrograms(object):

    def __init__(self, stationdriver):
        self.stationdriver = stationdriver
        self.lcu_interface = stationdriver.lcu_interface

    def getprogram(self, programname):
        return getattr(self, programname)

    def _fix_obsargs(self, kwargs_in):
        """Check and transform observational arguments to useful arguments."""
        kwargs_out = {}
        try:
            kwargs_out['pointing'] = modeparms.stdPointings(kwargs_in['pointsrc'])

        except KeyError:
            try:
                phi, theta, ref = kwargs_in['pointsrc'].split(',', 3)
                # FIXME Not always going to be correct
            except ValueError:
                raise ValueError("Error: %s invalid pointing syntax"\
                                 .format(kwargs_in['pointsrc']))
            else:
                kwargs_out['pointing'] = kwargs_in['pointsrc']
        if kwargs_in['integration'] > kwargs_in['duration_tot']:
            raise (ValueError, "integration {} is longer than duration_scan {}."
                               .format(kwargs_in['integration'],
                                       kwargs_in['duration_tot']))
        else:
            kwargs_out['integration'] = kwargs_in['integration']
            kwargs_out['duration_tot'] = kwargs_in['duration_tot']
        kwargs_out['freqbndobj'] = modeparms.FrequencyBand(kwargs_in['freqbndarg'])

        return kwargs_out

    def _streambeams_mltfreq(self, freqbndobj, pointing, recDuration=float('inf'),
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

    def do_bstnew(self, freqbndobj, integration, duration_tot, pointing):
        """Do a Beamlet STatistic (bst) recording on station. frequency in Hz.
        """
        duration_scan = duration_tot
        obsinfolist = []

        # Start beamforming
        (rcu_setup_cmd, beamctl_cmds) = self._streambeams_mltfreq(freqbndobj, pointing)

        # Record data
        rspctl_cmd = self.lcu_interface.rec_bst(integration, duration_scan)

        # beamlet statistics also generate empty *01[XY].dat so remove:
        self.lcu_interface.rm(self.lcu_interface.lcuDumpDir+"/*01[XY].dat")

        # Collect some obsinfo (calinfo will be added later):
        obsdatetime_stamp = self.stationdriver.get_data_timestamp()
        curr_obsinfo = dataIO.ObsInfo()
        curr_obsinfo.setobsinfo_fromparams('bst', obsdatetime_stamp, beamctl_cmds,
                                           rspctl_cmd)
        obsinfolist.append(curr_obsinfo)
        return obsinfolist

    def do_sstnew(self, freqbndobj, integration, duration_tot, pointing):
        """Run an sst static."""

        duration_scan =  duration_tot
        obsinfolist = []

        # Use this to do a system temperature measurement.
        SYS_TEMP_MEAS = False

        # Start beamforming
        rcu_setup_CMD, beamctl_cmds = \
            self._streambeams_mltfreq(freqbndobj, pointing)
        if SYS_TEMP_MEAS:
            lbbalnaoff_CMD = self.lcu_interface.turnoffLBA_LNAs()
            #beamctl_CMD += "\n"+lbbalnaoff_CMD

        # Get some metadata about operational settings:
        caltabinfo = ""    # No need for caltab info

        # Record data
        rspctl_cmd = self.lcu_interface.rec_sst(integration, duration_scan)
        obsdatetime_stamp = self.stationdriver.get_data_timestamp(-1)
        curr_obsinfo = dataIO.ObsInfo()
        curr_obsinfo.setobsinfo_fromparams('sst', obsdatetime_stamp,
                                           beamctl_cmds, rspctl_cmd,
                                           caltabinfo)
        obsinfolist.append(curr_obsinfo)
        return obsinfolist

    def do_xstnew(self, freqbndobj, integration, duration_tot, pointing,
                  duration_scan=None):
        """New: Run an xst statistic towards the given pointing. This corresponds to
        a crosscorrelation of all elements at the given frequency for
        integration seconds over a duration_scan of seconds."""

        caltabinfo = ""  # No need for caltab info for xst data
        obsinfolist = []
        nrsubbands = freqbndobj.nrsubbands()
        if duration_scan is None:
            if nrsubbands > 1:
                duration_scan = integration
            else:
                duration_scan = duration_tot
        # FIXME When duration_tot is too small for 1 rep this will fail badly.
        (rep, rst) = divmod(duration_tot, duration_scan*nrsubbands)
        rep = int(rep)
        # Start beamforming
        rcu_setup_CMD, beamctl_cmds = self._streambeams_mltfreq(freqbndobj, pointing)
        # Repeat rep times (the freq sweep)
        for itr in range(rep):
            # Start freq sweep
            for sb_rcumode in freqbndobj.sb_range:
                if ':' in sb_rcumode:
                    sblo, sbhi = sb_rcumode.split(':')
                    subbands = range(int(sblo),int(sbhi)+1)
                else:
                    subbands = [int(sb) for sb in sb_rcumode.split(',')]
                for subband in subbands:
                    # Record data
                    rspctl_cmd = self.lcu_interface.rec_xst(subband, integration,
                                                            duration_scan)
                    obsdatetime_stamp = self.stationdriver.get_data_timestamp(-1)
                    curr_obsinfo = dataIO.ObsInfo()
                    curr_obsinfo.setobsinfo_fromparams('xst', obsdatetime_stamp,
                                                       beamctl_cmds, rspctl_cmd,
                                                       caltabinfo)
                    obsinfolist.append(curr_obsinfo)
        return obsinfolist
