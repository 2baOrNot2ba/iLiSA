import sys
import os
import time
import datetime
import copy
import inspect
import ilisa.observations.modeparms as modeparms
import ilisa.observations.beamformedstreams.bfbackend as bfbackend

class BasicObsPrograms(object):

    def __init__(self, stationdriver):
        self.stationdriver = stationdriver
        self.lcu_interface = stationdriver.lcu_interface

    def getprogram(self, programname):
        programpointer = getattr(self, programname)
        defargstart = None
        # # FIXME Only return args without defaults
        # defarg = inspect.getargspec(programpointer).defaults
        # if defarg is not None:
        #     defargstart = -len(defarg)
        programargs = inspect.getargspec(programpointer).args[1:defargstart]
        return programpointer, programargs

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

    def do_bfs_OW(self, freqbndobj, duration_tot, pointsrc, bfdsesdumpdir,
                  starttime='NOW'):
        """Record BeamFormed Streams (BFS) with particular beamlet allocation."""

        band = freqbndobj.rcubands[0]
        ###
        bits = 8  # 8
        attenuation = None
        # Subbands allocation
        if band == '10_90' or band == '30_90':
            # LBA
            lanes = (0, 1)  # (0,1)
            beamletIDs = '0:243'  # '0:243'
            subbandNrs = '164:407'  # '164:407'
        elif band == '110_190':
            # HBAlo
            lanes = (0, 1, 2, 3)  # Normally (0,1,2,3) for all 4 lanes.
            beamletIDs = '0:487'
            subbandNrs = '12:499'
        elif band == '210_250':
            # HBAhi
            lanes = (0, 1)
            beamletIDs = '0:243'
            subbandNrs = '12:255'
        else:
            raise ValueError(
                "Wrong band: should be 10_90 (LBA), 110_190 (HBAlo) or 210_250 (HBAhi).")
        pointing = modeparms.normalizebeamctldir(pointsrc)

        # Wait until it is time to start
        pause = 5  # Sufficient?
        if starttime != "NOW":
            starttimestr = starttime.strftime("%Y-%m-%dT%H:%M:%S")
        else:
            starttimestr = datetime.datetime.utcnow().strftime("%Y-%m-%dT%H:%M:%S")
        st = self.stationdriver._waittoboot(starttimestr, pause)

        # From swlevel 0 it takes about 1:30min? to reach swlevel 3
        print("Booting @ {}".format(datetime.datetime.utcnow()))

        # Necessary since fork creates multiple instances of myobs and each one
        # will call it's __del__ on completion and __del__ shutdown...
        shutdown = self.stationdriver.halt_observingstate_when_finished
        self.stationdriver.halt_observingstate_when_finished = False
        self.stationdriver.exit_check = False

        # BEGIN Dummy or hot beam start: (takes about 10sec)
        # TODO: This seems necessary, otherwise beamctl will not start up next time,
        #       although it should not have to necessary.)
        print("Running warmup beam... @ {}".format(datetime.datetime.utcnow()))
        self.lcu_interface.run_beamctl(beamletIDs, subbandNrs, band, pointing)
        self.lcu_interface.rcusetup(bits,
                                    attenuation)  # setting bits also seems necessary
        self.lcu_interface.stop_beam()
        # END Dummy or hot start

        print("Pause {}s after boot.".format(pause))
        time.sleep(pause)

        # Real beam start:
        print("Now running real beam... @ {}".format(datetime.datetime.utcnow()))
        beamctl_cmds = self.lcu_interface.run_beamctl(beamletIDs, subbandNrs, band,
                                                     pointing)
        rcu_setup_cmd = self.lcu_interface.rcusetup(bits, attenuation)
        beamstart = datetime.datetime.utcnow()
        timeleft = st - beamstart
        if timeleft.total_seconds() < 0.:
            starttimestr = beamstart.strftime("%Y-%m-%dT%H:%M:%S")
        print("(Beam started) Time left before recording: {}".format(
            timeleft.total_seconds()))

        REC = True
        if REC == True:
            port0 = self.stationdriver.bf_port0
            stnid = self.lcu_interface.stnid
            bfbackend.rec_bf_streams(starttimestr, duration_tot, lanes, band,
                                     bfdsesdumpdir, port0, stnid)
        else:
            print("Not recording")
            time.sleep(duration_tot)
        sys.stdout.flush()
        self.lcu_interface.stop_beam()
        self.stationdriver.halt_observingstate_when_finished = shutdown
        return None
