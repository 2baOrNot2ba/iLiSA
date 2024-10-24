"""This package is for the parameters involved in observation modes.
"""
import os
import argparse
import math
import numpy
import datetime
import warnings

from ilisa.operations import DATETIMESTRFMT

ANTENNA_SETS = ['LBA_INNER', 'LBA_OUTER',
                'LBA_SPARSE_EVEN', 'LBA_SPARSE_ODD',
                'LBA_X', 'LBA_Y',
                'HBA_DUAL', 'HBA_JOINED', 'HBA_ZERO', 'HBA_ONE']
NQFREQ_NOM = 100.0e6  # Nominal Nyquist frequency in Hz
TotNrOfsb = 512  # Total number of subbands. (Subbands numbered 0:511)
nrofrcus = 192  # Number of RCUs
MIN_STATS_INTG = 1.0  # Minimum integration for statistics data in seconds.
MAX_NRLANES = 4
BASE_NR_BEAMLETS = 244  # Base number of total beamlets available
NRBEAMLETSBYBITS = {16:   BASE_NR_BEAMLETS,
                    8:  2*BASE_NR_BEAMLETS,
                    4:  4*BASE_NR_BEAMLETS}
ACC_DUR = TotNrOfsb*MIN_STATS_INTG  # Duration of one ACC in s
INTERV2ACCS = 7.0  # Interval between 2 ACCs in s
# SEPTON configurations:
#        1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24

elOn_step = \
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15, 0, 1, 2, 3, 4, 5, 6, 7,
     8, 9,10,11,12,13,14,15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15,
     0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15, 0, 1, 2, 3, 4, 5, 6, 7,
     8, 9,10,11,12,13,14,15, 0, 1, 2, 3, 4, 5, 6, 7, 8, 9,10,11,12,13,14,15]
elOn_gILT = \
    [15, 0,15, 3, 9,15,14, 2, 0, 3, 4,14,10, 8, 5,15,12, 0, 2,11, 3,12,12, 1,
      5, 4, 4, 8, 6, 3, 0, 5, 3,11, 3, 2, 8,15,13, 8, 3, 2, 9, 1,14, 8, 8, 0,
     12,13, 0,11,15, 3,12, 3,13, 3,10, 5, 0,10, 1, 6, 4,10, 3,15, 3,14, 0,12,
      0, 7, 0,12, 7, 3,13, 0, 7, 3,15, 4,14, 4, 3, 8, 4, 9,12, 0,14, 9, 3,11]
elOn_Generic_Int_201512 = \
    [ 0, 5, 3, 1, 8, 3,12,15,10,13,11, 5,12,12, 5, 2,10, 8, 0, 3, 5, 1, 4, 0,
     11, 6, 2, 4, 9,14,15, 3, 7, 5,13,15, 5, 6, 5,12,15, 7, 1, 1,14, 9, 4, 9,
      3, 9, 3,13, 7,14, 7,14, 2, 8, 8, 0, 1, 4, 2, 2,12,15, 5, 7, 6,10,12, 3,
      3,12, 7, 4, 6, 0, 5, 9, 1,10,10,11, 5,11, 7, 9, 7, 6, 4, 4,15, 4, 1,15]
# hex for elOn_Generic_Int_201512
# f0f39fe2034ea85fc02b3cc1544863053b328fd83291e880cd0bf3c3d3a50a164a3f3e0c070c73d073f4e43849c0e93b
elOn_same_el = 0
# elOn_same = [elOn_same_el for elemNr in range(stationcontrol.nrTiles)]
elemsOn = elOn_Generic_Int_201512  # elOn_same or elOn_step or elOn_gILT or ...


def beamctl_args2cmds(beamlets, subbands, band, anadigdir, rcus='all',
                      beamdurstr='', antset=''):
    """
    Create a beamctl command string from the given arguments

    Parameters
    ----------
    beamlets: str
        Beamlets CLI argument to beamctl.
    subbands: str
        Subbands CLI argument to beamctl.
    band: str
        Band CLI argument to beamctl.
    anadigdir: str
        Analogue or digital direction CLI argument to beamctl.
    rcus: str
        RCU CLI argument to beamctl.
    beamdurstr: str
        Beam duration CLI argument to beamctl.
    antset: str
        Antset CLI argument to beamctl.
    """
    if rcus == 'all':
        rcus = '0:191'
    if beamdurstr != '':
        beamdurstr = ',' + beamdurstr
    anadir = anadigdir
    digdir = anadigdir

    try:
        # See if band is actually old rcumode 3,5,7 etc
        band = rcumode2band(band)
    except ValueError:
        pass  # It's not an rcumode. Assume it's a proper band descriptor
    if not antset:
        antset = band2antset_eu(band)
    beamctl_cmd = ("beamctl --antennaset=" + antset + " --rcus=" + rcus
                   + " --band=" + band + " --beamlets=" + beamlets
                   + " --subbands=" + subbands
                   + " --anadir=" + anadir + beamdurstr
                   + " --digdir=" + digdir + beamdurstr)
    return beamctl_cmd


def parse_beamctl_args(beamctl_str):
    """
    Parse beamctl command arguments
    """
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

def rcusetup_args2cmds(bits, attenuation, mode=None):
    """
    Convert RCU arguments to command lines
    """
    rcusetup_cmds = []
    rcusetup_cmd = "rspctl --bitmode=" + str(bits)
    rcusetup_cmds.append(rcusetup_cmd)
    if attenuation:
        # NOTE attenuation only set when beamctl is runnning.
        rcusetup_cmd = "rspctl --rcuattenuation=" + str(attenuation)
        rcusetup_cmds.append(rcusetup_cmd)
    if mode:
        rcusetup_cmd = "rspctl --mode=" + str(mode)
        rcusetup_cmds.append(rcusetup_cmd)
    return rcusetup_cmds

def rspctl_stats_args2cmds(bsxtype, integration, duration, subband=0):
    """\
    Run rspctl statistics command
    """
    rspctl_cmds = []
    if bsxtype == 'xst':
        rspctl_cmd = "rspctl --xcsubband="+str(subband)
        rspctl_cmds.append(rspctl_cmd)
    if bsxtype == 'bst':
        stats_flag_and_val = '--statistics=beamlet'
    elif bsxtype == 'sst':
        stats_flag_and_val = '--statistics=subband'
    else:
        stats_flag_and_val = '--xcstatistics'
    rspctl_cmd = ("rspctl {}".format(stats_flag_and_val)
                  + " --integration={}".format(int(integration))
                  + " --duration={}".format(int(duration)))
    if bsxtype == 'bst':
        rspctl_cmd += " --select=0,1"
    rspctl_cmds.append(rspctl_cmd)
    return rspctl_cmds

def parse_rspctl_args(rspctl_cmds):
    """Parse rspctl command arguments

    Note that rspctl has persistent flags, i.e. multiple rspctl calls add up
    flags.
    """
    # TODO: Add the rest of the arguments
    rspctl_args = {}
    rspctl_parser = argparse.ArgumentParser(prog='rspctl')
    rspctl_parser.add_argument('--statistics')
    rspctl_parser.add_argument('--xcstatistics', action='store_true')
    rspctl_parser.add_argument('--select')
    rspctl_parser.add_argument('--integration')
    rspctl_parser.add_argument('--duration')
    rspctl_parser.add_argument('--xcsubband')
    rspctl_parser.add_argument('--directory')
    rspctl_parser.add_argument('--bitmode')
    rspctl_parser.add_argument('--mode')
    for rspctl_cmd in rspctl_cmds:
        rspctl_cmd_normal = rspctl_cmd.replace('=', ' ')
        argsdict = vars(rspctl_parser.parse_args(rspctl_cmd_normal.split()[1:]))
        argsdict = {k: v for (k, v) in argsdict.items() if v is not None}
        rspctl_args.update(argsdict)
    return rspctl_args


def parse_lofar_conf_files(filetext):
    """
    Parse LOFAR .conf files
    
    Input file text, output dict of content.
    LOFAR .conf files have lines with format:

    configtype.key = value
    """
    contentdict = {}
    for line in filetext.split('\n'):
        if line.startswith('#'):
            continue
        if '=' in line:
            cfgtypedotkey, value = line.split('=')
            cfgtypedotkey = cfgtypedotkey.strip()
            value = value.strip()
            cfgtype, key = cfgtypedotkey.split('.')
            if cfgtype not in contentdict:
                contentdict[cfgtype] = {}
            contentdict[cfgtype][key] = value
    return contentdict


class FreqSetup(object):
    """\
    Class that handles frequency bands and all things related to them.

    Examples:
        Set using RCU band name:
    >>> from ilisa.operations.modeparms import FreqSetup
    >>> fbb=FreqSetup('10_90')
    >>> fbb.__dict__
    {'arg': '10_90', 'rcumodes': [3], 'sb_range': ['51:461'], 'bits': 8, 'rcubands': ['10_90'], 'antsets': ['LBA_INNER'], 'rcusel': ['0:191']}

    Set using single frequency, either string or float:
    >>> fb1=FreqSetup('180e6')
    >>> fb1.__dict__
    {'arg': '180e6', 'rcumodes': [5], 'sb_range': ['410'], 'bits': 16, 'rcubands': ['110_190'], 'antsets': ['HBA_JOINED'], 'rcusel': ['0:191']}

    Set using a tuple given by (freqlo, freqhi, freqstep) (internal rep):
    >>> fbt3=FreqSetup((150e6,220e6,5.0e6))
    >>> fbt3.__dict__
    {'arg': (150000000.0, 220000000.0, 5000000.0), 'rcumodes': [5, 7], 'sb_range': ['256,281,306,331,356,381,406,431,456', '51,76,101'], 'bits': 16, 'rcubands': ['110_190', '210_250'], 'antsets': ['HBA_JOINED', 'HBA_JOINED'], 'rcusel': ['0:47,96:143', '48:95,144:191']}

    Set using slice style index string given by 'freqlo:freqhi':
    >>> fbs2=FreqSetup('35e6:1e6:210e6')
    >>> fbs2.__dict__
    {'arg': '35e6:1e6:210e6', 'rcumodes': [4, 5, 7], 'sb_range': ['179,184,189,194,199,204,209,214,219,224,229,234,239,244,249,254,259,264,269,274,279,284,289,294,299,304,309,314,319,324,329,334,339,344,349,354,359,364,369,374,379,384,389,394,399,404,409,414,419,424,429,434,439,444,449,454,459', '51,56,61,66,71,76,81,86,91,96,101,106,111,116,121,126,131,136,141,146,151,156,161,166,171,176,181,186,191,196,201,206,211,216,221,226,231,236,241,246,251,256,261,266,271,276,281,286,291,296,301,306,311,316,321,326,331,336,341,346,351,356,361,366,371,376,381,386,391,396,401,406,411,416,421,426,431,436,441,446,451,456,461', '51'], 'bits': 16, 'rcubands': ['30_90', '110_190', '210_250'], 'antsets': ['LBA_INNER', 'HBA_JOINED', 'HBA_JOINED'], 'rcusel': ['0:31,96:127', '32:63,128:159', '64:95,160:191']}

    or using full slice style index string given by 'freqlo:freqstep:freqhi':
    >>> fbs3=FreqSetup('80e6:10.0e6:220e6')
    >>> fbs3.__dict__
    {'arg': '80e6:10.0e6:220e6', 'rcumodes': [4, 5, 7], 'sb_range': ['410,461', '51,102,153,204,255,306,357,408,459', '51,102'], 'bits': 16, 'rcubands': ['30_90', '110_190', '210_250'], 'antsets': ['LBA_INNER', 'HBA_JOINED', 'HBA_JOINED'], 'rcusel': ['0:31,96:127', '32:63,128:159', '64:95,160:191']}

    """
    # rcumode mappings (rcumode 1,2 not considered here)
    rcumode_smpfrqs = {3: 200e6,
                       4: 200e6,
                       5: 200e6,
                       6: 160e6,
                       7: 200e6}
    rcumode_passbands = {3: ( 10e6,  90e6),
                         4: ( 30e6,  90e6),
                         5: (110e6, 190e6),
                         6: (170e6, 230e6),
                         7: (210e6, 250e6)}
    rcumode_bandnames = {3: '10_90',
                         4: '30_90',
                         5: '110_190',
                         6: '170_230',
                         7: '210_250'}
    rcumode_antsets = {3: 'LBA_INNER',
                       4: 'LBA_INNER',
                       5: 'HBA_JOINED',
                       6: 'HBA_JOINED',
                       7: 'HBA_JOINED'}
    rcumode_nqzones = {3: 0,
                       4: 0,
                       5: 1,
                       6: 2,
                       7: 2}
    # rcus must be split into at least the same number of partitions as number
    # of spws:
    subarr_rcusel = [  # 1 spw (all rcus)
                     ['0:191'],
                       # 2 spws 2 subarrs
                     ['0:47,96:143', '48:95,144:191'],
                       # 3 spws 3 subarrs
                     ['0:31,96:127', '32:63,128:159', '64:95,160:191']]

    combos = [[4],[3],[5],[7],[6],[4,5],[3,5],[5,7],[4,5,7],[3,5,7]]
    nrffts = 1024

    def __init__(self, freq_arg=None, chan_arg={}):
        self.freq_arg = freq_arg
        if not freq_arg: return
        if type(freq_arg) == dict:
            chan_arg = freq_arg
        if not chan_arg:
            if '_' in freq_arg:
                # freq_arg has rcuband format 'rcubandlo_rcubandhi'
                freqbins = self._band2freqbins(freq_arg)
            else:
                # Normal arg spec
                freqbins = self._freqslice2freqbins(freq_arg)
            self._rcumodes, self.subbands_spw = self._subband_hint(*freqbins)
        else:
            self._rcumodes = [chan_arg['rcumode']]
            self.subbands_spw = [chan_arg['sbs']]
        _beamlets, _bmltpntr, nrbeamlets = alloc_beamlets(self.subbands_spw)
        if chan_arg.get('bits'):
            self.bits = chan_arg.get('bits')
        else:
            self.bits = bits_support_nrbeamlets(nrbeamlets)
        self.rcubands = []
        self.antsets = []
        for rcumode in self._rcumodes:
            self.rcubands.append(self.rcumode_bandnames[rcumode])
            self.antsets.append(self.rcumode_antsets[rcumode])
        self.rcusel = self.subarr_rcusel[len(self._rcumodes) - 1]
        self.nrlanes = math.ceil(
            nrbeamlets/NRBEAMLETSBYBITS[self.bits]*MAX_NRLANES
        )

    def _band2freqbins(self, bandarg):
        """Convert band to frequency bins specification"""
        freqlo, freqhi, freqstep = 0.0, 0.0, None
        for spw in self.rcumode_bandnames:
                bandname =self.rcumode_bandnames[spw]
                if bandarg == bandname:
                    (freqlo, freqhi) = self.rcumode_passbands[spw]
                    freqstep = None
                    break
        return (freqlo, freqhi, freqstep)

    def _freqslice2freqbins(self, freqslicearg):
        """
        Convert frequency slice string to frequency bins specification

        Returns
        -------
        tuple : (freqlo, freqhi, freqstep)
            Where freqlo and freqhi are always floats. freqstep can be float or
            None.
        """
        freqlo, freqhi = 0.0, 0.0
        freqstep = None
        if type(freqslicearg) is tuple:
            if len(freqslicearg) < 3:
                (freqlo, freqhi) = freqslicearg
                freqstep = None
            else:
            # Assume already in right format
                (freqlo, freqhi, freqstep) = freqslicearg
        else:
            # Can be float so convert to str 
            freqbandarg = str(freqslicearg)
            # Check if freqbandarg has format:
            #   'freqlo:freqstep:freqhi'|'freqlo:freqhi'|'freqcntr'
            freqbandlst = freqbandarg.split(':')
            if len(freqbandlst) == 0 or len(freqbandlst) > 3:
                raise ValueError('Wrong number of freqband slice args')
            if len(freqbandlst) == 3:
                (freqlo, freqstep, freqhi) = tuple(map(float, freqbandlst))
                return (freqlo, freqhi, freqstep)
            # Must be only 1 or 2 args in freqbandtuple argument
            # so complete the incomplete freqband tuple
            freqlo = float(freqbandlst.pop(0))
            freqhi = freqlo
            if len(freqbandlst):
                freqhi = float(freqbandlst.pop())
            freqstep = None
            if len(freqbandlst):
                freqstep = float(freqbandlst.pop())
        return (freqlo, freqhi, freqstep)

    def _subband_hint(self, freqlo, freqhi, freqstep=None):
        """
        Calculate subbands in respective rcumodes that fulfill given freqrange
        """
        # Determine which rcumodes would support desired freqband:
        rcumodes = self.supportrcumodes(freqlo, freqhi)
        
        if not rcumodes:
            # TODO Consider extending to sbs out of passband but still measured
            return None, None
        sampfreq, sampres = self.getsampinfo(rcumodes[0])

        sb_lo = self.freq2sb(freqlo, sampfreq)
        sb_hi = self.freq2sb(freqhi, sampfreq)
        sbstep = 1  # Default if freqstep None and freqlo != freqhi
        if freqstep:
            sbstep = self.freqstep2sbstep(freqstep, sampres)
        elif freqlo == freqhi:
            sbstep = 0

        # Compute subband ranges
        sb_range = []
        nrrcumodes = len(rcumodes)
        if sbstep == 0:
            sb_range = ["{}".format(sb_lo)]
        else:
            if nrrcumodes == 1:
                if sbstep == 1:
                    sb_range = ["{}:{}".format(sb_lo, sb_hi)]
                else:
                    sb_list = range(sb_lo, sb_hi + 1, sbstep)
                    sb_range = [','.join([str(el) for el in sb_list])]
            elif nrrcumodes == 2:
                sb_hi0 = self.freq2sb(self.rcumode_passbands[rcumodes[0]][1],
                                      sampfreq)
                sb_lo1 = self.freq2sb(self.rcumode_passbands[rcumodes[1]][0],
                                      sampfreq)
                if sbstep == 1:
                    sb_range = ["{}:{}".format(sb_lo, sb_hi0)]
                    sb_range.append("{}:{}".format(sb_lo1, sb_hi))
                else:
                    sb_list = range(sb_lo, sb_hi0 + 1, sbstep)
                    sb_range = [','.join([str(el) for el in sb_list])]
                    sb_list = range(sb_lo1, sb_hi + 1, sbstep)
                    sb_range.append(','.join([str(el) for el in sb_list]))
            elif nrrcumodes == 3:
                sb_hi0 = self.freq2sb(self.rcumode_passbands[rcumodes[0]][1],
                                      sampfreq)
                sb_lo1 = self.freq2sb(self.rcumode_passbands[rcumodes[1]][0],
                                      sampfreq)
                sb_hi1 = self.freq2sb(self.rcumode_passbands[rcumodes[1]][1],
                                      sampfreq)
                sb_lo2 = self.freq2sb(self.rcumode_passbands[rcumodes[2]][0],
                                      sampfreq)
                if sbstep == 1:
                    sb_range = ["{}:{}".format(sb_lo, sb_hi0)]
                    sb_range.append("{}:{}".format(sb_lo1, sb_hi1))
                    sb_range.append("{}:{}".format(sb_lo2, sb_hi))
                else:
                    sb_list = range(sb_lo, sb_hi0 + 1, sbstep)
                    sb_range = [','.join([str(el) for el in sb_list])]
                    sb_list = range(sb_lo1, sb_hi1 + 1, sbstep)
                    sb_range.append(','.join([str(el) for el in sb_list]))
                    sb_list = range(sb_lo2, sb_hi + 1, sbstep)
                    sb_range.append(','.join([str(el) for el in sb_list]))

        return rcumodes, sb_range

    def getsampinfo(self, rcumode):
        """
        Get sampling info (sampl. freq, freq resolution for rcumode given
        """
        sampfreq = self.rcumode_smpfrqs[rcumode]
        sampres = sampfreq/self.nrffts
        return sampfreq, sampres

    def freq2sb(self, freq, sampfreq):
        """Map frequency to subband"""
        # Using round() but could have other truncation strategy...
        abs_sb = int(round(freq / sampfreq * self.nrffts))
        sb = abs_sb % int(self.nrffts/2)
        return sb

    def rcumode_sb2freq(self, rcumode, sb):
        abs_sb = self.rcumode_nqzones[rcumode]*self.nrffts/2+sb
        return abs_sb*self.rcumode_smpfrqs[rcumode]/float(self.nrffts)

    def freqstep2sbstep(self,freqstep, sbwidth):
        return int(math.floor(freqstep/sbwidth))

    def supportrcumodes(self, freqlo, freqhi):
        """
        Return an rcumode combo that supports  given freqlo-freqhi interval
        """
        # Will assumes 100 MHz Nyquist freq here:
        sb_lo, nq_lo = freq2sb(freqlo)
        sb_hi, nq_hi = freq2sb(freqhi)
        spw_lo = nqz2rcumode(nq_lo)
        spw_hi = nqz2rcumode(nq_hi)
        if spw_lo == spw_hi:
            return [spw_lo]
        else:
            return self.passbandsupport_spws(freqlo, freqhi)

    def passbandsupport_spws(self, freqlo, freqhi):
        """
        Return spw combo that passband-support given freqlo-freqhi interval
        """
        for combo in self.combos:
            if self.rcumode_passbands[combo[0]][0] <= freqlo \
               and freqhi <= self.rcumode_passbands[combo[-1]][1]:
                return combo
        return None

    def getrcubandnames(self):
        """Return list of rcu band names."""
        return list(self.rcumode_passbands.values())

    def edgefreqs(self, spw = 0):
        """
        Return a tuple of lowest and highest frequency in frequency band of spw
        """
        sbs = seqarg2list(self.subbands_spw[spw])
        freqlo = self.rcumode_sb2freq(self._rcumodes[spw], sbs[0])
        freqhi = self.rcumode_sb2freq(self._rcumodes[spw], sbs[-1])
        return freqlo, freqhi

    def nrsubbands(self):
        """Returns total number of subbands in this FreqSetup object."""
        nrsbs = 0
        for sbrange in self.subbands_spw:
            nrsbs += len(seqarg2list(sbrange))
        return nrsbs


def elementMap2str(elmap):
    elmapStr = ""
    for el in elmap:
        elmapStr = elmapStr+hex(el).lstrip('0').lstrip('x')
    return elmapStr


def str2elementMap2(elmapstr):
    elmap = []
    for el in elmapstr:
        elmap.append(int(el, 16))
    return elmap


def seqarg2list(seqarg):
    """
    Return a list that is a sequence of integers corresponding to the
    sequence given by the (extended) commandline argument string
    """
    arglist=[]
    for el in seqarg.split(','):
        els = el.split(':')
        seqlo, seqhi = int(els[0]), int(els[-1])
        if len(els) == 3:
            seqstep = int(els[1])
        else:
            seqstep = 1
        arglist.extend(range(seqlo, seqhi+1,seqstep))
    return arglist


def list2seqarg(lst):
    """\
    Convert a list of integers into a sequence argument

    Parameters
    ----------
    lst : list
        List of arguments to be converted.

    Returns
    -------
    seqarg : str
        Sequence argument str.
    """
    nrseq0 = lst.pop(0)
    arglst = [nrseq0]
    lastnrseq = nrseq0
    for arg in lst:
        if arg - lastnrseq == 1:
            if arglst[-1] != ':':
                arglst.append(':')
            if arg == lst[-1]:
                arglst.append(arg)
        else:
            if arglst[-1] == ':':
                arglst.append(lastnrseq)
            arglst.append(',')
            arglst.append(arg)
        lastnrseq = arg
    seqarg = ''.join([str(elem) for elem in arglst])
    lst.insert(0, nrseq0)
    return seqarg


def band2rcumode(band):
    """Map band to rcumode string (Inverse of rcumode2band())"""
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


def band2antset_eu(band):
    """
    Map band to antenna-set for EU stations

    Parameters
    ----------
    band: str
        Band name (can be: '10_90', '30_90', '110_190', '170_230', '210_250').

    Returns
    -------
    antset: str
        Antenna-set name.
        For EU stations can be: 'LBA_INNER' or 'HBA_JOINED'
        (subset of modeparms.ANTENNA_SETS)
    """
    if band == "10_90" or band == "30_90":
        antset = "LBA_INNER"
    elif band == "110_190" or band == "170_230" or band == "210_250":
        antset = "HBA_JOINED"
    else:
        raise ValueError("Undefined band: {}.".format(band))
    return antset


def _pol_label2nr(pol_label):
    """\
    Return id number of polarization label

    Parameters
    ----------
    pol_nr : int
        Polarization channel id number: 0 or 1

    Returns
    -------
    pol_label : str
        Polarization channel label: 'X' or 'Y'

    Raises
    ------
    ValueError
        If pol_label is not 'X' or 'Y'
    """
    if pol_label == 'X':
        pol_nr = 0
    elif pol_label == 'Y':
        pol_nr = 1
    else:
        raise ValueError("Polarization label {} is neither 'X' nor 'Y'"
                         .format(pol_label))
    return pol_nr


def _pol_nr2label(pol_nr):
    """\
    Return name of pol number

    Parameters
    ----------
    pol_nr : int
        Polarization channel id: 0 or 1

    Returns
    -------
    pol_label : str

    Raises
    ------
    ValueError
        If pol_nr is not 0 or 1
    """
    if pol_nr == 0:
        pol_label = 'X'
    elif pol_nr == 1:
        pol_label = 'Y'
    else:
        raise ValueError("Polarization number {} is neither 0 nor 1"
                         .format(pol_nr))
    return pol_label


def antpol2rcu(antnr, pol_label):
    """\
    Determine RCU# from antenna number and polarization channel

    Parameters
    ----------
    antnr : int
        Antenna number
    pol_label : str
        Polarization channel label: 'X' or 'Y'

    Returns
    -------
    rcunr : int
        RCU number
    """
    pol_nr = _pol_label2nr(pol_label)
    rcunr = 2 * antnr + pol_nr
    return rcunr


def rcu2antpol(rcunr):
    """\
    Convert RCU# to antenna # and polarization channel

    Parameters
    ----------
    rcunr : int
        RCU number

    Returns
    -------
    antnr : int
        Antenna number
    polchan : str
        Polarization channel: 'X' or 'Y'
    """
    antnr, pol_nr = divmod(rcunr, 2)
    pol_label = _pol_nr2label(pol_nr)
    return antnr, pol_label


def rcumode2sbfreqs(rcumode):
    """
    Get the frequencies (in Hz) of the subbands for the given rcumode

    Parameters
    ----------
    rcumode: int
        The RCUmode.

    Returns
    -------
    freqs: array of floats
        Frequencies in Hz. Array index corresponds to subband number.
    """
    nqzone = (int(rcumode)-3)/2
    # Note the endpoint=False here. Before it 2018-03-22 it was missing.
    freqs = numpy.linspace(nqzone * NQFREQ_NOM, (nqzone + 1) * NQFREQ_NOM,
                           TotNrOfsb, endpoint=False)
    return freqs


def rcumode2band(rcumode):
    """
    Map rcumode to band string as used in beamctl arguments

    Parameters
    ----------
    rcumode: int or str
        The RCU mode.
    Returns
    -------
    band: str
        The band name.
    """
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


def rcumode2antset_eu(rcumode):
    """
    Map rcumode to antenna set for EU stations.

    Antenna set is only meaningful in Dutch stations.
    But it is required in beamctl arguments.

    Parameters
    ----------
    rcumode: int
        The RCU mode.

    Returns
    -------
    antset: str
        Antenna set name.
    """
    rcumode = int(rcumode)
    if rcumode == 3 or rcumode == 4:
        antset = 'LBA_INNER'
    elif rcumode == 5 or rcumode == 6 or rcumode == 7:
        antset = 'HBA_JOINED'
    else:
        raise ValueError("Undefined rcumode: {}.".format(rcumode))
    return antset


def rcumode2nyquistzone(rcumode):
    """\
    Compute Nyquist zone for given RCU mode.

    Parameters
    ----------
    rcumode: int or str
        The RCU mode (a.k.a SPW).

    Returns
    -------
    nqzone: int
        Nyquist zone number.
    """
    nqzone = int((int(rcumode)-3)/2)
    return nqzone


def freq2sb(freq):
    """\
    Convert frequency in Hz to subband number and Nyquist zone

    Parameters
    ----------
    freq: float
        Frequency in Hz

    Returns
    -------
    sb: int
        Subband number
    nqzone: int
        Nyquist zone number
    """
    abs_sb = int(round(freq / NQFREQ_NOM * TotNrOfsb))
    sb = abs_sb % TotNrOfsb
    nqzone = abs_sb // TotNrOfsb
    return sb, nqzone


def sb2freq(sb, nqzone):
    """
    Convert subband in a given Nyquist zone to a frequency

    Parameters
    ----------
    sb: int or str
        Subband number
    nqzone: int or str
        Nyquist zone number

    Returns
    -------
    freq: float
        Frequency in Hz
    """
    freq = NQFREQ_NOM * (int(sb) / float(TotNrOfsb) + int(nqzone))
    return freq


def nqz2rcumode(nqzone, nqfreq=NQFREQ_NOM, filt_on=False):
    """\
    Convert Nyquist zone number and Nyquist frequency to RCU mode.

    Parameters
    ----------
    nqzone: int
        Nyquist zone number.
    nqfreq: float, optional
        Nyquist frequency in Hz. Default 100 MHz.
    filt_on: bool, optional.
        Is narrower band-pass filter on? Default False.

    Returns
    -------
    rcumode: int
        RCU mode.

    Raises
    ------
    ValueError
        If nqzone is out of range for nqfreq.
    """
    if nqfreq == NQFREQ_NOM:
        if nqzone == 0:
            if not filt_on:
                rcumode = 3
            else:
                rcumode = 4
        elif nqzone == 1:
            rcumode = 5
        elif nqzone == 2:
            rcumode = 7
        else:
            raise ValueError('nqzone out of range for 100MHz NQ')
    else:
        if nqzone == 1:
            rcumode = 6
        else:
            raise ValueError('nqzone out of range for non 100MHz NQ')
    return rcumode


def freq2bandarr(freq):
    """
    Frequency to band-array

    Parameters
    ----------
    freq: float
        Frequency in Hz.
    Returns
    -------
    bandarr : str
        Either LBA or HBA.
    """
    if freq < 100e6:
        bandarr = 'LBA'
    else:
        bandarr = 'HBA'
    return bandarr


def dt2mjd(dt):
    """
    Convert a python datetime to modified julian date

    Parameters
    ----------
    dt: datetime
        Datetime

    Returns
    -------
    mjd: float
        Modified Julian date
    """
    sigdec = 5
    ymd, hms = dt.date(), dt.time()
    dayfrac = (datetime.datetime.combine(datetime.date.min, hms)
                 - datetime.datetime.min).total_seconds()/86400.0
    ordinalfloat = ymd.toordinal() + dayfrac
    jd = round(ordinalfloat + 1721424.5, sigdec)
    mjd = round(jd - 2400000.5, sigdec)
    return mjd


def alloc_beamlets(subbands, bmlt_start=0):
    """
    Allocate subbands to beamlets sequentially over all spws

    Parameters
    ----------
    subbands : list
        List of subband sequence args.
    bmlt_start : int
        Beamlet number to start from.

    Returns
    -------
    beamlets : list | int
        List of beamlets.
    bmlt_pntr : int
        Beamlet number that would be allocated next.
    nrbmlts : int
        Number of beamlets allocated.
    """
    beamlets = []
    bmlt_pntr = bmlt_start
    singleton = False
    if type(subbands) is not list:
        singleton = True
        subbands = [subbands]
    for spw_sb in subbands:
        nrsbs = len(seqarg2list(spw_sb))
        bmltslice = str(bmlt_pntr)
        if nrsbs > 1:
            bmltslice += ':'+str(bmlt_pntr+nrsbs-1)
        beamlets += [bmltslice]
        bmlt_pntr += nrsbs
    nrbmlts = bmlt_pntr - bmlt_start
    if singleton:
        beamlets = beamlets.pop()
    return beamlets, bmlt_pntr, nrbmlts


def getlanes(subbands_spw, bits, nrlanes):
    """\
    Return a dict keyed on lanenr whose value is the beamlets allocated
    """
    max_bmlts_per_lane = round(NRBEAMLETSBYBITS[bits]/MAX_NRLANES)
    lanesplitbmlt = iter([(lanenr, (lanenr+1)*max_bmlts_per_lane-1)
                          for lanenr in range(nrlanes)])
    bmlts = []
    beamlets, _lstbmlt, _nrbmlts = alloc_beamlets(subbands_spw)
    for bmltarg in beamlets:
        bmlts.extend(seqarg2list(bmltarg))
    lanealloc = {0:[]}
    (lanenr, bmlt_hi) = next(lanesplitbmlt)
    for bmlt in bmlts:
        if bmlt <= bmlt_hi:
            lanealloc[lanenr].append(bmlt)
        else:
            lanealloc[lanenr+1] = []
            lanealloc[lanenr+1].append(bmlt)
            (lanenr, bmlt_hi) = next(lanesplitbmlt)
    return lanealloc


def bits_support_beamlets(beamlets):
    """
    Compute maximum number of bits that supports the beamlets given
    """
    # Compute total number of beamlets
    nrbeamlets = 0
    for spw_beamlet in beamlets:
        if ':' in spw_beamlet:
            bl_lo, bl_hi = spw_beamlet.split(':')
            nrbeamlets += int(bl_hi) - int(bl_lo) + 1
        else:
            nrbeamlets += len(spw_beamlet.split(','))
    # Compute maximum nr of bits that fulfills number of beamlets 
    return bits_support_nrbeamlets(nrbeamlets)


def bits_support_nrbeamlets(nrbeamlets):
    """
    Compute maximum number of bits that supports the given number of beamlets
    """
    # Compute maximum nr of bits that fulfills number of beamlets 
    bits = 16  # Default nr of bits
    for bits in NRBEAMLETSBYBITS.keys():
        if nrbeamlets <= NRBEAMLETSBYBITS[bits]:
            break
    return bits


def nrrcus_stnid(stnid):
    """
    Number of RCUs from station ID

    Parameters
    ---------
    stnid: str
        Station ID

    Returns
    -------
    nrrcus: int
        Number of RCUs
    """
    location_id = stnid[:2]
    if location_id == 'CS' or location_id == 'RS' or location_id == 'FI':
        nrrcus = 96
    else:
        # EU station
        nrrcus = 192
    return nrrcus


def timestr2datetime(timestr, dt_format=DATETIMESTRFMT):
    """\
    Convert into a python datetime

    Parameters
    ----------
    timestr: str or datetime.datetime
        Date-Time string in ISO-like format '%Y-%m-%dT%H:%M:%S'
        OR 'NOW' or 'ASAP', which imply the current UT datetime.
        If timestr is a datetime object, it is returned as is.
    dt_format : str
        Format string to parse input timestr
        (syntax as defined in python datetime package).

    Returns
    -------
    dattim: datetime.datetime
        Python datetime object corresponding to input.
    """
    if type(timestr) == datetime.datetime: return timestr
    if timestr == 'NOW' or timestr == 'ASAP':
        # Set time to nearest rounded second from now:
        dattim = datetime.datetime.utcnow()
        dattim = dattim.replace(microsecond=0)
        dattim += datetime.timedelta(seconds=1)
    else:
        try:
            dattim = datetime.datetime.strptime(timestr, dt_format)
        except:
            raise RuntimeError("Wrong datetime format.")
    return dattim


def as_asapdatetime(timestr, dt_format=DATETIMESTRFMT):
    """\
    Convert timestr to datetime but keep ASAP
    """
    if timestr == 'ASAP': return 'ASAP'
    return timestr2datetime(timestr, dt_format)


def astimestr(datim, dt_format=DATETIMESTRFMT):
    """\
    Convert datetime object or 'ASAP' to datetime string or 'ASAP' resp.

    Parameters
    ----------
    datim : datetime or 'ASAP'
        The datetime or 'ASAP' to be converted.
    dt_format : str
        python datetime format.

    Returns
    -------
    timestr : str
        Either 'ASAP' or the string formatted datim.
    """
    if datim == 'ASAP':
        timestr = 'ASAP'
    else:
        timestr = datim.strftime(dt_format)
    return timestr


def normalizetimestr(timestr, dt_pformat=DATETIMESTRFMT,
                     dt_fformat=DATETIMESTRFMT):
    """\
    Normalize iLiSA time str

    Parameters
    ----------
    timestr :  str or datetime.datetime
        The date-time to convert. If `timestr` is a string that is datetime
        formatted, then `dt_pformat` is used to parse it. If it is the str
        'ASAP' it is evaluated as the date-time now. Special case is if
        `timestr` is a datetime.datetime.

    Returns
    -------
    timestr_nrm : str
        Date-time string in the format specified by `dt_fformat`.
    """
    if type(timestr) != datetime.datetime:
        datim = timestr2datetime(timestr, dt_pformat)
    else:
        datim = timestr
    timestr_nrm = astimestr(datim, dt_fformat)
    return timestr_nrm


def is_ldattype(obsdatatype):
    """
    Test if a string 'obsdatatype' is a type of LOFAR data

    Parameters
    ----------
    obsdatatype: str
        Acronym for an LDAT type. Can be: 'acc', 'bst', 'sst', 'xst' or 'bfs'.

    Returns
    -------
    _isldattype: bool
        If True, `obsdatatype` is an LDAT-type, else it is not.
    """
    _isldattype = False
    if (obsdatatype == 'acc' or
            obsdatatype == 'bst' or
            obsdatatype == 'bst-357' or
            obsdatatype == 'sst' or
            obsdatatype == 'xst' or
            obsdatatype == 'xst-SEPTON' or
            obsdatatype == 'bfs'):
        _isldattype = True
    return _isldattype


def parse_ldattype_datetime_from_filename(filenamepath):
    """
    Parse LDAT-type and datetime from LDAT filename

    The data files produced on the LCU (LDAT) have the format

        20221214_101419_xst.dat

    which consists of a date and UT followed by a lofar data type with
    file extension '.dat'.

    Parameters
    ----------
    filename: str
        LDAT filename

    Returns
    -------
    dattim: datetime
        Datetime corresponding to ldat filename.
    ldattype: str
        LDAT type
    """
    filename = os.path.basename(os.path.normpath(filenamepath))
    filenamebase, ext = os.path.splitext(filename)
    if ext != '.dat':
        warnings.warn(
            "LDAT name not in expected (extension was {} rather than '-dat')"
                .format(ext))
    Ymd, HMS, ldattype = filenamebase.split('_')
    if not is_ldattype(ldattype):
        warnings.warn("'{}' is not valid LDAT type.".format(ldattype))
    dattim = datetime.datetime.strptime(Ymd+'_'+HMS, '%Y%m%d_%H%M%S')
    return dattim, ldattype


def format_filename_from_ldattype_datetime(dattim, ldattype):
    """
    Return the LDAT formatted file name for LDAAT-type and datetime

    Parameters
    ----------
    dattim: datetime
        Date-time of file.
    ldattype:
        LDAT-type, e.g. 'bst','sst' or 'xst'.

    Returns
    -------
    filename: str
        The file name.
    """
    if not is_ldattype(ldattype):
        warnings.warn("'{}' is not valid LDAT type.".format(ldattype))
    Ymd_HMS = dattim.strftime('%Y%m%d_%H%M%S')
    filename = "{}_{}.dat".format(Ymd_HMS, ldattype)
    return filename


def hmsstr2deltatime(hms):
    """\
    Convert hms (hours, minutes, seconds) string to python deltatime

    Parameters
    ----------
    hms: str
        Duration string in format '%Hh%Mm%Ss', '%Hh', '%Mm', '%Ss',
        or if none of these it will interpret hms as seconds.

    Returns
    -------
    delta_t: datetime.deltatime
        Python deltatime object corresponding to input.
    """
    parse_format = ''
    if 'h' in hms and 'm' in hms and 's' in hms:
        parse_format = '%Hh%Mm%Ss'
    elif 'h' in hms and 'm' not in hms and 's' not in hms:
        parse_format = '%Hh'
    elif 'h' not in hms and 'm' in hms and 's' not in hms:
        parse_format = '%Mm'
    elif 'h' not in hms and 'm' not in hms and 's' in hms:
        parse_format = '%Ss'
    if parse_format != '':
        dt = datetime.datetime.strptime(hms, parse_format)
        delta_t = datetime.timedelta(hours=dt.hour, minutes=dt.minute,
                                     seconds=dt.second)
    else:
        delta_t = datetime.timedelta(seconds=float(eval(hms)))
    return delta_t


def modelogic(freqspec, pointing, duration_tot, integration, bsx_type, bfs, acc,
              allsky):
    """Determine if the mode parameters given make sense"""
    if acc and not pointing:
        warnings.warn('ACC requires beamctl')
    # TODO add more conditions


def _xtract_bsx(rec_type):
    """Parse rec_type to bsx_type"""
    bsx_type = None
    if not rec_type or rec_type == 'None':
        # rec_type None means running a beam with no recording
        bsx_type = None
    elif rec_type == 'bst' or rec_type == 'sst' or rec_type == 'xst':
        bsx_type = rec_type
    elif rec_type == 'bfs' or rec_type == 'tbb' or rec_type == 'dmp':
        # 'dmp' is for just recording without setting up a beam.
        pass
    else:
        raise RuntimeError('Unknown rec_type {}'.format(rec_type))
    return bsx_type