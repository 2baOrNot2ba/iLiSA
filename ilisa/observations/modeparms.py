"""This package is for the parameters involved in observation modes.
"""
import argparse
import math
import numpy
import datetime

rcusbsep = "+"
Nqfreq = 100.0e6  # Nyquist frequency in Hz
TotNrOfsb = 512  # Total number of subbands. (Subbands numbered 0:511)
nrofrcus = 192  # Number of RCUs
MIN_STATS_INTG = 1.0  # Minimum integration for statistics data in seconds.
BASE_NR_BEAMLETS = 244

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


def parse_rspctl_args(rspctl_strs):
    """Parse rspctl command arguments.
    Note that rspctl has persistent flags, i.e. multiple rspctl calls add up flags."""
    # TODO: Add the rest of the arguments
    rspctl_args ={}
    rspctl_parser = argparse.ArgumentParser()
    rspctl_parser.add_argument('--statistics')
    rspctl_parser.add_argument('--xcstatistics', action='store_true')
    rspctl_parser.add_argument('--integration')
    rspctl_parser.add_argument('--duration')
    rspctl_parser.add_argument('--xcsubband')
    rspctl_parser.add_argument('--directory')
    rspctl_parser.add_argument('--bitmode')
    rspctl_strs = rspctl_strs.lstrip('; ')
    for rspctl_line in rspctl_strs.split('\n'):
        for rspctl_str in rspctl_line.split(';'):
            rspctl_str_normalized = rspctl_str.replace('=', ' ')
            argsdict = vars(rspctl_parser.parse_args(rspctl_str_normalized.split()[1:]))
            argsdict = { k:v for (k,v) in argsdict.items() if v is not None }
            rspctl_args.update(argsdict)
    return rspctl_args

def parse_lofar_conf_files(filetext):
    """Parse LOFAR .conf files. Input file text, output dict of content.
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


class FrequencyBand(object):
    """Class that handles frequency bands and all things related to them.

    Examples:
        Set using RCU band name:
    >>> from ilisa.observations.modeparms import FrequencyBand
    >>> fbb=FrequencyBand('10_90')
    >>> fbb.__dict__
    {'antsets': ['LBA_INNER'],
     'arg': '10_90',
     'beamlets': ['0:410'],
     'bits': 8,
     'rcubands': ['10_90'],
     'rcumodes': [3],
     'rcusel': ['0:191'],
     'sb_range': ['51:461']}

    Set using single frequency, either string or float:
    >>> fb1=FrequencyBand('180e6')
    >>> fb1.__dict__
    {'antsets': ['HBA_JOINED'],
     'arg': '180e6',
     'beamlets': ['0'],
     'bits': 16,
     'rcubands': ['110_190'],
     'rcumodes': [5],
     'rcusel': ['0:191'],
     'sb_range': ['410']}

    Set using a tuple given by (freqlo, freqhi):
    >>> fbt2=FrequencyBand((40e6,120e6))
    >>> fbt2.__dict__
    {'antsets': ['LBA_INNER', 'HBA_JOINED'],
     'arg': (40000000.0, 120000000.0),
     'beamlets': ['0:256', '257:308'],
     'bits': 8,
     'rcubands': ['30_90', '110_190'],
     'rcumodes': [4, 5],
     'rcusel': ['0:47,96:143', '48:95,144:191'],
     'sb_range': ['205:461', '51:102']}

    Set using a tuple given by (freqlo, freqstep, freqhi):
    >>> fbt3=FrequencyBand((150e6,5.0e6,220e6))
    >>> fbt3.__dict__
    {'antsets': ['HBA_JOINED', 'HBA_JOINED'],
     'arg': (150000000.0, 5000000.0, 220000000.0),
     'beamlets': ['0:8', '9:11'],
     'bits': 16,
     'rcubands': ['110_190', '210_250'],
     'rcumodes': [5, 7],
     'rcusel': ['0:47,96:143', '48:95,144:191'],
     'sb_range': ['256,281,306,331,356,381,406,431,456', '51,76,101']}

    Set using slice style index string given by 'freqlo:freqhi':
    >>> fbs2=FrequencyBand('35e6:1e6:210e6')
    >>> fbs2.__dict__
    {'antsets': ['LBA_INNER', 'HBA_JOINED', 'HBA_JOINED'],
     'arg': '35e6:210e6',
     'beamlets': ['0:282', '283:693', '694:694'],
     'bits': 4,
     'rcubands': ['30_90', '110_190', '210_250'],
     'rcumodes': [4, 5, 7],
     'rcusel': ['0:31,96:127', '32:63,128:159', '64:95,160:191'],
     'sb_range': ['179:461', '51:461', '51:51']}

     or using full slice style index string given by 'freqlo:freqstep:freqhi':
     >>> fbs3=FrequencyBand('80e6:10.0e6:220e6')
     >>> fbs3.__dict__
    {'antsets': ['LBA_INNER', 'HBA_JOINED', 'HBA_JOINED'],
     'arg': '80e6:10e6:220e6',
     'beamlets': ['0:1', '2:10', '11:12'],
     'bits': 16,
     'rcubands': ['30_90', '110_190', '210_250'],
     'rcumodes': [4, 5, 7],
     'rcusel': ['0:31,96:127', '32:63,128:159', '64:95,160:191'],
     'sb_range': ['410,461', '51,102,153,204,255,306,357,408,459', '51,102']}

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
    # rcus must be split into at least the same number of partitions as number of spws:
    subarr_rcusel = [['0:191'],                                         # 1 spw (all rcus)
                     ['0:47,96:143', '48:95,144:191'],                  # 2 spws 2 subarrs
                     ['0:31,96:127', '32:63,128:159', '64:95,160:191']  # 3 spws 3 subarrs
                     ]

    combos = [[4],[3],[5],[7],[6],[4,5],[3,5],[5,7],[4,5,7],[3,5,7]]

    nrffts = 1024
    nrbeamletsbybits = {16:   BASE_NR_BEAMLETS,
                        8:  2*BASE_NR_BEAMLETS,
                        4:  4*BASE_NR_BEAMLETS}
    nrlanes = 4

    def __init__(self, arg=None):
        self.arg = arg
        if arg is not None:
            if type(arg) is not tuple:
                freqrange = self._freqbandarg2tuple(arg)
            else:
                freqrange = arg
            self.rcumodes, self.rcubands, self.antsets, self.sb_range, \
                self.bits, self.beamlets, self.rcusel \
                = self.find_obsctlparams(freqrange)

    def _freqbandarg2tuple(self, freqbandarg):
        """Convert freqbandarg to a tuple. The tuple is either (freqlo, freqstep, freqhi)
        or (freqlo, freqhi) or (freqcntr,)."""
        if type(freqbandarg) is str:
            # Check if freqbandarg has a rcuband format 'rcubandlo_rcubandhi'
            if '_' in freqbandarg:
                band = None
                for rcumode in self.rcumode_bandnames.keys():
                    if freqbandarg == self.rcumode_bandnames[rcumode]:
                        band = freqbandarg
                        freqlo, freqhi = self.rcumode_passbands[rcumode]
                        break
                if band is None:
                    raise ValueError("Unknown rcu band name: {}.".format(freqbandarg))
                return (freqlo, freqhi)
            # Check if freqbandarg has format:
            #   'freqlo:freqstep:freqhi'|'freqlo:freqhi'|'freqcntr'
            try:
                freqrange = eval('({},)'.format(freqbandarg.replace(':',',')))
            except:
                raise ValueError("Wrong frequency band slice format.")
            return freqrange
        elif type(freqbandarg) is float:
            freqcntr = freqbandarg
            return (freqcntr,)
        else:
            raise ValueError("Wrong frequency band format.")

    def find_obsctlparams(self, freqrange):
        """Find and set observation ctrl parameters related to freqband."""
        nrfrqbndterms = len(freqrange)
        if nrfrqbndterms == 3:
            freqlo, freqstep, freqhi = freqrange
        elif nrfrqbndterms == 2:
            freqlo, freqhi = freqrange
        elif nrfrqbndterms == 1:
            (freqcntr,) = freqrange
            freqlo, freqhi = (freqcntr, freqcntr)

        # Determine which rcumodes would support desired freqband:
        rcumodes = self.supportrcumodes(freqlo, freqhi)
        rcubands = []
        antsets = []
        if not rcumodes:
            # TODO Consider extending to sbs outside of passband but still measured
            return None, None, None, None, None, None, None
        sampfreq, sampres = self.getsampinfo(rcumodes[0])

        for rcumode in rcumodes:
            rcubands.append(self.rcumode_bandnames[rcumode])
            antsets.append(self.rcumode_antsets[rcumode])
        if nrfrqbndterms == 1:
            sb_lo = self.freq2sb(freqcntr, sampfreq)
            sb_hi = sb_lo
            sbstep = 0
        else:
            sb_lo = self.freq2sb(freqlo, sampfreq)
            sb_hi = self.freq2sb(freqhi, sampfreq)
            if nrfrqbndterms == 3:
                sbstep = self.freqstep2sbstep(freqstep, sampres)
            else:
                sbstep = 1

        # Compute subband ranges
        nrrcumodes = len(rcumodes)
        if sbstep == 0:
            sb_range = ["{}".format(sb_lo)]
            beamlets = ['0']
        else:
            if nrrcumodes == 1:
                if sbstep == 1:
                    sb_range = ["{}:{}".format(sb_lo, sb_hi)]
                    nrsbs0 = sb_hi - sb_lo + 1
                else:
                    sb_list = range(sb_lo, sb_hi + 1, sbstep)
                    sb_range = [','.join([str(el) for el in sb_list])]
                    nrsbs0 =  len(sb_list)
                beamlets = ["0:{}".format(nrsbs0-1)]
            elif nrrcumodes == 2:
                sb_hi0 = self.freq2sb(self.rcumode_passbands[rcumodes[0]][1],
                                      sampfreq)
                sb_lo1 = self.freq2sb(self.rcumode_passbands[rcumodes[1]][0],
                                      sampfreq)
                if sbstep == 1:
                    sb_range = ["{}:{}".format(sb_lo, sb_hi0)]
                    nrsbs0 = sb_hi0 - sb_lo + 1
                    sb_range.append("{}:{}".format(sb_lo1, sb_hi))
                    nrsbs1 = sb_hi - sb_lo1 + 1
                else:
                    sb_list = range(sb_lo, sb_hi0 + 1, sbstep)
                    sb_range = [','.join([str(el) for el in sb_list])]
                    nrsbs0 = len(sb_list)
                    sb_list = range(sb_lo1, sb_hi + 1, sbstep)
                    sb_range.append(','.join([str(el) for el in sb_list]))
                    nrsbs1 = len(sb_list)
                beamlets = ["0:{}".format(nrsbs0 - 1)]
                beamlets.append("{}:{}".format(nrsbs0, nrsbs0 + nrsbs1 - 1))
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
                    nrsbs0 = sb_hi0 - sb_lo + 1
                    sb_range.append("{}:{}".format(sb_lo1, sb_hi1))
                    nrsbs1 = sb_hi1 - sb_lo1 + 1
                    sb_range.append("{}:{}".format(sb_lo2, sb_hi))
                    nrsbs2 = sb_hi - sb_lo2 + 1
                else:
                    sb_list = range(sb_lo, sb_hi0 + 1, sbstep)
                    sb_range = [','.join([str(el) for el in sb_list])]
                    nrsbs0 = len(sb_list)
                    sb_list = range(sb_lo1, sb_hi1 + 1, sbstep)
                    sb_range.append(','.join([str(el) for el in sb_list]))
                    nrsbs1 = len(sb_list)
                    sb_list = range(sb_lo2, sb_hi + 1, sbstep)
                    sb_range.append(','.join([str(el) for el in sb_list]))
                    nrsbs2 = len(sb_list)
                beamlets = ["0:{}".format(nrsbs0 - 1)]
                beamlets.append("{}:{}".format(nrsbs0, nrsbs0 + nrsbs1 - 1))
                beamlets.append("{}:{}".format(nrsbs0+nrsbs1, nrsbs0+nrsbs1+nrsbs2-1))
        nrbeamlets = 0
        for rcumode_beamlet in beamlets:
            if ":" in rcumode_beamlet:
                bl_lo, bl_hi = rcumode_beamlet.split(':')
                nrbeamlets += int(bl_hi) -int(bl_lo) + 1
            else:
                nrbeamlets += len(rcumode_beamlet.split(','))
        if nrbeamlets > self.nrbeamletsbybits[4]:
            raise ValueError("Requires too many beamlets ({} > {})."
                             .format(nrbeamlets, self.nrbeamletsbybits[4]))
        for bits in self.nrbeamletsbybits.keys():
            # Choose largest bit size that fulfills nrbeamlets
            if nrbeamlets <= self.nrbeamletsbybits[bits]:
                break
        rcusel = self.subarr_rcusel[nrrcumodes-1]
        return rcumodes, rcubands, antsets, sb_range, bits, beamlets, rcusel

    def getsampinfo(self,rcumode):
        sampfreq = self.rcumode_smpfrqs[rcumode]
        sampres = sampfreq/self.nrffts
        return sampfreq, sampres

    def freq2sb(self, freq, sampfreq):
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
        for combo in self.combos:
            if self.rcumode_passbands[combo[0]][0] <= freqlo \
               and freqhi <= self.rcumode_passbands[combo[-1]][1]:
                return combo
        return None

    def getrcubandnames(self):
        """Return list of rcu band names."""
        return list(self.rcumode_passbands.values())

    def edgefreqs(self, spw = 0):
        """Return a tuple of lowest and highest frequency in frequency band of spw."""
        sbs = seqarg2list(self.sb_range[spw])
        freqlo = self.rcumode_sb2freq(self.rcumodes[spw], sbs[0])
        freqhi = self.rcumode_sb2freq(self.rcumodes[spw], sbs[-1])
        return freqlo, freqhi

    def nrsubbands(self):
        """Returns total number of subbands in this FrequencyBand object."""
        nrsbs = 0
        for sbrange in self.sb_range:
            nrsbs += len(seqarg2list(sbrange))
        return nrsbs

    def getlanes(self):
        """Return a dict keyed on lanenr whose value is the beamlets allocated."""
        bmlts_per_lane = self.nrbeamletsbybits[self.bits]/self.nrlanes
        lanesplitblmt = iter([(lanenr, (lanenr+1)*bmlts_per_lane-1) for lanenr in
                              range(self.nrlanes)])
        bmlts = []
        for bmltarg in self.beamlets:
            bmlts.extend(seqarg2list(bmltarg))
        lanealloc = {0:[]}
        (lanenr, bmlt_hi) = next(lanesplitblmt)
        for bmlt in bmlts:
            if bmlt <= bmlt_hi:
                lanealloc[lanenr].append(bmlt)
            else:
                lanealloc[lanenr+1] = []
                lanealloc[lanenr+1].append(bmlt)
                (lanenr, bmlt_hi) = next(lanesplitblmt)
        return lanealloc

    def get_maxbeamletsbybits(self, bits=None):
        """Return the maximum number of beamlets depending on bit depth."""
        if not bits:
            try:
                bits = self.bits
            except NameError:
                raise ValueError("Need to specify bits.")
        return self.nrbeamletsbybits[bits]


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
    """Return a list that is a sequence of integers corresponding to the sequence given
    by the (extended) commandline argument string."""
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


def seqlists2slicestr(seqlists):
    """Convert a sequence list to slice format"""
    # Instead of comma separated list format (e.g. 202,204,206),
    # try to construct subband slice syntax (e.g. 202:2:206), if
    # possible, to avoid file names that are potentially longer than 255
    # (not allowed in linux filesyst)
    def seqlist2slice(seqlist):
        seqlistcanon = []
        for seqel in seqlist.split(','):
            seqel = [int(el) for el in seqel.split(':')]
            seq = range(seqel[0], seqel[-1]+1)
            seqlistcanon.extend(seq)
        seqsteps = set(numpy.diff(seqlistcanon))
        if len(seqsteps) > 1:
            raise ValueError('Subband spec too complicated.')
        elif len(seqsteps) == 0:
            slicestr = "{}".format(seqlistcanon[0])
        else:
            seqstep = seqsteps.pop()
            seqstepstr = str(seqstep) + ':' if seqstep > 1 else ''
            slicestr = "{}:{}{}".format(seqlistcanon[0], seqstepstr, seqlistcanon[-1])
        return slicestr

    if type(seqlists) is list:
        slicestrlist = []
        for seqlist in seqlists:
            seqstr = seqlist2slice(seqlist)
            slicestrlist.append(seqstr)
        slicestr = rcusbsep.join(slicestrlist)
    else:
        slicestr = seqlist2slice(seqlists)
    return slicestr


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


def rcumode2sbfreqs(rcumode):
    """Get the frequencies (in Hz) of the subbands for the given rcumode.
    Returns an array of frequencies where index is subband number."""
    NZ = (int(rcumode)-3)/2
    # Note the endpoint=False here. Before it 2018-03-22 it was missing.
    freqs = numpy.linspace(NZ*Nqfreq, (NZ+1)*Nqfreq, TotNrOfsb, endpoint=False)
    return freqs


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


def rcumode2antset(rcumode):
    """Map rcumode to antennaset, which is used in beamctl arguments.
    Assumption is that one wants to use as many of antennas in field as
    possible. (This function may soon be deprecated.)
    """
    # NOTE new/more antennasets are now available.
    rcumode = int(rcumode)
    if rcumode == 3 or rcumode == 4:
        antset = 'LBA_INNER'
    elif rcumode == 5 or rcumode == 6 or rcumode == 7:
        antset = 'HBA_JOINED'
    else:
        raise ValueError("Undefined rcumode: {}.".format(rcumode))
    return antset


def rcumode2nyquistzone(rcumode):
    nz = int((int(rcumode)-3)/2)
    return nz


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

def dt2mjd(dt):
    """Convert a python datetime to modified julian date."""
    sigdec = 5
    ymd, hms = dt.date(), dt.time()
    dayfrac = (datetime.datetime.combine(datetime.date.min, hms)
                 - datetime.datetime.min).total_seconds()/86400.0
    ordinalfloat = ymd.toordinal() + dayfrac
    jd = round(ordinalfloat + 1721424.5, sigdec)
    mjd = round(jd - 2400000.5, sigdec)
    return mjd
