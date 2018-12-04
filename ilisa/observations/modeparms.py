"""This package is for the parameters involved in observation modes."""
import math
import numpy

from ilisa.observations.stationcontrol import Nqfreq, TotNrOfsb


class FrequencyBand(object):
    """Class that handles frequency bands and all things related to them.

    Examples:
        Set using RCU band name:
    >>> from ilisa.observations.modeparms import FrequencyBand
    >>> fbb=FrequencyBand('10_90')
    >>> fbb.__dict__
    {'antsets': ['LBA_INNER'],
     'beamlets': ['0:410'],
     'bits': 8,
     'rcubands': ['10_90'],
     'rcumodes': [3],
     'sb_range': ['51:461']}

    Set using single frequency, either string or float:
    >>> fb1=FrequencyBand('180e6')
    >>> fb1.__dict__
    {'antsets': ['HBA_JOINED'],
     'beamlets': ['0'],
     'bits': 16,
     'rcubands': ['110_190'],
     'rcumodes': [5],
     'sb_range': ['410']}

    Set using a tuple given by (freqlo, freqhi):
    >>> fbt2=FrequencyBand((40e6,120e6))
    >>> fbt2.__dict__
    {'antsets': ['LBA_INNER', 'HBA_JOINED'],
     'beamlets': ['0:256', '257:308'],
     'bits': 8,
     'rcubands': ['30_90', '110_190'],
     'rcumodes': [4, 5],
     'sb_range': ['205:461', '51:102']}

    Set using a tuple given by (freqlo, freqstep, freqhi):
    >>> fbt3=FrequencyBand((150e6,5.0e6,220e6))
    >>> fbt3.__dict__
    {'antsets': ['HBA_JOINED', 'HBA_JOINED'],
     'beamlets': ['0:8', '9:11'],
     'bits': 8,
     'rcubands': ['110_190', '210_250'],
     'rcumodes': [5, 7],
     'sb_range': ['256,281,306,331,356,381,406,431,456', '51,76,101']}

    Set using slice style index string given by 'freqlo:freqstep:freqhi':
    >>> fbs3=FrequencyBand('35e6:1e6:210e6')
    >>> fbs3.__dict__


     Set using slice style index string given by 'freqlo:freqhi':
     >>> fbs2=FrequencyBand('80e6:10.0e6:220e6')
     >>> fbs2.__dict__
    {'antsets': ['LBA_INNER', 'HBA_JOINED', 'HBA_JOINED'],
     'beamlets': ['0:1', '2:10', '11:12'],
     'bits': 16,
     'rcubands': ['30_90', '110_190', '210_250'],
     'rcumodes': [4, 5, 7],
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
    subarr_rcusel = [{0: '0:191'},
                     {0: '0:47,96:143',
                      1: '48:95,144:191'},
                     {0: '0:31,96:127',
                      1: '32:63,128:159',
                      2: '64:95,160:191'}]

    combos = [[4],[3],[5],[7],[6],[4,5],[3,5],[5,7],[4,5,7],[3,5,7]]

    nrffts = 1024
    nrbeamletsbybits = {16: 244, 8: 488, 4: 976}

    def __init__(self, arg):
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
            return None, None, None, None, None, None
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
            raise ValueError("Requires too many beamlets.")
        for bits in self.nrbeamletsbybits.keys():
            # Choose largest bit size that fulfills nrbeamlets
            if nrbeamlets <= self.nrbeamletsbybits[bits]:
                break
        rcusel = [self.subarr_rcusel[nrrcumodes-1][idx] for idx in range(nrrcumodes)]
        return rcumodes, rcubands, antsets, sb_range, bits, beamlets, rcusel

    def getsampinfo(self,rcumode):
        sampfreq = self.rcumode_smpfrqs[rcumode]
        sampres = sampfreq/self.nrffts
        return sampfreq, sampres

    def freq2sb(self, freq, sampfreq):
        # Using round() but could have other truncation strategy...
        abs_sb = int(round(freq / sampfreq * self.nrffts))
        sb = abs_sb % (self.nrffts/2)
        return sb

    def rcumode_sb2freq(self, rcumode, sb):
        abs_sb = self.rcumode_nqzones[rcumode]*self.nrffts/2+sb
        return abs_sb*self.rcumode_smpfrqs[rcumode]/float(self.nrffts)

    def freqstep2sbstep(self,freqstep, sbwidth):
        return int(math.floor(freqstep/sbwidth))

    def supportrcumodes(self,freqlo,freqhi):
        for combo in self.combos:
            if self.rcumode_passbands[combo[0]][0] <= freqlo \
               and freqhi <= self.rcumode_passbands[combo[-1]][1]:
                return combo
        return None

    def getrcubandnames(self):
        """Return list of rcu band names."""
        return list(self.rcumode_passbands.values())


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


def stdPointings(directionterm='?'):
    """Find beamctl direction string based on direction term.

    Parameters
    ----------
    directionterm : str, '?', '', or None
        Source name or term for a direction. E.g. 'Z' for Zenith
        or 'CasA' for Cassiopeia A.

    Returns
    -------
    beamctldir : str
        Argument suitable for beamctl direction arguments --anadir and
        --digdir. If directionterm is None it will return all the direction terms it
        knows. If it is '', then it will return an empty direction ',,'.

    Raises
    ------
    KeyError
        Thrown if directionterm is not understood.
    """
    term2beamstr = {  # 1e-6 rad < 1arcsec
        '':     ',,',
        'N':    str(0*math.pi/2)+',0.,AZELGEO',
        'E':    str(1*math.pi/2)+',0.,AZELGEO',
        'S':    str(2*math.pi/2)+',0.,AZELGEO',
        'W':    str(3*math.pi/2)+',0.,AZELGEO',
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
    if directionterm is '?':
        return term2beamstr.keys()
    if directionterm in term2beamstr:
        return term2beamstr[directionterm]
    elif directionterm is None:
        return None
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
