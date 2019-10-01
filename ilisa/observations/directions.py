import os
import math

import numpy
import yaml

import ilisa.observations


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


def pointing_str2tuple(beamctldirarg):
    """Convert a beamctl direction string into direction tuple.

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
        'CasA': '6.123487,1.026515,J2000', #3C-461
        'CygA': '5.233660,0.710940,J2000', #3C-405
        'TauA': '1.459672,0.384225,J2000', #3C-144
        'VirA': '3.276086,0.216265,J2000', #3C-274
        'Sun':  '0.,0.,SUN',
        'Jupiter': '0.,0.,JUPITER',
        'Moon': '0.,0.,MOON',
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
    beamctldir = pointing_str2tuple(gendirstr)
    if beamctldir is None:
        try:
            beamctldirstr = stdPointings(gendirstr)
        except:
            raise ValueError('General direction term {} unknown.'.format(gendirstr))
    else:
        beamctldirstr = gendirstr
    return beamctldirstr


def lookupsource(src_name):
    """Lookup the pointing direction for the name of a source.

    Parameters
    ----------
    src_name : str
        Name of source.

    Returns
    -------
    pointing : tuple
        Length 3 tuple with (az, el, ref).

    """
    user_srcs_dir = os.path.join(ilisa.observations.user_data_dir, 'sources')
    user_srcs_files = os.listdir(user_srcs_dir)

    with open(os.path.join(user_srcs_dir, user_srcs_files.pop())) as usf:
        user_srcs = yaml.load(usf)
    try:
        pointing = user_srcs[src_name]
    except KeyError:
        raise RuntimeError('Source {} not found in {}.'.format(src_name,
                                                               user_srcs_files))
    return pointing