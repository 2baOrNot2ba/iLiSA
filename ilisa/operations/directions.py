import os
import math
import numpy
import re
import yaml

import ilisa.operations


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
    """
    Convert a beamctl direction string into direction tuple

    Inverse of pointing_tuple2str().

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
        beamctldirargs = beamctldirarg.split(',')
        angle1str, angle2str, refsys = [arg.strip() for arg in beamctldirargs]
        angle1 = gen_angle2rad(angle1str)
        angle2 = gen_angle2rad(angle2str)
    except:
        return None
    # TODO should check that refsys is one of the valid refsys strings.
    dirtuple = (angle1, angle2, refsys)
    return dirtuple


def pointing_tuple2str(dirtuple):
    """
    Convert a direction tuple to a string
    
    Inverse of pointing_str2tuple().

    Parameters
    ----------
    dirtuple : tuple or None
        Direction tuple defined as (angle1: float, angle2: float, refsys: str)
        or None if beamctldirarg was not correct format.

    Returns
    -------
    beamctldirarg : str
        String with format 'angle1,angle2,refsys'. Returns None if input is incorrect value.
    
    Examples
    --------
    >>> from ilisa.operations.directions import pointing_tuple2str,
    ...                                           pointing_str2tuple
    >>> pointing_tuple2str((1.0,2.0,'AZELGEO'))
    '1.0,2.0,AZELGEO'
    >>> pointing_str2tuple(pointing_tuple2str((1.0,2.0,'AZELGEO')))
    (1.0, 2.0, 'AZELGEO')
    """
    beamctldirarg = None
    if type(dirtuple) is tuple and len(dirtuple) == 3:
        dir_str_tuple = (str(dirtuple[0]), str(dirtuple[1]), dirtuple[2])
        beamctldirarg = ",".join(dir_str_tuple)
    return beamctldirarg


def check_directionstr(beamctldirarg):
    """\
    Check direction string is correctly formatted

    Parameters
    ----------
    beamctldirarg : str
        Direction argument to beamctl process. Format e.g. 'az,el,J2000'.
        A value of None, is interpreted as valid and means no direction.

    Returns
    -------
    bool
        True if beamctldirarg is a valid beamctl direction argument, else False.
    """
    if beamctldirarg is None:
        return True
    if pointing_tuple2str(pointing_str2tuple(beamctldirarg)):
        return True
    else:
        return False


def std_pointings(directionterm='?'):
    """
    Find beamctl direction string based on direction term.

    Parameters
    ----------
    directionterm : str, '?', '', or None
        Source name or term for a direction. E.g. 'Z' for Zenith
        or 'Cas-A' for Cassiopeia A.

    Returns
    -------
    beamctldir : str
        Argument suitable for beamctl direction arguments --anadir and
        --digdir. If directionterm is None it will return all the direction
        terms it knows. If it is '', then it will return an empty
        direction ',,'.

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
        'NE':   str(0.5 * math.pi / 2) + ',0.,AZELGEO',
        'SE':   str(1.5 * math.pi / 2) + ',0.,AZELGEO',
        'SW':   str(2.5 * math.pi / 2) + ',0.,AZELGEO',
        'NW':   str(3.5 * math.pi / 2) + ',0.,AZELGEO',
        'CasA':  '6.123487,1.026515,J2000',  # Alternative to 'Cas-A'
        'Cas-A': '6.123487,1.026515,J2000',  # 3C-461
        'Cyg-A': '5.233660,0.710940,J2000',  # 3C-405
        'Tau-A': '1.459672,0.384225,J2000',  # 3C-144
        'Vir-A': '3.276086,0.216265,J2000',  # 3C-274
        'Sun':  '0.,0.,SUN',
        'Moon': '0.,0.,MOON',
        'Jupiter': '0.,0.,JUPITER',
        'Uranus': '0.,0.,URANUS',
        'NCP':  '0.,'+str(math.pi/2)+',ITRF'
    }
    if directionterm == '?':
        return term2beamstr.keys()
    if directionterm in term2beamstr:
        return term2beamstr[directionterm]
    elif directionterm is None:
        return None
    else:
        pointing = lookupsource(directionterm)
        if not pointing:
            raise KeyError('Requested source {} unknown.'.format(directionterm))
        return pointing


def directionterm2tuple(directionterm):
    """
    Get direction tuple from (str) term

    Parameters
    ----------
    directionterm : str
        Named direction
    
    Returns
    -------
    directiontuple : (float, float, str)
        Direction given as angle1, angle2 and reference system respectively.
    """
    return pointing_str2tuple(std_pointings(directionterm))


def normalizebeamctldir(gendirstr):
    """
    Parse a general direction string.

    Parameters
    ----------
    gendirstr : str
        This could be one of the direction terms or it could a beamctl
        direction string.

    Returns
    -------
    beamctldirstr : str
        The beamctl direction string corresponding to input gendirstr.
        If gendirstr is
        None, then beamctldirstr is None.
    """
    if gendirstr == 'None' or gendirstr is None:
        return None
    beamctldir = pointing_str2tuple(gendirstr)
    if beamctldir is None:
        try:
            beamctldirstr = std_pointings(gendirstr)
        except:
            raise ValueError(
                        'General direction term {} unknown.'.format(gendirstr))
    else:
        beamctldirstr = pointing_tuple2str(beamctldir)
    return beamctldirstr


def lookupsource(src_name):
    """
    Lookup the pointing direction for the name of a source stored in user
    created files.

    Parameters
    ----------
    src_name : str
        Name of source.

    Returns
    -------
    direction : tuple
        Length 3 tuple with (az, el, ref).

    """
    src_database = {}
    user_srcs_dir = os.path.join(ilisa.operations.USER_DATA_DIR,
                                 'source_catalogs')
    user_srcs_files = os.listdir(user_srcs_dir)

    for user_srcs_file in user_srcs_files:
        with open(os.path.join(user_srcs_dir, user_srcs_file)) as usf:
            srcs_cat = yaml.safe_load(usf)
        srcs_data = srcs_cat['sources']
        for source_entry in srcs_data:
            source_names = source_entry['name']
            _direction_str = source_entry['direction']
            if type(source_names) is not list:
                source_names = [source_names]
            for _source_name in source_names:
                # No conversion:
                src_database[_source_name] = _direction_str
    # Lookup src_name in src_database else return None
    direction = src_database.get(src_name)
    return direction


def _req_calsrc_proc(req_calsrc, allsky, pointingstr):
    # TODO absorb this function into rest of module
    # Helper function Determine phaseref 
    if req_calsrc is not None:
        pntstr = std_pointings(req_calsrc)
    elif allsky:
        pntstr = std_pointings('Z')
    else:
        pntstr = pointingstr
    phaseref = tuple(pntstr.split(','))
    return phaseref


def read_sched_srccatre(filename):
    """Read a sched keyin free formatted source catalogue."""
    with open(filename) as f:
        srccat_raw = f.read()
    srccat_bl = re.sub('!.*\n', '\n', srccat_raw)  # Remove comments on lines
    srccat = re.sub('^\s*\n', '', srccat_bl)  # Remove blank lines
    srclistentries_split = re.split('/', srccat, re.MULTILINE)
    srclistentries_el = [line.replace('\n', ' ') for line in
                         srclistentries_split]
    srclistentries = [line for line in srclistentries_el if line != '']
    keywords = ['source', 'equinox', 'flux', 'fluxref', 'ra', 'raerr', 'dec',
                'decerr',
                'calcode', 'remarks']
    keywordsre = '(' + ')('.join(keywords) + ')'
    valre = '[^' + keywordsre + '($)]+'
    equinox_default = None
    for srclistentry in srclistentries:
        src_ent = {}
        equinox = None
        for kw in keywords:
            kifreefrmt = r"{0}\s*[\s=]\s*({1})".format(kw, valre)
            m=re.search(kifreefrmt, srclistentry, re.MULTILINE | re.IGNORECASE)
            if m:
                val = m.group(1).rstrip()
                if kw == 'equinox':
                    equinox =  val
                if kw == 'source':
                    # TODO Be a bit more careful here
                    val = val.replace("'","")
                    val = val.split(',')
                src_ent[kw] = val
        if equinox:
            equinox_default = equinox
        else:
            equinox = equinox_default
        src_ent['equinox'] = equinox
        print(src_ent)


def gen_angle2rad(anglestr):
    """\
    Parse angle string to radians

    Parameters
    ----------
    anglestr : str
        Angle string in format accepted by casacore.

    Returns
    -------
    angle : float
        Angle corresponding to input in radians.
    """
    from casacore.quanta import quantity

    to_unit = 'rad'
    ang = quantity(anglestr)
    return ang.get_value(to_unit)


import sys


if __name__=="__main__":
    print(gen_angle2rad(sys.argv[1]), 'rad')
