import datetime

import numpy as np
from casacore.measures import measures
from casacore.quanta import quantity


def ITRF2lonlat(x_itrf, y_itrf, z_itrf):
    """\
    Convert ITRF cartesian position coordinates to WGS84 latitude and longitude

    Parameters
    ----------
    x_itrf: float
        X coordinate of position in ITRF system in meters.
    y_itrf: float
        Y coordinate of position in ITRF system in meters.
    z_itrf: float
        Z coordinate of position in ITRF system in meters.

    Returns
    -------
    lon: float
        Longitude in degrees
    lat: float
        Latitude in degrees
    hgt: float
        Height above surface in meters.

    Examples
    --------
    >>> from ilisa.calim.geodesy import ITRF2lonlat
    >>> ITRF2lonlat(3370286.88256, 712053.913283, 5349991.484)
    11.929671631184405, 57.39876274671682, 41.63424105290324)
    """
    dm = measures()
    posstr = ["{}m".format(crd) for crd in (x_itrf, y_itrf, z_itrf)]
    p_itrf = dm.position('itrf', *posstr)
    p_wgs84 = dm.measure(p_itrf, 'wgs84')
    _ref = dm.get_ref(p_wgs84)
    lon = quantity(p_wgs84['m0']).get('deg').get_value()
    lat = quantity(p_wgs84['m1']).get('deg').get_value()
    hgt = quantity(p_wgs84['m2']).get('m').get_value()
    return lon, lat, hgt


def lonlat2ITRF(lon, lat, h='0.0m'):
    """\
    Convert a longtitude and latitude specified position to cartesian ITRF


    Parameters
    ----------
    lon : str
        Longitude with units string.
    lat : str
        Latitude with units string.
    h : str
        Height above geodetic with unit string.

    Returns
    -------
    x, y, z : float
        Cartesian ITRF position.

    Examples
    --------
    >>> from ilisa.calim.geodesy import lonlat2ITRF
    >>> lonlat2ITRF('11.929671631184405deg', '57.39876274671682deg', '41.634m')
    3370286.882433161, 712053.913256202, 5349991.4837967735
    """
    dm = measures()
    p_wgs84 = dm.position('wgs84', lon, lat, h)
    p_itrf = dm.measure(p_wgs84, 'itrf')
    ilon = quantity(p_itrf['m0']).get('m').get_value()
    ilat = quantity(p_itrf['m1']).get('m').get_value()
    irad = quantity(p_itrf['m2']).get('m').get_value()
    x = irad * np.cos(ilon) * np.cos(ilat)
    y = irad * np.sin(ilon) * np.cos(ilat)
    z = irad * np.sin(ilat)
    return x, y, z


def utcpos2lmst(utctim, lonlat):
    """\
    Convert UTC to local mean sidereal time

    Parameters
    ----------
    utctim : datetime
        A UTC date-time of an epoch.
    lonlat : tuple
        A tuple of len 3 (lon, lat, hgt), or 2 (lon, lat) where hgt=0.

    Returns
    -------
    lsmt_timdel : tuple
    """
    # Convert utctim instant to ISO datetime str
    if type(utctim) is datetime.datetime:
        dattim_str = utctim.isoformat()
    else:
        dattim_str = utctim
    me = measures()
    if len(lonlat) < 3:
        lonlat += ('0m',)
    lonlathgt = []
    for _e in lonlat[:2]:
        if type(_e) is float:
            _e = str(_e) + 'deg'
        lonlathgt.append(_e)
    if type(lonlat[2]) is float or type(lonlat[2]) is int:
        lonlathgt.append(str(lonlat[2])+'m')
    pos = me.position('WGS84', *lonlathgt)
    me.doframe(pos)
    epoch = me.epoch('UTC', quantity(dattim_str))
    lmst = me.measure(epoch, 'LMST')
    fracday = lmst['m0']['value'] % 1
    # hms_str = quantity(str(fracday)+'d').to_time().formatted()
    lsmt_timdel = datetime.timedelta(days=fracday)
    return lsmt_timdel
