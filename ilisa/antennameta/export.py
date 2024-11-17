#!/usr/bin/env python
"""The module provides methods to export LOFAR antennafield file metadata to
other formats such as CASA .cfg files.
"""
import os
import datetime
import argparse
import numpy as np
from casacore.measures import measures
from casacore.quanta import quantity
from ilisa.antennameta.antennafieldlib import parseAntennaField, \
     parseiHBADeltasfile, getArrayBandParams, list_stations, BANDARRS

CASA_CFG_DTYPE = [('X', float), ('Y', float), ('Z', float), ('Diam', float),
                  ('Name', 'U7')]
CASA_CFG_FMT = '%12f %12f %12f %4.1f %s'
CASA_CFG_DEST = os.path.join(os.path.dirname(__file__), 'share/simmos/')
ALIGNMENT_DEST = os.path.join(os.path.dirname(__file__), 'share/alignment/')


def _get_casacfg_header(tier, bandarr=None, stnid=None, coordsys=None):
    """Make a casa config header."""
    header = "observatory=LOFAR\n"
    columnlabels = "X Y Z Diam Name"
    if tier == 'rot':
        header += ("station="+stnid+"\n"
                   + "arrayband="+bandarr+"\n"
                   + "coordsys=XYZ"+"\n"
                   + "Rotation matrix (alignment station frame w.r.t. ITRF)"+"\n")
        columnlabels = "x_hat y_hat z_hat"
    else:
        if tier == 'station':
            header += ("station=" + stnid + "\n"
                       + "arrayband=" + bandarr + "\n")
        elif tier == 'tile' and bandarr == 'HBA':
            header += ("station=" + stnid + "\n"
                       + "arrayband=" + bandarr + "\n")
        if coordsys == 'ITRF':
            header += "coordsys=XYZ" + "\n"
        elif coordsys.lower().startswith('loc'):
            header += "coordsys=LOC (local tangent plane)" + "\n"
    header += ("\n"
               + "Created with {}\n".format("iLiSA")
               + "Created at {}\n".format(datetime.datetime.utcnow())
               + "\n"
               + columnlabels)
    return header


def output_arrcfg_station(stnid, bandarr, coordsys, output='default'):
    """Output a station array (given by stnid and bandarr) configuration in
    a CASA simmos .cfg format.
    """
    _pos, stn_rot, stn_relpos, _intile_pos = getArrayBandParams(stnid, bandarr)
    header = _get_casacfg_header('station', bandarr, stnid, coordsys)
    nrelems = stn_relpos.shape[0]
    outtable = np.zeros(nrelems, dtype=CASA_CFG_DTYPE)
    if bandarr == 'LBA':
        diam = 1.5
    else:  # bandarr == 'HBA'
        diam = 3.0
    outtable['X'] = np.squeeze(stn_relpos[:, 0])
    outtable['Y'] = np.squeeze(stn_relpos[:, 1])
    outtable['Z'] = np.squeeze(stn_relpos[:, 2])
    outtable['Diam'] = diam
    outtable['Name'] = ['ANT'+str(elem) for elem in range(nrelems)]
    if output == 'default':
        output = os.path.join(CASA_CFG_DEST, stnid+"_"+bandarr+'.cfg')
    np.savetxt(output, outtable, fmt=CASA_CFG_FMT, header=header)


def output_arrcfg_tile(stnid, coordsys):
    """Output a station's HBA tile array configuration in
    CASA simobs .cfg format.
    """
    bandarr = 'HBA'
    hbadeltas = np.asarray(parseiHBADeltasfile(stnid))
    header = _get_casacfg_header('tile', bandarr, stnid, coordsys)
    nrelems = len(hbadeltas)
    diam = 0.5
    outtable = np.zeros(nrelems, dtype=CASA_CFG_DTYPE)
    outtable['X'] = hbadeltas[:, 0]
    outtable['Y'] = hbadeltas[:, 1]
    outtable['Z'] = hbadeltas[:, 2]
    outtable['Diam'] = diam
    outtable['Name'] = ['ELM'+str(elem) for elem in range(nrelems)]
    filename = os.path.join(CASA_CFG_DEST, stnid+'_tiles.cfg')
    np.savetxt(filename, outtable, fmt=CASA_CFG_FMT, header=header)


def output_rotmat_station(stnid, bandarr, output='default'):
    """
    Save a station bandarray's rotation matrix.
    :param stnid: Station ID.
    :param bandarr: 'HBA' or 'LBA'
    :param output: Name of output file. If set to 'default', will use default
                   name convention, i.e. '<stnid>_<bandarr>.txt'.
    """
    stnpos, stnrot, stnrelpos, stnintilepos = getArrayBandParams(stnid,
                                                                 bandarr)
    if output == 'default':
        output = os.path.join(ALIGNMENT_DEST, '{}_{}.txt'.format(stnid,
                                                                 bandarr))
    header = _get_casacfg_header('rot', bandarr, stnid)
    np.savetxt(output, stnrot, fmt="%12f %12f %12f", header=header)


def max_stn_baselines():
    """Compute max baseline length of each station."""
    stn_id_list = list_stations()
    the_maxbaselines = {}
    for stnnr, stnid in enumerate(stn_id_list):
        the_maxbaselines[stnid] = {}
        ant_fld = parseAntennaField(stnid)
        if 'HBA0' in ant_fld:
            nrelems = len(ant_fld['HBA']['REL_POS_X'])
            ant_fld['HBA0']['REL_POS_X'] = ant_fld['HBA']['REL_POS_X'][:nrelems//2]
            ant_fld['HBA1']['REL_POS_X'] = ant_fld['HBA']['REL_POS_X'][nrelems//2:]
            del ant_fld['HBA']
        for bandarr in ant_fld.keys():
            stn_rel_pos = np.array(ant_fld[bandarr]['REL_POS_X'])
            nrelems = stn_rel_pos.shape[0]
            absuvwelemnr = []
            for elemnr in range(nrelems):
                absuvwelemnr.append(np.amax(np.linalg.norm(
                    stn_rel_pos - stn_rel_pos[elemnr, :], ord=2, axis=1)))
            the_maxbaselines[stnid][bandarr] = max(absuvwelemnr)
    return the_maxbaselines


def print_maxbaselines():
    """print maximum baseline for all stations."""
    m = max_stn_baselines()
    for k in m.keys():
        for b in m[k].keys():
            print(k, b, m[k][b])


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
    >>> from ilisa.antennameta.export import ITRF2lonlat
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
    >>> from ilisa.antennameta.export import lonlat2ITRF
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


def cli_export():
    """Export AntennaField data for all or selected stations to casa CSV files.
    """
    stn_id_list = list_stations()
    stn_id_list.sort(key=lambda _stnid: _stnid[2:])

    parser = argparse.ArgumentParser()
    parser.add_argument('-s', '--stationid',
                        help="""Station ID to process.
                        Choose from {}
                        If not given, do all.
                        """.format(stn_id_list))
    parser.add_argument('-b', '--bandarr',
                        help="""Band array to process.
                                Choose from {}.
                                If not given, do all.
                                """.format(BANDARRS))
    parser.add_argument('-c', '--coordsys', type=str, default=None,
                        help="Choose from coord. sys.: LOCAL, ITRF (Def: None)")
    args = parser.parse_args()

    if args.stationid is not None:
        stn_id_list = [args.stationid]
    bandarrs = BANDARRS
    if args.bandarr is not None:
        bandarrs = [args.bandarr]
    coordsys = args.coordsys
    if coordsys is None:
        if args.stationid is None:
            coordsys = 'ITRF'
        else:
            coordsys = 'LOCAL'

    the_maxbaselines = max_stn_baselines()
    for bandarr in bandarrs:
        header = _get_casacfg_header('ILT', bandarr, coordsys='ITRF')
        #        if bandarr == 'LBA':
        #            diam = 55.5
        #        else:
        #            diam = 63.3
        outtable = np.zeros(0, dtype=CASA_CFG_DTYPE)
        for stnnr, stnid in enumerate(stn_id_list):
            if stnid == 'NenuFAR' and bandarr == 'HBA':
                continue
            row = np.zeros(1, dtype=CASA_CFG_DTYPE)
            ant_fld = parseAntennaField(stnid)
            position = ant_fld[bandarr]['POSITION']
            row['X'], row['Y'], row['Z'] = position
            try:
                row['Diam'] = the_maxbaselines[stnid][bandarr]
            except KeyError:
                row['Diam'] = the_maxbaselines[stnid]['HBA0']
            row['Name'] = stnid
            outtable = np.append(outtable, row)
            output_arrcfg_station(stnid, bandarr, coordsys)
            output_rotmat_station(stnid, bandarr)
        # output array cfg_for ILT for this bandarr:
        filename = os.path.join(CASA_CFG_DEST, "ILT_" + bandarr + '.cfg')
        np.savetxt(filename, outtable, fmt=CASA_CFG_FMT, header=header)


if __name__ == '__main__':
    cli_export()
