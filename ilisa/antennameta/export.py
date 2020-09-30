#!/usr/bin/env python
"""The module provides methods to export LOFAR antennafield file metadata to
other formats such as CASA .cfg files.
"""
import os
import datetime
import argparse
import numpy as np
from ilisa.antennameta.antennafieldlib import parseAntennaField, parseiHBADeltasfile, \
    getArrayBandParams, list_stations, BANDARRS

CASA_CFG_DTYPE = [('X', float), ('Y', float), ('Z', float), ('Diam', float),
                  ('Name', 'S7')]
CASA_CFG_FMT = '%12f %12f %12f %4.1f %s'
CASA_CFG_DEST = os.path.join(os.path.dirname(__file__), 'share/simmos/')
ALIGNMENT_DEST = os.path.join(os.path.dirname(__file__), 'share/alignment/')


def _get_casacfg_header(tier, bandarr=None, stnid=None):
    """Make a casa config header."""
    header = "observatory=LOFAR\n"
    columnlabels = "X Y Z Diam Name"
    if tier == 'ILT':
        header += "coordsys=XYZ" + "\n"
    elif tier == 'station':
        header += ("station="+stnid+"\n"
                   + "arrayband="+bandarr+"\n"
                   + "coordsys=XYZ"+"\n")
    elif tier == 'tile' and bandarr == 'HBA':
        header += ("station="+stnid+"\n"
                   + "arrayband="+bandarr+"\n"
                   + "coordsys=XYZ"+"\n")
    elif tier == 'rot':
        header += ("station="+stnid+"\n"
                   + "arrayband="+bandarr+"\n"
                   + "coordsys=XYZ"+"\n"
                   + "Rotation matrix (alignment station frame w.r.t. ITRF)"+"\n")
        columnlabels = "x_hat y_hat z_hat"
    else:
        raise ValueError("Tier {} not valid.".format(tier))
    header += ("\n"
               + "Created with {}\n".format("iLiSA")
               + "Created at {}\n".format(datetime.datetime.utcnow())
               + "\n"
               + columnlabels)
    return header


def output_arrcfg_station(stnid, bandarr, output='default'):
    """Output a station array (given by stnid and bandarr) configuration in
    a CASA simmos .cfg format.
    """
    _pos, stn_rot, stn_relpos, _intile_pos = getArrayBandParams(stnid, bandarr)
    header = _get_casacfg_header('station', bandarr, stnid)
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


def output_arrcfg_tile(stnid):
    """Output a station's HBA tile array configuration in
    CASA simobs .cfg format.
    """
    bandarr = 'HBA'
    hbadeltas = np.asarray(parseiHBADeltasfile(stnid))
    header = _get_casacfg_header('tile', bandarr, stnid)
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
    :param output: Name of output file. If set to 'default', will use default name
                   convention, i.e. '<stnid>_<bandarr>.txt'.
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
        if ant_fld['HBA0']['POSITION']:
            nrelems = len(ant_fld['HBA']['REL_POS'])
            ant_fld['HBA0']['REL_POS'] = ant_fld['HBA']['REL_POS'][:nrelems/2]
            ant_fld['HBA1']['REL_POS'] = ant_fld['HBA']['REL_POS'][nrelems/2:]
            del ant_fld['HBA']
        else:
            del ant_fld['HBA0']
            del ant_fld['HBA1']
        for bandarr in ant_fld.keys():
            stn_rel_pos = np.array(ant_fld[bandarr]['REL_POS'])
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


def main():
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
    args = parser.parse_args()

    if args.stationid is not None:
        stn_id_list = [args.stationid]
    bandarrs = BANDARRS
    if args.bandarr is not None:
        bandarrs = [args.bandarr]

    the_maxbaselines = max_stn_baselines()
    for bandarr in bandarrs:
        header = _get_casacfg_header('ILT', bandarr)
        #        if bandarr == 'LBA':
        #            diam = 55.5
        #        else:
        #            diam = 63.3
        outtable = np.zeros(0, dtype=CASA_CFG_DTYPE)
        for stnnr, stnid in enumerate(stn_id_list):
            if not (stnid == 'NenuFAR' and bandarr == 'HBA'):
                print('Doing {} {}'.format(stnid, bandarr))
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
                output_arrcfg_station(stnid, bandarr)
                output_rotmat_station(stnid, bandarr)
            if bandarr == 'HBA' and stnid != 'NenuFAR':
                output_arrcfg_tile(stnid)
        # output array cfg_for ILT for this bandarr:
        filename = os.path.join(CASA_CFG_DEST, "ILT_" + bandarr + '.cfg')
        np.savetxt(filename, outtable, fmt=CASA_CFG_FMT, header=header)


if __name__ == '__main__':
    main()
