#!/usr/bin/env python
"""The module provides methods to export LOFAR antennafield file metadata to
other formats such as CASA .cfg files.
"""
import os
import datetime
import argparse
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from ilisa.antennameta.antennafieldlib import parseAntennaField, getArrayBandParams, list_stations

BANDARRS = ['LBA', 'HBA']
CASA_CFG_DTYPE = [('X',float),('Y',float),('Z',float),('Diam',float),('Name','S5')]
CASA_CFG_FMT = '%12f %12f %12f %4.1f %s'
CASA_CFG_DEST = os.path.join(os.path.dirname(__file__),'share/simmos/')
ALIGNMENT_DEST = os.path.join(os.path.dirname(__file__),'share/alignment/')


def _get_casacfg_header(tier, bandarr=None, stnid=None):
    header = "observatory=LOFAR\n"
    columnlabels = "X Y Z Diam Name"
    if tier == 'ILT':
        header += ("coordsys=XYZ"+"\n"
                  )
    elif tier == 'station':
        header += ("station="+stnid+"\n"
                  +"arrayband="+bandarr+"\n"
                  +"coordsys=XYZ"+"\n"
                  )
    elif tier == 'tile' and bandarr == 'HBA':
        header += ("station="+stnid+"\n"
                  +"arrayband="+bandarr+"\n"
                  +"coordsys=XYZ"+"\n"
                  )
    elif tier == 'rot':
        header += ("station="+stnid+"\n"
                  +"arrayband="+bandarr+"\n"
                  +"coordsys=XYZ"+"\n"
                  +"Rotation matrix (alignment station frame w.r.t. ITRF)"+"\n"
                  )
        columnlabels = "x_hat y_hat z_hat"
    else:
        print("Tier {} not valid.".format(tier))
        raise
    header += ("\n"
              #+"Created by {}\n".format(os.path.basename(__file__))
              +"Created with {}\n".format("LoStationMeta")
              +"Created at {}\n".format(datetime.datetime.utcnow())
              +"\n"
              +columnlabels
              )
    return header


def save_arrcfg_station(bandarr, stnid):
    """Output a station array configuration in a CASA simmos .cfg format.
    """
    stnPos, stnRot, stnRelPos, stnIntilePos = getArrayBandParams(stnid, bandarr)
    header = _get_casacfg_header('station', bandarr, stnid)
    nrelems = stnRelPos.shape[0]
    outtable = np.zeros(nrelems, dtype=CASA_CFG_DTYPE)
    if bandarr == 'LBA':
        diam = 1.5
    elif bandarr == 'HBA':
        diam = 3.0
    outtable['X'] = np.squeeze(stnRelPos[:,0])
    outtable['Y'] = np.squeeze(stnRelPos[:,1])
    outtable['Z'] = np.squeeze(stnRelPos[:,2])
    outtable['Diam'] = diam
    outtable['Name'] = ['ANT'+str(elem) for elem in range(nrelems)]
    filename = os.path.join(CASA_CFG_DEST, stnid+"_"+bandarr+'.cfg')
    np.savetxt(filename, outtable, fmt=CASA_CFG_FMT, header=header)


def save_arrcfg_tile(stnid):
    """Output a station's HBA tile array configuration in CASA simobs .cfg format.
    """
    bandarr = 'HBA'
    hbadeltas = np.asarray(parseiHBADeltasfile(stnid))
    header = _get_casacfg_header('tile', bandarr, stnid)
    nrelems = len(hbadeltas)
    diam = 0.5
    outtable = np.zeros(nrelems, dtype=CASA_CFG_DTYPE)
    outtable['X'] = hbadeltas[:,0]
    outtable['Y'] = hbadeltas[:,1]
    outtable['Z'] = hbadeltas[:,2]
    outtable['Diam'] = diam
    outtable['Name'] = ['ELM'+str(elem) for elem in range(nrelems)]
    filename = os.path.join(CASA_CFG_DEST, stnid+'_tiles.cfg')
    np.savetxt(filename, outtable, fmt=CASA_CFG_FMT, header=header)


def save_all_rotmats():
    """Save all LOFAR station rotation matrices. The rotation matrices define
    the alignment of the stations."""
    stnid_list = sorted(list_stations())
    for stnid in stnid_list:
        for bandarr in BANDARRS:
            stnpos, stnrot, stnrelpos, stnintilepos = \
                         getArrayBandParams(stnid, bandarr)
            print("{} {}".format(stnid, bandarr))
            filename = os.path.join(ALIGNMENT_DEST, '{}_{}.txt'.format(stnid,
                                                                       bandarr))
            header = _get_casacfg_header('rot', bandarr, stnid)
            np.savetxt(filename, stnrot, fmt="%12f %12f %12f", header=header)


def save_all_arrcfg_ILT():
    """Save ILT's array configuration in a CASA sim .cfg format, for both LBA
       and HBA.
    """
    stnId_list = list_stations()
    stnId_list.sort(key=lambda stnid: stnid[2:])
    nrstns = len(stnId_list)
    the_maxbaselines = maxbaselinelen()
    for bandarr in BANDARRS:
        header = _get_casacfg_header('ILT', bandarr)
#        if bandarr == 'LBA':
#            diam = 55.5
#        else:
#            diam = 63.3
        outtable = np.zeros(nrstns, dtype=CASA_CFG_DTYPE)
        for stnnr, stnid in enumerate(stnId_list):
            AntFld = parseAntennaField(stnid)
            position = AntFld[bandarr]['POSITION']
            outtable['X'][stnnr], outtable['Y'][stnnr], outtable['Z'][stnnr] = position
            try:
                outtable['Diam'][stnnr] = the_maxbaselines[stnid][bandarr]
            except KeyError:
                outtable['Diam'][stnnr] = the_maxbaselines[stnid]['HBA0']
            outtable['Name'][stnnr] = stnid
        filename = os.path.join(CASA_CFG_DEST, "ILT_"+bandarr+'.cfg')
        np.savetxt(filename, outtable, fmt=CASA_CFG_FMT, header=header)


def save_all_stations_arrcfg():
    """Save all station array configurations.
    """
    stnId_list = list_stations()
    for bandarr in BANDARRS:
        for stnid in stnId_list:
            print stnid, bandarr
            save_arrcfg_station(bandarr, stnid)


def save_all_tile_arrcfg():
    """Save all tile array configurations.
    """
    stnId_list = list_stations()
    bandarr = 'HBA'
    for stnid in stnId_list:
        save_arrcfg_tile(stnid)


def maxbaselinelen():
    """Compute max baseline length of each station."""
    stnId_list = list_stations()
    the_maxbaselines = {}
    for stnnr, stnid in enumerate(stnId_list):
        the_maxbaselines[stnid] = {}
        AntFld = parseAntennaField(stnid)
        if AntFld['HBA0']['POSITION']:
            nrelems = len(AntFld['HBA']['REL_POS'])
            AntFld['HBA0']['REL_POS'] = AntFld['HBA']['REL_POS'][:nrelems/2]
            AntFld['HBA1']['REL_POS'] = AntFld['HBA']['REL_POS'][nrelems/2:]
            del AntFld['HBA']
        else:
            del AntFld['HBA0']
            del AntFld['HBA1']
        for bandarr in AntFld.keys():
            stnRelPos = np.array(AntFld[bandarr]['REL_POS'])
            nrelems = stnRelPos.shape[0]
            absuvwelemnr = []
            for elemnr in range(nrelems):
                absuvwelemnr.append(np.amax(np.linalg.norm(stnRelPos - stnRelPos[elemnr,:], ord=2, axis=1)))
            #print stnid, bandarr, nrelems, max(absuvwelemnr)
            the_maxbaselines[stnid][bandarr] = max(absuvwelemnr)
    return the_maxbaselines


def plot_arrayconfiguration(stnid, bandarr, coordsys):
    """Plot array configurations"""
    tier = 'station'
    if stnid == 'ILT':
        tier = 'ILT'
        stnId_list = list_stations()
        pos = []
        names = []
        for stnid in stnId_list:
            names.append(stnid)
            stnpos, stnrot, stnrelpos, stnintilepos = \
                    getArrayBandParams(stnid, bandarr)
            pos.append(np.asarray(stnpos).squeeze().tolist())
        pos = np.array(pos)
    else:
        if bandarr == 'tile':
            tier = 'tile'
            bandarr = 'HBA'
        stnpos, stnrot, stnrelpos, stnintilepos = getArrayBandParams(\
                                                                 stnid, bandarr)
        if coordsys == 'local':
            stnrelpos = stnrelpos * stnrot
            stnintilepos = stnintilepos * stnrot
        if tier == 'tile':
            pos = np.asarray(stnintilepos)
            nameprefix = 'elem'
        else:
            pos = np.asarray(stnrelpos)
            nameprefix = 'ant'
        names = [nameprefix+str(elem) for elem in range(pos.shape[0])]
    # Plot
    fig = plt.figure()
    ax = fig.gca(projection='3d')
    ax.plot(pos[:,0],pos[:,1],pos[:,2],'*')
    for idx, name in enumerate(names):
        ax.text(pos[idx,0],pos[idx,1],pos[idx,2],name)
    ax.set_aspect('equal')
    plt.show()


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("stnid")
    parser.add_argument("bandarr")
    parser.add_argument("coordsys", nargs='?', default='abs')
    args = parser.parse_args()
    plot_arrayconfiguration(args.stnid, args.bandarr, args.coordsys)


if __name__ == '__main__':
    save_all_arrcfg_ILT()
    #save_all_stations_arrcfg()
    #save_all_tile_arrcfg()
    #save_all_rotmats()
    #main()
#    m= maxbaselinelen()
#    for k in m.keys():
#        for b in m[k].keys():
#            print k, b, m[k][b]
