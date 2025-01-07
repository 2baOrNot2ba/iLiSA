#!/usr/bin/env python
"""A module to read LOFAR antenna field files.
The ROTATION_MATRIX field is the rotation matrix of the station;
it maps the the local station coordinates to the ITRF coordinates
such that

    r_ITRF = ROTATION_MATRIX * r_stn where r are 3-element column vectors.
"""
import argparse
import os
from os.path import dirname
import numpy as np
from matplotlib import pyplot as plt
from ilisa.operations.modeparms import str2elementMap2
from ilisa.calim.geodesy import ITRF2lonlat

__version__ = 0.3

STATICMETADATA = os.path.join(dirname(__file__),'share/StaticMetaData/')
ANTENNAFIELDDIR = STATICMETADATA
IHBADELTASDIR = os.path.join(dirname(__file__),'share/iHBADeltas/')
IHBADELTASDIR = STATICMETADATA

COMMENT_CHAR = '#'
AFfileNameType = 3

# LOFAR band array's elements diameter or 1D extent (approximately) 
ELEMENT_DIAMETER = {'LBA': 2.0,
                    'HBA': 5.0}  # meters
BANDARRS = ['LBA', 'HBA']
HBASUBARRS = ['HBA0', 'HBA1']

def _getAntennaFieldFile(stationName):
    antenna_field_dir = ANTENNAFIELDDIR
    if AFfileNameType==2 or AFfileNameType==3:
       basename = stationName+'-'+'AntennaField'+'.conf'
    else:
       basename='AntennaField'+stationName+'.conf'
    filepath = os.path.join(antenna_field_dir, basename)
    return filepath


def _getiHBADeltafile(stationName):
    """Get file path to iHBADelta for given stationName. """
    basename = stationName+'-'+'iHBADeltas'+'.conf'
    filepath = os.path.join(IHBADELTASDIR, basename)
    return filepath


def parseAntennaField(stationName):
    filepath = _getAntennaFieldFile(stationName)
    return parseAntennaFieldFile(filepath)


def parseAntennaFieldFile(filename, AFfileNameType=3):
    """Parse LOFAR AntennaField file by name and return data as a dict.
    Note that this reads in both LBA and HBA parameters.
    """
    # All bandarrays and hba-sub-bandarrays will have attributes
    _all_ba_min = ['POSITION' , 'NORMAL_VECTOR', 'ROTATION_MATRIX']
    antflddata = {}
    for antband in BANDARRS:
        antflddata[antband] = {k: [] for k in _all_ba_min}
        # All bandarrays but not hba-sub-bandarrays will have attribute
        antflddata[antband]['REL_POS_X'] = []
        antflddata[antband]['REL_POS_Y'] = []
    try:
        f = open(filename)
    except IOError:
        print("Error: "+filename+" does not exist.")
        raise
    line = f.readline()
    while line:
        if COMMENT_CHAR in line:
            line, comment = line.split(COMMENT_CHAR, 1)
        if ("HBA0" in line or "HBA1" in line) and ('HBA0' not in antflddata):
            for antband in HBASUBARRS:
                antflddata[antband] = {k: [] for k in _all_ba_min}
        if "HBA0" in line:
            antband = "HBA0"
        elif "HBA1" in line:
            antband = "HBA1"
        elif "HBA" in line:
            antband = "HBA"
        elif "LBA" in line:
            antband = "LBA"
        else:
            line=f.readline()
            continue
        where, rest = line.split(antband, 1)
        where = where.strip()
        if where == '':
            # Read absolute position of station origin
            line = f.readline()
            elementposshape, elementposLine = line.split(' ',1)
            elementposLine = elementposLine.strip('[] \n').split()
            position = [float(v) for v in elementposLine]
            antflddata[antband]['POSITION'] = position
            if antband != "HBA0" and antband != "HBA1":
                # Read relative position of each element
                line = f.readline()
                dimstr,rest = line.split('[',1)
                shp = dimstr.split('x')
                if shp[0][0] == '(':
                    (idxbeg,idxend) = shp[0].strip('() ').split(',')
                    nrrows = int(idxend)+1
                else:
                    nrrows = int(shp[0])
                for elementNr in range(0,nrrows):
                    line = f.readline()
                    vals = line.split()
                    posxpol = [float(v) for v in vals[0:3]]
                    posypol = [float(v) for v in vals[3:6]]
                    antflddata[antband]['REL_POS_X'].append(posxpol)
                    antflddata[antband]['REL_POS_Y'].append(posypol)
                    # Note: Skip ypol as it is identical to xpol
                # Read ending ']' line
                _ = f.readline()
        elif where == 'NORMAL_VECTOR':
            line = f.readline()
            elementposshape, elementposLine = line.split(' ',1)
            elementposLine = elementposLine.strip('[] \n').split()
            nrmv = [float(v) for v in elementposLine]
            antflddata[antband][where] = nrmv
        elif where == 'ROTATION_MATRIX':
            line = f.readline()
            elementposshape, elementposLine = line.split(' ',1)
            elementposLine = elementposLine.strip('[] \n').split()
            elementposLine = line.split()
            dimstr, rest = line.split('[',1)
            shp = dimstr.split('x')
            for xyz in range(3):
              line = f.readline()
              rowstr = line.split()
              row = [float(v) for v in rowstr]
              antflddata[antband][where].append(row)
            # Read ending ']' line
            _ = f.readline()
        line = f.readline()
    return antflddata


def parseiHBADeltasfile(stationName):
    """Parse iHBADelta file."""
    iHBADeltasdata = []
    filepath = _getiHBADeltafile(stationName)
    f = open(filepath)
    line=f.readline()
    HBADeltaLine = False
    while line:
        if COMMENT_CHAR in line:
            line, comment = line.split(COMMENT_CHAR, 1)
        if "HBADeltas" in line:
            HBADeltaLine = True
            break
        else:
            line=f.readline()
            continue
    if HBADeltaLine:
        line = f.readline()
        nrelems, nrdims = line.strip().rstrip('[').split('x')
        if nrelems[0] == '(':
            nrelems = int(nrelems.strip('() ').split(',')[1])+1
            nrdims = int(nrdims.strip('() ').split(',')[1])+1
        else:
            nrelems, nrdims = int(nrelems), int(nrdims)
        for elemnr in range(nrelems):
            line = f.readline()
            xpos, ypos, zpos = line.lstrip().split()
            iHBADeltasdata.append([float(xpos), float(ypos), float(zpos)])
            #elempos[elemnr,0], elempos[elemnr,1], elempos[elemnr,2] =\
            #                 float(xpos), float(ypos), float(zpos)
    else:
        raise RuntimeError("Error: iHBADeltas file is corrupt.")
    return iHBADeltasdata


def getHBAsepton(stnName, HBAsingleelems):
    """Get configuration parameters for an HBA station in SEPTON mode. SEPTON
    stands for single element per tile ON. In this mode the HBA can be seen as
    an instantenously omni-directional interferometer capable of producing
    allsky snapshots.
    """
    stnPos, stnRot, stnRelPos, stnIntilePos = getArrayBandParams(stnName, 'HBA')
    nrtiles = len(HBAsingleelems)
    relpospertile = np.zeros((nrtiles,3))
    for tilenr in range(nrtiles):
        relpospertile[tilenr] = stnIntilePos[HBAsingleelems[tilenr]]
    stnRelPos += relpospertile
    return stnPos, stnRot, stnRelPos


def getArrayBandParams(stnid, arrband):
    """Get configuration parameters for an array band of the station.
    Array band can be HBA or LBA.
    
    Note: Core stations (CS) antenna field files specify rotation and normal
          vector for HBA0 & HBA1 subarrays while their relative positions
          are given specified as HBA.
    Note: Normal vector can be obtain from last column in rotation matrix.
    """
    antfld = parseAntennaField(stnid)
    stnLoc = stnid[0:2]
    subarr = None
    if arrband == 'LBA':
        subarr = 'LBA'
        hbadeltas = [0.,0.,0.]
    elif arrband.startswith('HBA'):
        # Set subarr appropriately
        if stnLoc == 'CS':
            if arrband == 'HBA0' or arrband == 'HBA1':
                subarr = arrband
            elif arrband == 'HBA':
                # Arbitrarily choose HBA0 (HBA1 could also be used)
                subarr = 'HBA0' 
        elif arrband == 'HBA':
            subarr = arrband
        arrband = 'HBA'
        hbadeltas = parseiHBADeltasfile(stnid)
    else:
        raise ValueError("Array band '{}' invalid. Only 'LBA', 'HBA' are valid,"
                         " except for 'CS' stations for which 'HBA0', 'HBA1'"
                         " are also valid".format(arrband))
    # Use matrices in order to facilitate coord transformation.
    # axis=0
    stnpos = np.matrix(antfld[arrband]['POSITION']).T
    stnrot = np.matrix(antfld[subarr]['ROTATION_MATRIX'])
    stnrelpos = np.matrix(antfld[arrband]['REL_POS_X'])
    stnintilepos = np.matrix(hbadeltas)
    return stnpos, stnrot, stnrelpos, stnintilepos


def get_antset_params(stnid, antset):
    """Get configuration parameters for the antenna_set of the station.

    Parameters
    ----------
    stnid: str
        Station ID
    antset: str
        Antenna set specifier. E.g. LBA_INNER.

    Returns
    -------
    stnpos : array_like
        ITRF position vector of station's antset center in meters.
    stnrot : array_like
        Rotation matrix
    stnrelpos : array_like
        Positions of dual-pol antenna elements. (For HBA these refer to tiles)
        Axis 0 is antenna element number
        & axis 1 is cartesian 'x','y','z' coords,
        while values are in meters.
    stnintilepos : array_like
        Position of antennas within tile.
    """
    antfld = parseAntennaField(stnid)
    if antset.startswith('LBA'):
        arrband = 'LBA'
        hbadeltas = [0., 0., 0.]
    elif antset.startswith('HBA'):
        arrband = 'HBA'
        hbadeltas = parseiHBADeltasfile(stnid)
    stnpos = np.atleast_2d(antfld[arrband]['POSITION']).T
    stnrot = np.array(antfld[arrband]['ROTATION_MATRIX'])
    stnrelpos_x = np.array(antfld[arrband]['REL_POS_X'])
    stnrelpos_y = np.array(antfld[arrband]['REL_POS_Y'])
    stnintilepos = np.array(hbadeltas)
    nrhbaelems = len(antfld['HBA']['REL_POS_X'])
    if arrband == 'LBA' and nrhbaelems < 96:
        if '_' in antset:
            _LBA, antsettype = antset.split('_', 1)
        else:
            raise RuntimeError('Need to specify antset_type')
        if antsettype == 'INNER':
            stnrelpos = stnrelpos_x[:48]
        elif antsettype == 'OUTER':
            stnrelpos = stnrelpos_x[48:]
        if antsettype == 'X':
            # FIXME:
            stnrelpos = stnrelpos_x
        elif antsettype == 'Y':
            # FIXME:
            stnrelpos = stnrelpos_y
        elif antsettype == 'SPARSE_EVEN':
            stnrelpos = stnrelpos_x[0::2]
            stnrelpos = stnrelpos.reshape((2,24,3)).swapaxes(0,1).reshape((48,3))
        elif antsettype == 'SPARSE_ODD':
            stnrelpos = stnrelpos_x[1::2]
            stnrelpos = stnrelpos.reshape((2,24,3)).swapaxes(0,1).reshape((48,3))
    else:
        stnrelpos = stnrelpos_x
    return stnpos, stnrot, stnrelpos, stnintilepos


def antset_lonlat(stnid, antset):
    """\
    Get longitude, latitude position of antenna-set

    Parameters
    ----------
    stnid : str
        LOFAR station ID
    antset : str
        LOFAR station antenna-set, e.g. 'LBA' or 'HBA'
    Returns
    -------
    stnpos_lonlat : tuple
        Tuple with (longitude, latitude, height)
    """
    stnpos, stnrot, stnrelpos, stnintilepos = get_antset_params(stnid, antset)
    stnpos_lonlat = ITRF2lonlat(*stnpos.squeeze())
    return stnpos_lonlat


def list_stations(antenna_field_dir=ANTENNAFIELDDIR):
    """List all the available LOFAR station-ids.
    """
    dirlist = os.listdir(antenna_field_dir)
    stnId_list = []
    for f in dirlist:
        splitres = f.split("-AntennaField.conf")
        if len(splitres) == 2 and splitres[1] == '':
            stnId = splitres[0]
            stnId_list.append(stnId)
    return stnId_list


def get_tier_layouts(stnid, bandarr, coordsys='local'):
    """\
    Plot or print different kinds of array configurations

    Parameters
    ----------
    stnid : str
        The Station ID. Can be a LOFAR station 'SE607' or the international
        LOFAR 'ILT'.
    bandarr : str
        Band-array name: 'LBA', 'HBA', 'tile' or hexadecimal.
        'tile' refers to the layout within an HBA tile of 16 elements.
        The hexadecimal specifies a single-element-per-tile-on (SEPTON)
        mode element layout, where each hexadecimal digit refers to the ON
        element within each consecutive HBA tile.
    coordsys : str
        Coordinate system to use for output. Choose from 'local' or 'ITRF'.

    Returns
    -------
    pos : array
        (nrants, 3)  array of X,Y,Z ITRF positions.
    names : list
        List of unique name for each element in array.
    xyzlbls : tuple
        Labels for X,Y,Z position components.
    """
    if stnid == 'ILT':
        stns = list_stations()
        pos = []
        names = []
        for stn in stns:
            if stn == 'NenuFAR' and bandarr == 'HBA':
                continue
            if stn == 'UK902' and bandarr == 'LBA':
                continue
            names.append(stn)
            stnpos, stnrot, stnrelpos, stnintilepos = getArrayBandParams(
                stn, bandarr)
            pos.append(np.asarray(stnpos).squeeze().tolist())
        pos = np.array(pos)
    else:
        tier = None
        if bandarr == 'tile':
            tier = 'tile'
            bandarr = 'HBA'
        elif len(bandarr) == 96:
            tier = 'septon'
            elemmap = str2elementMap2(bandarr)
            bandarr = 'HBA'
        stnpos, stnrot, stnrelpos, stnintilepos = getArrayBandParams(
            stnid, bandarr)
        if coordsys == 'local':
            stnrelpos = stnrelpos * stnrot
            stnintilepos = stnintilepos * stnrot
        if tier == 'tile':
            _rel_pos = np.asarray(stnintilepos)
            nameprefix = 'elem'
        elif tier == 'septon':
            _rel_pos = np.asarray(stnrelpos)
            for tilenr, elemnr in enumerate(elemmap):
                _rel_pos[tilenr, :] += np.asarray(stnintilepos[elemnr]).squeeze()
            nameprefix = 'ant'
        else:
            if bandarr == 'HBA':
                nameprefix = 'tile'
            else:
                nameprefix = 'ant'
            _rel_pos = np.asarray(stnrelpos)
        pos = _rel_pos
        if coordsys == 'ITRF':
            pos += stnpos.T
        names = [nameprefix+str(elem) for elem in range(pos.shape[0])]

    # Set coord. labels:
    if coordsys == 'local':
        xlbl = 'p'
        ylbl = 'q'
        zlbl = 'r'
    else:
        # ITRF
        xlbl = 'X'
        ylbl = 'Y'
        zlbl = 'Z'

    return pos, names, (xlbl, ylbl, zlbl)


def print_layout_pos(pos, names):
    """\
    Print out layout positions
    """
    for idx in range(len(names)):
        print('{} {:10.2f} {:10.2f} {:10.2f}'.format(names[idx],
                                                     pos[idx, 0],
                                                     pos[idx, 1],
                                                     pos[idx, 2]))


def plot_layout_pos(pos, names, projection='3d', xyzlbls=('x','y','z'),
                    title='NoTitle'):
    """\
    Plot layout positions
    """
    (xlbl, ylbl, zlbl) = xyzlbls
    fig = plt.figure()

    ax = fig.gca(projection=projection)
    fontsize = 10
    if len(title) > 96:
        fontsize = 6
    ax.set_title(title, fontsize=fontsize)
    ax.set_xlabel(xlbl + ' [m]')
    ax.set_ylabel(ylbl + ' [m]')
    if projection == '3d':
        ax.set_zlabel(zlbl)
        ax.plot(pos[:, 0], pos[:, 1], pos[:, 2], '*')
    else:
        ax.plot(pos[:, 0], pos[:, 1], '*')
    for idx, name in enumerate(names):
        if projection == '3d':
            ax.text(pos[idx, 0], pos[idx, 1], pos[idx, 2], '   '+name,
                    fontsize=4)
        else:
            ax.text(pos[idx, 0], pos[idx, 1], '   ' + name, fontsize=4)
            ax.axis('equal')

    # Plot mean of positions
    geommean = np.mean(pos, axis=0)
    if projection != '3d':
        geommean = geommean[:2]
    ax.plot(*geommean, 'r+')
    ax.text(*geommean, 'Barycenter', fontsize=5)

    plt.show()


def cli_showarrconfs():
    stnId_list = list_stations()

    parser = argparse.ArgumentParser()
    parser.add_argument('view',
                        help="""Choose: 'print' or 'plot'.""")
    parser.add_argument('stnid',
                        help="""Station ID to process.
                        Choose from ILT or {}.
                        """.format(stnId_list))
    parser.add_argument('bandarr',
                        help="""Band array to process.
                        Choose from 'tile', '{}' or a SEPTON hex string."""
                        .format("', '".join(BANDARRS)))
    parser.add_argument('-c', '--coordsys', default='local',
                        help="""Coordinate system to use.
                        Choose from 'local' (default) or 'ITRF'.""")
    args = parser.parse_args()
    pos, names, xyzlbls = get_tier_layouts(args.stnid, args.bandarr, args.coordsys)
    if args.coordsys == 'local':
        projection = None
    bandarr_lbl = 'Band-array: '
    if len(args.bandarr) == 96:
        bandarr_lbl = 'SEPTON:\n'
    title = "Array configuration of {}\n{}{}\nCoord.sys.: {}"\
        .format(args.stnid, bandarr_lbl, args.bandarr, args.coordsys)
    if args.view == 'plot':
        if args.coordsys == 'ITRF':
            projection = '3d'
        else:
            projection = None
        plot_layout_pos(pos, names, projection, xyzlbls, title)
    else:
        print_layout_pos(pos, names)


if __name__ == '__main__':
    cli_showarrconfs()
