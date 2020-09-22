#!/usr/bin/env python
"""A module to read LOFAR antenna field files.
The ROTATION_MATRIX field is the rotation matrix of the station;
it maps the the local station coordinates to the ITRF coordinates
such that

    r_ITRF = ROTATION_MATRIX * r_stn where r are 3-element column vectors.
"""

import numpy as np
import os
from os.path import dirname

__version__ = 0.2

STATICMETADATA = os.path.join(dirname(__file__),'share/StaticMetaData/')
ANTENNAFIELDDIR = STATICMETADATA
IHBADELTASDIR = os.path.join(dirname(__file__),'share/iHBADeltas/')
IHBADELTASDIR = STATICMETADATA

COMMENT_CHAR = '#'
AFfileNameType = 3


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
    AntFldData={'LBA': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]},
                'HBA': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[],'REL_POS':[]},
                'HBA0': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[]},
                'HBA1': {'NORMAL_VECTOR': [], 'ROTATION_MATRIX':[],'POSITION':[]}
               }
    
    try:
        f = open(filename)
    except IOError:
        print("Error: "+filename+" does not exist.")
        raise
    line = f.readline()
    while line:
        if COMMENT_CHAR in line:
            line, comment = line.split(COMMENT_CHAR, 1)
        if "HBA0" in line:
            AntBand = "HBA0"
        elif "HBA1" in line:
            AntBand = "HBA1"
        elif "HBA" in line:
            AntBand = "HBA"
        elif "LBA" in line:
            AntBand = "LBA"
        else:
            line=f.readline()
            continue
        where, rest = line.split(AntBand, 1)
        where = where.strip()
        if where == '':
            # Read absolute position of station origin
            line = f.readline()
            elementposshape, elementposLine = line.split(' ',1)
            elementposLine = elementposLine.strip('[] \n').split()
            position = [float(v) for v in elementposLine]
            AntFldData[AntBand]['POSITION'] = position
            if AntBand != "HBA0" and AntBand != "HBA1":
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
                    AntFldData[AntBand]['REL_POS'].append(posxpol) 
                    # Note: Skip ypol as it is identical to xpol
                # Read ending ']' line
                line=f.readline()
        elif where == 'NORMAL_VECTOR':
            line = f.readline()
            elementposshape, elementposLine = line.split(' ',1)
            elementposLine = elementposLine.strip('[] \n').split()
            nrmv = [float(v) for v in elementposLine]
            AntFldData[AntBand][where]=nrmv
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
                  AntFldData[AntBand][where].append(row)
              # Read ending ']' line
              line = f.readline()
        line = f.readline()
    return AntFldData


def parseiHBADeltasfile(stationName):
    """Parse iHBADelta file."""
    iHBADeltasdata = []
    filepath = _getiHBADeltafile(stationName)
    f = open(filepath)
    line=f.readline()
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
    errmess = "Array band not valid. Only 'LBA', 'HBA' are valid, "\
            + "except for 'CS' stations for which 'HBA0', 'HBA1' are also valid."
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
        else:
            raise ValueError(errmess)
        arrband = 'HBA'
        hbadeltas = parseiHBADeltasfile(stnid)
    else:
        raise ValueError(errmess)
    # Use matrices in order to facilitate coord transformation.
    # axis=0
    stnpos = np.matrix(antfld[arrband]['POSITION']).T
    stnrot = np.matrix(antfld[subarr]['ROTATION_MATRIX'])
    stnrelpos = np.matrix(antfld[arrband]['REL_POS'])
    stnintilepos = np.matrix(hbadeltas)
    return stnpos, stnrot, stnrelpos, stnintilepos


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


BANDARRS = ['LBA', 'HBA']
