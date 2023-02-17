import argparse
import numpy as np
from matplotlib import pyplot as plt

from ilisa.antennameta.antennafieldlib import list_stations, \
    getArrayBandParams, BANDARRS, parseAntennaField


def inspect_arrayconfiguration(view, stnid, bandarr, coordsys='local'):
    """Plot different kinds of array configurations."""
    if stnid == 'ILT':
        stns = list_stations()
        pos = []
        names = []
        for stn in stns:
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
        stnpos, stnrot, stnrelpos, stnintilepos = getArrayBandParams(
            stnid, bandarr)
        if coordsys == 'local':
            stnrelpos = stnrelpos * stnrot
            stnintilepos = stnintilepos * stnrot
        if tier == 'tile':
            pos = np.asarray(stnintilepos)
            nameprefix = 'elem'
        else:
            if bandarr == 'HBA':
                nameprefix = 'tile'
            else:
                nameprefix = 'ant'
            pos = np.asarray(stnrelpos)

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

    if view == 'print':
        for idx in range(len(names)):
            print(names[idx], pos[idx, 0], pos[idx, 1], pos[idx, 2])
    else:
        # Plot
        fig = plt.figure()
        projection = '3d'
        if coordsys == 'local':
            projection = None
        ax = fig.gca(projection=projection)
        ax.set_title("Array configuration of {} {} in coordsys {}"
                     .format(stnid, bandarr, coordsys))
        ax.set_xlabel(xlbl)
        ax.set_ylabel(ylbl)
        if projection == '3d':
            ax.set_zlabel(zlbl)
            ax.plot(pos[:, 0], pos[:, 1], pos[:, 2], '*')
        else:
            ax.plot(pos[:, 0], pos[:, 1], '*')
        for idx, name in enumerate(names):
            if projection == '3d':
                ax.text(pos[idx, 0], pos[idx, 1], pos[idx, 2], '   '+name,
                        fontsize=6)
                ax.set_zlim(-0.5, 0.5)
            else:
                ax.text(pos[idx, 0], pos[idx, 1], '   ' + name, fontsize=6)
                ax.axis('equal')
        plt.show()


def printantennafieldfile(stnid):
    afd = parseAntennaField(stnid)
    print(afd)


def print_arrayconfiguration(stnid, bandarr):
    stnpos, stnrot, stnrelpos, stnintilepos = \
            getArrayBandParams(stnid, bandarr)
    print("Position:\n{}".format(stnpos))
    print("Rotation:\n{}".format(stnrot))
    print("Relative positions:\n{}".format(stnrelpos))
    print("In-tile positions:\n{}".format(stnintilepos))


if __name__ == '__main__':
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
                        Choose from 'tile' or {}.""".format(BANDARRS))
    parser.add_argument('-c', '--coordsys',
                        help="""Coordinate system to use.
                        Choose from 'local' or 'ITRF'.
                        """)
    args = parser.parse_args()
    inspect_arrayconfiguration(args.view, args.stnid, args.bandarr,
                               args.coordsys)
