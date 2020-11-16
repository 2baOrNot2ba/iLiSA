#!/usr/bin/env python
import sys
import argparse
import ilisa.antennameta.export as export


def main(stnid, bandarr, quantity='all', output='default'):
    """Print out a station array antenna's metadata."""
    if quantity=='cfg' or quantity=='all':
        print("Quantity: Array configuration:")
        export.output_arrcfg_station(stnid, bandarr, output=output)
        # export.plot_arrayconfiguration(stnid, bandarr, "local")
    if quantity=='rot' or quantity=='all':
        print("Quantity: Array rotation matrix:")
        export.output_rotmat_station(stnid, bandarr, output=output)
    if quantity != 'cfg' and quantity != 'rot' and quantity != 'all':
        errmess = "Quantity {} not valid. Choose: cfg, rot, or all.".format(quantity)
        raise ValueError, errmess


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("stnid")
    parser.add_argument("bandarr")
    parser.add_argument("quantity", nargs='?', default='all')
    args = parser.parse_args()
    main(args.stnid, args.bandarr, quantity=args.quantity, output=sys.stdout)

