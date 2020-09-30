#!/usr/bin/env python
"""Module for handling the LOFAR antenna delay files.
  
"""
import sys
import os
import numpy as np


STATICMETADATA = os.path.join(os.path.dirname(__file__), 'share/StaticMetaData/')
CABLEDELAYDIR = STATICMETADATA
TABCOLUMNS = ('RCU', 'LBL_len', 'LBL_delay', 'LBH_len', 'LBH_delay', 'HBA_len',
              'HBA_delay')
TABFORMATS = ('i4', 'f4',     'f4',       'f4',     'f4',       'f4',     'f4')


def _load_stn_cabledelay_file(f):
    lengths_n_delays = np.loadtxt(f, dtype={
                       'names': TABCOLUMNS,
                       'formats': TABFORMATS})
    return lengths_n_delays


def _stnid2filename(stnid):
    filename = stnid+'-'+'CableDelays.conf'
    return os.path.join(CABLEDELAYDIR, filename)


def _get_units(quantitystr):
    arr, qtype = quantitystr.split('_')
    if qtype == 'len':
        unit = 'm'
    elif qtype == 'delay':
        unit = 'ns'
    else:
        raise ValueError("Unknown quantity type")
    return unit


def get_stn_cabledelays(stnid):
    f = _stnid2filename(stnid)
    lengths_n_delays = _load_stn_cabledelay_file(f)
    return lengths_n_delays


def main():
    stnid = sys.argv[1]
    quantity = sys.argv[2]
    if quantity not in TABCOLUMNS:
        raise ValueError("Choose one of the following quantities: {}".format(TABCOLUMNS[1:]))
    unit = _get_units(quantity)
    lengths_n_delays = get_stn_cabledelays(stnid)
    print("RCU [#] {} [{}]".format(quantity, unit))
    for row, rcu in enumerate(lengths_n_delays['RCU']):
        print("{}       {}".format(rcu, lengths_n_delays[quantity][row]))


if __name__ == '__main__':
    main()
