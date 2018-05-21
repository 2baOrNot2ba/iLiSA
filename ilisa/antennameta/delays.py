#!/usr/bin/env python
"""Module for handling the LOFAR antenna delay files.
  
"""
import sys
import os
import numpy as np
#import pkg_resources
#my_data = pkg_resources.resource_filename(__name__, "share/StaticMetaData")


STATICMETADATA = os.path.join(os.path.dirname(__file__),'share/StaticMetaData/')
CABLEDELAYDIR = STATICMETADATA
TABCOLUMNS = ('RCU','LBL_len','LBL_delay','LBH_len','LBH_delay','HBA_len','HBA_delay')
TABFORMATS = ('i4', 'f4',     'f4',       'f4',     'f4',       'f4',     'f4')

def _load_stn_cabledelay_file(f):
    lengthsNdelays = np.loadtxt(f, dtype={
                       'names': TABCOLUMNS,
                       'formats': TABFORMATS})
    return lengthsNdelays


def _stnid2filename(stnid):
    filename = stnid+'-'+'CableDelays.conf'
    return os.path.join(CABLEDELAYDIR,filename)

def _get_units(quantitystr):
    arr, qtype = quantitystr.split('_')
    if qtype == 'len':
        unit = 'm'
    elif qtype == 'delay':
        unit = 'ns'
    else:
        raise ValueError, "Unknown quantity type"
    return unit


def get_stn_cabledelays(stnid):
    f = _stnid2filename(stnid)
    lengthsNdelays = _load_stn_cabledelay_file(f)
    return lengthsNdelays


if __name__ == '__main__':
    stnid = sys.argv[1]
    quantity = sys.argv[2]
    if quantity not in TABCOLUMNS:
        raise ValueError, "Choose one of the following quantities: {}".format(TABCOLUMNS[1:])
    unit = _get_units(quantity)
    lengthsNdelays = get_stn_cabledelays(stnid)
    print("RCU [#] {} [{}]".format(quantity,unit))
    for row, rcu in enumerate(lengthsNdelays['RCU']):
        print("{}       {}".format(rcu, lengthsNdelays[quantity][row]))

