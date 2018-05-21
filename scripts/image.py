#!/usr/bin/env python
import sys
import ilisa.observations.imaging as imaging


if __name__ == "__main__":

    cvcpath = sys.argv[1]
    if '_acc' in cvcpath:
        cvctype = 'acc'
    else:
        cvctype = 'xst'
    cubeindex = float(sys.argv[2])
    #stnid = sys.argv[3] #"SE607" 
    calsrcdes = sys.argv[3]
    imaging.cvcimage(cvctype, cvcpath, cubeindex, calsrcdes)

