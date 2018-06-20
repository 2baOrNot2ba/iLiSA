#!/usr/bin/env python
import sys
import ilisa.observations.imaging as imaging


if __name__ == "__main__":

    cvcpath = sys.argv[1]
    cubeindex = float(sys.argv[2])
    calsrcdes = sys.argv[3]
    imaging.cvcimage(cvcpath, cubeindex, calsrcdes)

