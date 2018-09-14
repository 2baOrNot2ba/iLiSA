#!/usr/bin/env python
import sys
import os
import argparse
import numpy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import ilisa.observations.stationcontrol as stationcontrol
import ilisa.observations.dataIO as dataIO


def plotbst(bstff):
    BSTdata, obsfileinfo = dataIO.readbstfolder(bstff)
    sbs = obsfileinfo['subbands']
    intg = obsfileinfo['integration']
    dur = obsfileinfo['duration']
    freqinband = stationcontrol.rcumode2sbfreqs(obsfileinfo['rcumode'])
    freqs = freqinband[obsfileinfo['sblo']:(obsfileinfo['sbhi']+1)]
    ts = numpy.arange(0., dur, intg)
    plt.subplot(211)
    plt.pcolormesh(freqs/1e6, ts/3600, BSTdata['X'],
                   norm=colors.LogNorm())
    plt.title('X-pol')
    plt.subplot(212)
    plt.pcolormesh(freqs/1e6, ts/3600, BSTdata['Y'],
                   norm=colors.LogNorm())
    plt.title('Y-pol')
    plt.xlabel('frequency [MHz]')
    plt.ylabel('Time [h]')
    plt.show()


def plotsst(sstff, freqreq):
    SSTdata, obsfolderinfo = dataIO.readsstfolder(sstff)
    SSTdata = numpy.array(SSTdata)
    freqs = stationcontrol.rcumode2sbfreqs(obsfolderinfo['rcumode'])
    if  freqreq:
        sbreq = int(numpy.argmin(numpy.abs(freqs-freqreq)))
        show = 'persb'
    else:
        show = 'mean'
    intg = obsfolderinfo['integration']
    dur = obsfolderinfo['duration']
    ts = numpy.arange(0., dur, intg)
    if show == 'mean':
        meandynspec = numpy.mean(SSTdata, axis=0)
        res = meandynspec
        plt.pcolormesh(freqs/1e6, ts/3600, res, norm=colors.LogNorm())
        plt.colorbar()
        plt.title('Mean (over RCUs) dynamicspectrum')
        plt.xlabel('frequency [MHz]')
        plt.ylabel('Time [h]')
    elif show == 'persb':
        res = SSTdata[:,:,sbreq]
        resX = res[0::2,:]
        resY = res[1::2,:]
        plt.subplot(211)
        plt.plot(ts/3600,numpy.transpose(resX))
        plt.xlabel('Time [h]')
        plt.title('X pol')
        plt.subplot(212)
        plt.pcolormesh(ts/3600, numpy.arange(96), resY, norm=colors.LogNorm())
        plt.ylabel('rcu [nr]')
        plt.xlabel('Time [h]')
        plt.title('Y pol')
        plt.suptitle('Freq {} MHz'.format(freqs[sbreq]/1e6))
    plt.show()


def plotxst(xstff):
    xstobj = dataIO.CVCfiles(xstff)
    XSTdata = xstobj.getdata()
    obsfileinfo = xstobj.getobsfileinfo()
    sb = obsfileinfo['subband']
    intg = int(obsfileinfo['integration'])
    dur = int(obsfileinfo['duration'])
    
    freq = stationcontrol.sb2freq(sb,stationcontrol.rcumode2NyquistZone(obsfileinfo['rcumode']))
    ts = numpy.arange(0., dur, intg)
    
    for tidx in range(XSTdata.shape[0]):
      print "Kill plot window for next plot..."
      plt.imshow(numpy.abs(XSTdata[tidx,...]), norm=colors.LogNorm(), interpolation='none')
      plt.title("Time (from start) {}s".format(ts[tidx]))
      plt.xlabel('RCU [#]')
      plt.ylabel('RCU [#]')
      plt.show()


def whichst(bsxff):
    suf = bsxff.split('_')[-1]
    return suf


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument('bsxff', help="bst, sst or xst filefolder")
    parser.add_argument("freq", type=float, nargs='?', default=None)
    args = parser.parse_args()
    bsxff = os.path.normpath(args.bsxff)
    suf = whichst(bsxff)
    if suf=='bst':
        plotbst(bsxff)
    elif suf=='sst':
        plotsst(bsxff, args.freq)
    elif suf=='xst':
        plotxst(bsxff)
    else:
        raise RuntimeError, "Not a bst, sst, or xst filefolder"

