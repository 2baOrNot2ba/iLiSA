#!/usr/bin/env python
import sys
import os
import argparse
import numpy
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import ilisa.observations.stationcontrol as stationcontrol
import ilisa.observations.dataIO as dataIO
import ilisa.observations.imaging as imaging


def plotbst(bstff):
    BSTdata, obsfileinfo = dataIO.readbstfolder(bstff)

    intg = obsfileinfo['integration']
    dur = obsfileinfo['duration']
    freqs = obsfileinfo['frequencies']
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
        if res.shape[0] > 1:
            plt.pcolormesh(freqs/1e6, ts/3600, res, norm=colors.LogNorm())
            plt.colorbar()
            plt.title('Mean (over RCUs) dynamicspectrum')
            plt.xlabel('Frequency [MHz]')
            plt.ylabel('Time [h]')
        else:
            # Only one integration so show it as 2D spectrum
            plt.plot(freqs/1e6, res[0,:])
            plt.yscale('log')
            plt.xlabel('Frequency [MHz]')
            plt.ylabel('Power [arb. unit]')
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
    obsfileinfo = xstobj.getobsfolderinfo()
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


def plotacc(args):
    accff = os.path.normpath(args.dataff)
    freq = args.freq
    dataobj = dataIO.CVCfiles(accff)
    data = dataobj.getdata()
    nrfiles = dataobj.getnrfiles()
    if nrfiles>1:
        data = data[0]
    sb, nqzone = stationcontrol.freq2sb(freq)
    while sb<512:
        print(sb,stationcontrol.sb2freq(sb,nqzone)/1e6)
        plt.pcolormesh(numpy.abs(data[sb]))
        plt.show()
        sb += 1


def plot_bsxst(args):
    bsxff = os.path.normpath(args.dataff)
    suf = whichst(bsxff)
    if suf=='acc':
        plotacc(args)
    if suf=='bst' or suf=='bst-357':
        plotbst(bsxff)
    elif suf=='sst':
        plotsst(bsxff, args.freq)
    elif suf=='xst' or suf=='xst-SEPTON':
        plotxst(bsxff)
    else:
        raise RuntimeError, "Not a bst, sst, or xst filefolder"

def image(args):
    """Image visibility-type data."""
    imaging.cvcimage(args.cvcpath, args.cubeindex, args.phaseref)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    parser_plot = subparsers.add_parser('plot', help='plot help')
    parser_plot.set_defaults(func=plot_bsxst)
    parser_plot.add_argument('dataff', help="acc, bst, sst or xst filefolder")
    parser_plot.add_argument('freq', type=float, nargs='?', default=None)

    parser_image = subparsers.add_parser('image', help='image help')
    parser_image.set_defaults(func=image)
    parser_image.add_argument('cvcpath', help="acc or xst filefolder")
    parser_image.add_argument('cubeindex', type=float, default=None)
    parser_image.add_argument('phaseref', type=str, default='Z')

    args = parser.parse_args()
    args.func(args)
