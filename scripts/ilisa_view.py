#!/usr/bin/env python
import os
import argparse
import numpy
import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as colors

import ilisa.observations.dataIO as dataIO
import ilisa.observations.imaging as imaging
import ilisa.observations.modeparms as modeparms


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
    starttime = obsfolderinfo['datetime']
    SSTdata = numpy.array(SSTdata)
    freqs = modeparms.rcumode2sbfreqs(obsfolderinfo['rcumode'])
    if  freqreq:
        sbreq = int(numpy.argmin(numpy.abs(freqs-freqreq)))
        show = 'persb'
    else:
        show = 'mean'
    intg = obsfolderinfo['integration']
    dur = obsfolderinfo['duration']
    ts = [starttime + datetime.timedelta(seconds=td) for td in numpy.arange(0., dur, intg)]
    if show == 'mean':
        meandynspec = numpy.mean(SSTdata, axis=0)
        res = meandynspec
        if res.shape[0] > 1:
            plt.pcolormesh(freqs/1e6, starttime + tds, res, norm=colors.LogNorm())
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
        ampVStime = True
        res = SSTdata[:,:,sbreq]
        resX = res[0::2,:]
        resY = res[1::2,:]
        plt.subplot(211)
        if ampVStime:
            plt.plot(ts, numpy.transpose(resX))
            #plt.gcf().autofmt_xdate()
        else:
            plt.pcolormesh(ts, numpy.arange(96), resX, norm=colors.LogNorm())
        plt.title('X pol')
        plt.subplot(212)
        if ampVStime:
            plt.plot(ts, numpy.transpose(resY))
            plt.gcf().autofmt_xdate()
        else:
            plt.pcolormesh(ts, numpy.arange(96), resY, norm=colors.LogNorm())
            plt.ylabel('rcu [nr]')
        plt.xlabel('Time [UT]')
        plt.title('Y pol')
        plt.suptitle('Freq {} MHz'.format(freqs[sbreq]/1e6))
    plt.show()


def plotxst(xstff):
    xstobj = dataIO.CVCfiles(xstff)
    XSTdataset = xstobj.getdata()
    for sbstepidx in range(len(XSTdataset)):
        obsinfo = xstobj.stnsesinfo.obsinfos[sbstepidx]
        intg = float(obsinfo.rspctl_cmd['integration'])
        dur = float(obsinfo.rspctl_cmd['duration'])

        freq = obsinfo.get_recfreq()
        ts = numpy.arange(0., dur, intg)
        XSTdata = XSTdataset[sbstepidx]
        for tidx in range(XSTdata.shape[0]):
            print("Kill plot window for next plot...")
            plt.imshow(numpy.abs(XSTdata[tidx,...]), norm=colors.LogNorm(),
                       interpolation='none')
            plt.title("Time (from start {}) {}s @ freq={} MHz".format(obsinfo.starttime,
                                                                      ts[tidx], freq/1e6))
            plt.xlabel('RCU [#]')
            plt.ylabel('RCU [#]')
            plt.show()


def whichst(bsxff):
    suf = bsxff.split('_')[-1]
    return suf


def plotacc(args):
    accff = os.path.normpath(args.dataff)
    dataobj = dataIO.CVCfiles(accff)
    obsfolderinfo = dataobj.getobsfolderinfo()
    data = dataobj.getdata()
    nrfiles = dataobj.getnrfiles()
    if nrfiles>1:
        data = data[0]
    if args.freq is None:
        args.freq = 0.0
    sb, nqzone = modeparms.freq2sb(args.freq)
    while sb<512:
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
        raise RuntimeError("Not a bst, sst, or xst filefolder")


def image(args):
    """Image visibility-type data."""
    cvcobj = dataIO.CVCfiles(args.cvcpath)
    CVCdataset = cvcobj.getdata()
    for fileidx in range(args.filenr, cvcobj.getnrfiles()):
        integration = cvcobj.stnsesinfo.get_integration()
        CVCdata = CVCdataset[fileidx]
        if CVCdata.ndim == 2:
            intgs = 1
        else:
            intgs = CVCdata.shape[0]
        for tidx in range(args.sampnr, intgs):
            ll, mm, skyimages, t, freq, stnid, phaseref = \
                imaging.cvcimage(cvcobj, fileidx, tidx, args.phaseref)
            imaging.plotskyimage(ll, mm, skyimages, t, freq, stnid, phaseref, integration)


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
    parser_image.add_argument('-p','--phaseref', type=str, default=None)
    parser_image.add_argument('-n','--filenr', type=int, default=0)
    parser_image.add_argument('-s','--sampnr', type=int, default=0)

    args = parser.parse_args()
    args.func(args)
