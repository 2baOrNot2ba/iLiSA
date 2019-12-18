#!/usr/bin/env python
import os
import argparse
import numpy
import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates

import ilisa.observations.dataIO as dataIO
import ilisa.observations.imaging as imaging
import ilisa.observations.modeparms as modeparms


def plotbst(bstff):
    """Plot BST data."""
    BSTdata, obsfileinfo = dataIO.readbstfolder(bstff)
    starttime = obsfileinfo['datetime']
    intg = obsfileinfo['integration']
    dur = obsfileinfo['duration']
    freqs = obsfileinfo['frequencies']
    ts = numpy.arange(0., dur, intg)
    ts = [starttime+datetime.timedelta(seconds=t) for t in ts]
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True)
    bstxplt = ax1.pcolormesh(ts, freqs/1e6, BSTdata['X'].T,
                             norm=colors.LogNorm())
    fig.colorbar(bstxplt, ax=ax1)
    ax1.set_ylabel('Frequency [MHz]')
    ax1.set_title('X-pol')
    bstyplt = ax2.pcolormesh(ts, freqs/1e6, BSTdata['Y'].T,
                             norm=colors.LogNorm())
    fig.colorbar(bstyplt, ax=ax2)
    ax2.set_title('Y-pol')
    fig.autofmt_xdate()

    ax2.xaxis.set_major_formatter( mdates.DateFormatter('%H:%M:%S'))
    ax2.set_xlabel('Time [UT]')
    ax2.set_ylabel('Frequency [MHz]')
    plt.suptitle('BST X,Y @ {}'.format(starttime))
    plt.show()


def plotsst(sstff, freqreq):
    """Plot SST data."""
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
            plt.pcolormesh(freqs/1e6, ts, res, norm=colors.LogNorm())
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
    """Plot XST data."""
    colorscale = None  # Colorscale for xst data plot (default None)
    if colorscale == 'log':
        normcolor = colors.LogNorm()
    else:
        normcolor = None
    xstobj = dataIO.CVCfiles(xstff)
    XSTdataset = xstobj.getdata()
    file_ids = xstobj.scanrecinfo.get_file_ids()
    for sbstepidx in range(len(XSTdataset)):
        obsinfo = xstobj.scanrecinfo.obsinfos[file_ids[sbstepidx]]
        intg = obsinfo.integration
        dur = obsinfo.duration_scan

        freq = obsinfo.get_recfreq()
        ts = numpy.arange(0., dur, intg)
        XSTdata = XSTdataset[sbstepidx]
        for tidx in range(XSTdata.shape[0]):
            print("Kill plot window for next plot...")
            plt.imshow(numpy.abs(XSTdata[tidx,...]), norm=normcolor,
                       interpolation='none')
            plt.title(
"""Time (from start {}) {}s
@ freq={} MHz""".format(obsinfo.get_starttime(), ts[tidx], freq/1e6))
            plt.xlabel('RCU [#]')
            plt.ylabel('RCU [#]')
            plt.colorbar()
            plt.show()


def plotacc(accff):
    """Plot of ACC folder files."""
    dataobj = dataIO.CVCfiles(accff)
    data = dataobj.getdata()
    if args.freq is None:
        args.freq = 0.0
    sb, nqzone = modeparms.freq2sb(args.freq)
    for fileidx in range(0, dataobj.getnrfiles()):
        while sb<512:
            fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True)
            absdatplt = ax1.pcolormesh(numpy.abs(data[fileidx][sb]))
            ax1.set_title('Abs value')
            #ax1.set_xlabel('RCU [#]')
            ax1.set_ylabel('RCU [#]')
            fig.colorbar(absdatplt, ax=ax1)
            angdatplt = ax2.pcolormesh(numpy.angle(data[fileidx][sb]), cmap=plt.get_cmap('hsv'))
            ax2.set_title('Phase value')
            ax2.set_xlabel('RCU [#]')
            ax2.set_ylabel('RCU [#]')
            fig.colorbar(angdatplt, ax=ax2)
            plt.suptitle('Station element covariance. Time: {}UT, SB: {}'\
                      .format(dataobj.samptimeset[fileidx][sb], sb))
            fig.show()
            inpres = raw_input('Press return for next plot... ')
            sb += 1


def _plot_bsxst(args):
    if lofar_datatype =='acc':
        plotacc(args.dataff)
    if lofar_datatype=='bst' or lofar_datatype=='bst-357':
        plotbst(args.dataff)
    elif lofar_datatype=='sst':
        plotsst(args.dataff, args.freq)
    elif lofar_datatype=='xst' or lofar_datatype=='xst-SEPTON':
        plotxst(args.dataff)
    else:
        raise RuntimeError("Not a bst, sst, or xst filefolder")


def image(args):
    """Image visibility-type data."""
    if lofar_datatype != 'acc' and lofar_datatype != 'xst':
        raise RuntimeError("Datafolder '{}'\n not ACC or XST type data."\
                           .format(args.dataff))
    cvcobj = dataIO.CVCfiles(args.dataff)
    CVCdataset = cvcobj.getdata()
    for fileidx in range(args.filenr, cvcobj.getnrfiles()):
        integration = cvcobj.scanrecinfo.get_integration()
        CVCdata = CVCdataset[fileidx]
        if CVCdata.ndim == 2:
            intgs = 1
        else:
            intgs = CVCdata.shape[0]
        for tidx in range(args.sampnr, intgs):
            ll, mm, skyimages, polrep, t, freq, stnid, phaseref = \
                imaging.cvcimage(cvcobj, fileidx, tidx, args.phaseref,
                                 pbcor=args.correctpb)
            imaging.plotskyimage(ll, mm, skyimages, polrep, t, freq, stnid, phaseref,
                                 integration, pbcor=args.correctpb)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    parser_plot = subparsers.add_parser('plot', help='plot help')
    parser_plot.set_defaults(func=_plot_bsxst)
    parser_plot.add_argument('dataff', help="acc, bst, sst or xst filefolder")
    parser_plot.add_argument('freq', type=float, nargs='?', default=None)

    parser_image = subparsers.add_parser('image', help='image help')
    parser_image.set_defaults(func=image)
    parser_image.add_argument('dataff', help="acc or xst filefolder")
    parser_image.add_argument('-p','--phaseref', type=str, default=None)
    parser_image.add_argument('-n','--filenr', type=int, default=0)
    parser_image.add_argument('-s','--sampnr', type=int, default=0)
    parser_image.add_argument('-c','--correctpb', help="Correct for primary beam",
                              action="store_true")

    args = parser.parse_args()
    args.dataff = os.path.normpath(args.dataff)
    lofar_datatype = dataIO.datafolder_type(args.dataff)
    args.func(args)
