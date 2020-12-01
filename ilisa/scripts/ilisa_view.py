import os
import argparse
import numpy
import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates

import ilisa.observations.dataIO as dataIO
import ilisa.calim.imaging as imaging
import ilisa.observations.modeparms as modeparms


def plotbst(bstff, pol_stokes=True):
    """Plot BST data."""
    BSTdata, obsfileinfo = dataIO.readbstfolder(bstff)
    starttime = obsfileinfo['datetime']
    intg = obsfileinfo['integration']
    dur = obsfileinfo['duration']
    freqs = obsfileinfo['frequencies']
    pointing = obsfileinfo['pointing']

    ts = numpy.arange(0., dur, intg)
    ts = [starttime+datetime.timedelta(seconds=t) for t in ts]
    fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True)
    data2plot_p_name, data2plot_p = BSTdata['X'].T, 'X-pol'
    data2plot_q_name, data2plot_q = BSTdata['Y'].T, 'Y-pol'
    if pol_stokes:
        data2plot_p_name = 'Stokes I'
        data2plot_p = BSTdata['X'].T + BSTdata['Y'].T
        data2plot_q_name = 'Stokes Q'
        data2plot_q = BSTdata['X'].T - BSTdata['Y'].T
    bstxplt = ax1.pcolormesh(ts, freqs/1e6, data2plot_p,
                             norm=colors.LogNorm())
    fig.colorbar(bstxplt, ax=ax1)
    ax1.set_ylabel('Frequency [MHz]')
    ax1.set_title('{}, pointing {},{},{}'.format(data2plot_p_name, *pointing))
    bstyplt = ax2.pcolormesh(ts, freqs/1e6, data2plot_q,
                             norm=colors.LogNorm())
    fig.colorbar(bstyplt, ax=ax2)
    ax2.set_title('{}, pointing {},{},{}'.format(data2plot_q_name, *pointing))
    fig.autofmt_xdate()

    ax2.xaxis.set_major_formatter( mdates.DateFormatter('%H:%M:%S'))
    ax2.set_xlabel('Time [UT]')
    ax2.set_ylabel('Frequency [MHz]')
    plt.suptitle('BST X,Y @ {}'.format(starttime))
    plt.show()


def plotsst(sstff, freqreq):
    """Plot SST data."""
    SSTdata, obsfolderinfo = dataIO.readsstfolder(sstff)
    scanrecinfo = dataIO.ScanRecInfo()
    scanrecinfo.read_scanrec(sstff)
    starttime = obsfolderinfo['datetime']
    SSTdata = numpy.array(SSTdata)
    freqs = modeparms.rcumode2sbfreqs(obsfolderinfo['rcumode'])
    if freqreq:
        sbreq = int(numpy.argmin(numpy.abs(freqs-freqreq)))
        show = 'persb'
    else:
        show = 'mean'
    intg = obsfolderinfo['integration']
    dur = obsfolderinfo['duration']
    ts = [starttime + datetime.timedelta(seconds=td) for td in 
                                                numpy.arange(0., dur, intg)]
    if show == 'mean':
        meandynspec = numpy.mean(SSTdata, axis=0)
        res = meandynspec
        if res.shape[0] > 1:
            plt.pcolormesh(freqs/1e6, ts, res, norm=colors.LogNorm())
            plt.colorbar()
            plt.title('Mean (over RCUs) dynamicspectrum\n'
                      +'Starttime: {} Station: {}'
                      .format(starttime, scanrecinfo.stnid))
            plt.xlabel('Frequency [MHz]')
            plt.ylabel('Time [h]')
        else:
            # Only one integration so show it as 2D spectrum
            plt.plot(freqs/1e6, res[0, :])
            plt.yscale('log')
            plt.xlabel('Frequency [MHz]')
            plt.ylabel('Power [arb. unit]')
    elif show == 'persb':
        ampVStime = True
        res = SSTdata[:, :, sbreq]
        resX = res[0::2, :]
        resY = res[1::2, :]
        plt.subplot(211)
        if ampVStime:
            plt.plot(ts, numpy.transpose(resX))
            # plt.gcf().autofmt_xdate()
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
    obs_ids = xstobj.scanrecinfo.get_obs_ids()
    # Assume freq sweep over nr of files, so filenr is also sb. 
    for sbstepidx in range(xstobj.getnrfiles()):
        obsinfo = xstobj.scanrecinfo.obsinfos[obs_ids[sbstepidx]]
        intg = obsinfo.integration
        dur = obsinfo.duration_scan

        freq = obsinfo.get_recfreq()
        ts = numpy.arange(0., dur, intg)
        xstdata = xstobj.covcube_fb(sbstepidx, crlpolrep=None)
        for tidx in range(xstdata.shape[0]):
            print("Kill plot window for next plot...")
            plt.imshow(numpy.abs(xstdata[tidx,...]), norm=normcolor,
                       interpolation='none')
            plt.title("""Time (from start {}) {}s
                      @ freq={} MHz""".format(obsinfo.get_starttime(),
                      ts[tidx], freq/1e6))
            plt.xlabel('RCU [#]')
            plt.ylabel('RCU [#]')
            plt.colorbar()
            plt.show()


def plotacc(accff, freqreq=None):
    """Plot of ACC folder files."""
    dataobj = dataIO.CVCfiles(accff)
    if freqreq is None:
        freqreq = 0.0
    sb, _nqzone = modeparms.freq2sb(freqreq)
    for fileidx in range(0, dataobj.getnrfiles()):
        filecvc = dataobj.covcube_fb(fileidx, crlpolrep=None)
        while sb < 512:
            fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True)
            absdatplt = ax1.pcolormesh(numpy.abs(filecvc[sb]))
            ax1.set_title('Abs value')
            ax1.set_ylabel('RCU [#]')
            fig.colorbar(absdatplt, ax=ax1)
            angdatplt = ax2.pcolormesh(numpy.angle(filecvc[sb]),
                                       cmap=plt.get_cmap('hsv'))
            ax2.set_title('Phase value')
            ax2.set_xlabel('RCU [#]')
            ax2.set_ylabel('RCU [#]')
            fig.colorbar(angdatplt, ax=ax2)
            plt.suptitle('Station element covariance. Time: {}UT, SB: {}'\
                                .format(dataobj.samptimeset[fileidx][sb], sb))
            plt.show()
            sb += 1


def plot_bsxst(args):
    lofar_datatype = dataIO.datafolder_type(args.dataff)
    if lofar_datatype =='acc':
        plotacc(args.dataff, args.freq)
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
    polrep = 'stokes'
    lofar_datatype = dataIO.datafolder_type(args.dataff)
    fluxperbeam = not args.fluxpersterradian
    if lofar_datatype != 'acc' and lofar_datatype != 'xst':
        raise RuntimeError("Datafolder '{}'\n not ACC or XST type data."
                           .format(args.dataff))
    cvcobj = dataIO.CVCfiles(args.dataff)
    calibrated = False
    if cvcobj.scanrecinfo.calibrationfile:
        calibrated = True
    stnid = cvcobj.scanrecinfo.get_stnid()
    for fileidx in range(args.filenr, cvcobj.getnrfiles()):
        integration = cvcobj.scanrecinfo.get_integration()
        cvpol = cvcobj.covcube_fb(fileidx)
        intgs = cvpol.shape[-3]
        for tidx in range(args.sampnr, intgs):
            t = cvcobj.samptimeset[fileidx][tidx]
            freq = cvcobj.freqset[fileidx][tidx]
            skyimages, ll, mm, phaseref = \
                imaging.cvc_image(cvcobj, fileidx, tidx, args.phaseref,
                                  polrep=polrep, pbcor=args.correctpb,
                                  fluxperbeam=fluxperbeam)
            imaging.plotskyimage(ll, mm, skyimages, polrep, t, freq, stnid,
                                 integration, phaseref, calibrated,
                                 pbcor=args.correctpb, maskhrz=False,
                                 fluxperbeam=fluxperbeam)

def nfimage(args):
    """
    Near field image.
    """
    polrep = 'S0'
    lofar_datatype = dataIO.datafolder_type(args.dataff)
    if lofar_datatype != 'acc' and lofar_datatype != 'xst':
        raise RuntimeError("Datafolder '{}'\n not ACC or XST type data."
                           .format(args.dataff))
    cvcobj = dataIO.CVCfiles(args.dataff)
    stnid = cvcobj.scanrecinfo.get_stnid()
    for fileidx in range(args.filenr, cvcobj.getnrfiles()):
        integration = cvcobj.scanrecinfo.get_integration()
        cvpol = cvcobj.covcube_fb(fileidx)
        intgs = cvpol.shape[-3]
        for tidx in range(args.sampnr, intgs):
            xx, yy, nfimages = imaging.nearfield_grd_image(cvcobj, fileidx,
                                                           tidx)
            t = cvcobj.samptimeset[fileidx][tidx]
            freq = cvcobj.freqset[fileidx][tidx]
            imaging.plotskyimage(xx, yy, nfimages, polrep, t, freq, stnid,
                                 integration, maskhrz=False)


def main():
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(help='sub-command help')

    parser_plot = subparsers.add_parser('plot', help='plot help')
    parser_plot.set_defaults(func=plot_bsxst)
    parser_plot.add_argument('dataff', help="acc, bst, sst or xst filefolder")
    parser_plot.add_argument('freq', type=float, nargs='?', default=None)

    parser_image = subparsers.add_parser('image', help='image help')
    parser_image.set_defaults(func=image)
    parser_image.add_argument('dataff', help="acc or xst filefolder")
    parser_image.add_argument('-p', '--phaseref', type=str, default=None)
    parser_image.add_argument('-n', '--filenr', type=int, default=0)
    parser_image.add_argument('-s', '--sampnr', type=int, default=0)
    parser_image.add_argument('-c', '--correctpb',
                              help="Correct for primary beam",
                              action="store_true")
    parser_image.add_argument('-f', '--fluxpersterradian',
                              help="Normalize flux per sterradian",
                              action="store_true")

    parser_image = subparsers.add_parser('nfimage', help='nearfield image')
    parser_image.set_defaults(func=nfimage)
    parser_image.add_argument('dataff', help="acc or xst filefolder")
    parser_image.add_argument('-n', '--filenr', type=int, default=0)
    parser_image.add_argument('-s', '--sampnr', type=int, default=0)

    args = parser.parse_args()
    args.dataff = os.path.normpath(args.dataff)
    args.func(args)


if __name__ == "__main__":
    main()
