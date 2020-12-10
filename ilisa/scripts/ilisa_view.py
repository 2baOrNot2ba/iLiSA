import os
import argparse
import numpy
import datetime
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates

import ilisa.monitorcontrol.data_io as dataIO
import ilisa.monitorcontrol.modeparms as modeparms


def plotbst(bstff, pol_stokes=True):
    """Plot BST data."""
    BSTdata, obsfileinfo = dataIO.readbstfolder(bstff)
    stnid = obsfileinfo['station_id']
    starttime = obsfileinfo['datetime']
    intg = obsfileinfo['integration']
    dur = obsfileinfo['duration_scan']
    freqs = obsfileinfo['frequencies']
    pointing = obsfileinfo['pointing']

    ts = numpy.arange(0., dur, intg)
    ts = [starttime+datetime.timedelta(seconds=t) for t in ts]
    fig, (ax_p, ax_q) = plt.subplots(2, 1, sharex=True, sharey=True)
    data2plot_p_name, data2plot_p = BSTdata['X'].T, 'X-pol'
    data2plot_q_name, data2plot_q = BSTdata['Y'].T, 'Y-pol'
    if pol_stokes:
        stokes_norm = True
        data2plot_p_name = 'Stokes I'
        data2plot_p = BSTdata['X'].T + BSTdata['Y'].T
        data2plot_q_name = '(antenna) Stokes Q'
        data2plot_q = BSTdata['X'].T - BSTdata['Y'].T
        data2plot_q_unit = 'Signed flux [arb. units]'
        if stokes_norm:
            data2plot_q = data2plot_q / data2plot_p
            data2plot_q_name = 'Stokes q'
            data2plot_q_unit = 'Signed relative flux [%]'
            bstplt_q = ax_q.pcolormesh(ts, freqs / 1e6, data2plot_q,
                                       cmap='RdBu_r',
                                       norm=colors.SymLogNorm(linthresh=1e-3))
        else:
            bstplt_q = ax_q.pcolormesh(ts, freqs / 1e6, data2plot_q,
                                       cmap='RdBu_r',
                                       norm=colors.SymLogNorm(linthresh=1e2))

    bstplt_p = ax_p.pcolormesh(ts, freqs/1e6, data2plot_p,
                              norm=colors.LogNorm())
    cbar_p = fig.colorbar(bstplt_p, ax=ax_p)
    cbar_p.set_label('Flux [arb. units]')
    ax_p.set_ylabel('Frequency [MHz]')
    ax_p.set_title('{}'.format(data2plot_p_name))

    cbar_q = fig.colorbar(bstplt_q, ax=ax_q)
    cbar_q.set_label(data2plot_q_unit)
    ax_q.set_title('{}'.format(data2plot_q_name))
    fig.autofmt_xdate()

    ax_q.xaxis.set_major_formatter( mdates.DateFormatter('%H:%M:%S'))
    ax_q.set_xlabel('Datetime [UT]  Starts: {}'.format(starttime))
    ax_q.set_ylabel('Frequency [MHz]')

    supertitle = ('{} BST Intg: {}s Dur: {}s'.format(stnid, intg, dur)
                  + ' Pointing: {},{},{}'.format(*pointing))
    plt.suptitle(supertitle)
    plt.show()


def plotsst(sstff, freqreq):
    """Plot SST data."""
    SSTdata, obsfolderinfo = dataIO.readsstfolder(sstff)
    scanrecinfo = dataIO.ScanRecInfo()
    scanrecinfo.read_scanrec(sstff)
    starttime = obsfolderinfo['datetime']
    SSTdata = numpy.array(SSTdata)
    print(SSTdata.shape)
    #freqs = modeparms.rcumode2sbfreqs(obsfolderinfo['rcumode'])
    freqs = obsfolderinfo['frequencies']
    sbreq = None
    if freqreq:
        sbreq = int(numpy.argmin(numpy.abs(freqs-freqreq)))
        show = 'persb'
    else:
        show = 'mean'
    intg = obsfolderinfo['integration']
    dur = obsfolderinfo['duration_scan']
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


def main():
    parser = argparse.ArgumentParser()

    parser.set_defaults(func=plot_bsxst)
    parser.add_argument('dataff', help="acc, bst, sst or xst filefolder")
    parser.add_argument('freq', type=float, nargs='?', default=None)

    args = parser.parse_args()
    args.dataff = os.path.normpath(args.dataff)
    args.func(args)


if __name__ == "__main__":
    main()
