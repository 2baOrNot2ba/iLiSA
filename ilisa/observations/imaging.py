"""Provides support direct imaging of LOFAR stand-alone data such as ACC and
XST."""
from __future__ import print_function
import sys
import datetime
import casacore.measures
import casacore.quanta.quantity
import numpy
from scipy.constants import speed_of_light
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import ilisa.antennameta.antennafieldlib as antennafieldlib
import ilisa.antennameta.calibrationtables as calibrationtables
import ilisa.observations.modeparms
import ilisa.observations.stationcontrol as stationcontrol
import ilisa.observations.observing as observing
import ilisa.observations.dataIO as dataIO

c = speed_of_light


def fov(freq):
    if freq < 100e6:
        return 1.0
    else:
        return 0.6*100e6/freq/2.0


def phaseref_xstpol(xstpol, obstime0, freq, stnPos, antpos, pointing):
    """ """
    (pntRA, pntDEC, pntref) = pointing
    pos_ITRF_X = str(stnPos[0,0])+'m'
    pos_ITRF_Y = str(stnPos[1,0])+'m'
    pos_ITRF_Z = str(stnPos[2,0])+'m'
    obsme = casacore.measures.measures()
    where = obsme.position("ITRF", pos_ITRF_X, pos_ITRF_Y, pos_ITRF_Z)
    what = obsme.direction(pntref, pntRA+"rad", pntDEC+"rad")
    when = obsme.epoch("UTC", obstime0.isoformat('T'))
    obsme.doframe(where)
    obsme.doframe(when)
    obsme.doframe(what)
    nrant = antpos.shape[0]
    UVWxyz = numpy.zeros((nrant,3))
    for idx in range(nrant):
        bl = obsme.baseline("ITRF",
            *[str(comp)+'m' for comp in numpy.asarray(antpos[idx,:]).squeeze()])
        UVWxyz[idx,:] = numpy.asarray(obsme.to_uvw(bl)["xyz"].get_value('m'))
    lambda0 = c/freq
    phasefactors = numpy.exp(-2.0j*numpy.pi*UVWxyz[:,2]/lambda0)
    PP = numpy.einsum('i,k->ik',phasefactors, numpy.conj(phasefactors))
    xstpupol = \
      numpy.array([[PP*xstpol[0,0,...].squeeze(), PP*xstpol[0,1,...].squeeze()],
                  [PP*xstpol[1,0,...].squeeze(), PP*xstpol[1,1,...].squeeze()]])
    return xstpupol, UVWxyz


# def phaseref_accs(acc, sbobstimes, freqs, stnPos, antpos, pointing):
#     """ """
#     accphasedup = numpy.zeros(acc.shape, dtype=complex)
#     print("sb", end='')
#     for sb in range(acc.shape[0]):
#         print(sb, end='')
#         xst = acc[sb,...].squeeze()
#         xstphasedup, UVWxyz = phaseref_xst(xst, sbobstimes[sb], freqs[sb],
#                                        stnPos, antpos, pointing)
#         accphasedup[sb,...] = xstphasedup
#     return accphasedup


def phaseref_accpol(accpol, sbobstimes, freqs, stnPos, antpos, pointing):
    """Phase up a spectral cube of visibilities (ACC) to pointing direction."""
    accphasedup = numpy.zeros(accpol.shape, dtype=complex)
    (pntRA, pntDEC, pntref) = pointing
    pos_ITRF_X = str(stnPos[0,0])+'m'
    pos_ITRF_Y = str(stnPos[1,0])+'m'
    pos_ITRF_Z = str(stnPos[2,0])+'m'
    # Set up casacore measures object.
    obsme = casacore.measures.measures()
    where = obsme.position("ITRF", pos_ITRF_X, pos_ITRF_Y, pos_ITRF_Z)
    what = obsme.direction(pntref, pntRA+"rad", pntDEC+"rad")
    obsme.doframe(where)
    obsme.doframe(what)
    # Set up baselines of the array
    nrant = antpos.shape[0]
    bl = []
    for antnr in range(nrant):
        bl.append( obsme.baseline("ITRF",
          *[str(comp)+'m' for comp in numpy.asarray(antpos[antnr,:]).squeeze()])
        )
    # Step through subbands & phase up autocovariance matrix to point direction
    for sb in range(len(sbobstimes)):
        sys.stdout.write("\rsb: {}/511".format(sb))
        sys.stdout.flush()
        when = obsme.epoch("UTC", sbobstimes[sb].isoformat('T'))
        obsme.doframe(when)
        UVWxyz = numpy.zeros((nrant,3))
        for antnr in range(nrant):
            UVWxyz[antnr,:] = numpy.asarray(
                                  obsme.to_uvw(bl[antnr])["xyz"].get_value('m'))
        lambda0 = c/freqs[sb]
        # Compute phase factors. w-direction (component 2) is towards pointing
        phasefactors = numpy.exp(-2.0j*numpy.pi*UVWxyz[:,2]/lambda0)
        PP = numpy.einsum('i,j->ij',phasefactors, numpy.conj(phasefactors))
        accphasedup[0,0,sb,:,:] = PP*accpol[0,0,sb,:,:]
        accphasedup[0,1,sb,:,:] = PP*accpol[0,1,sb,:,:]
        accphasedup[1,0,sb,:,:] = PP*accpol[1,0,sb,:,:]
        accphasedup[1,1,sb,:,:] = PP*accpol[1,1,sb,:,:]
    print()
    return accphasedup


def xst2skyim_stn2Dcoord(xstpol, stn2Dcoord, freq, include_autocorr=True, allsky=False):

    if not include_autocorr:
        for indi in range(2):
            for indj in range(2):
                numpy.fill_diagonal(xstpol[indi,indj,:,:],0.0)
    posU, posV = stn2Dcoord[0,:].squeeze(), stn2Dcoord[1,:].squeeze()
    lambda0 = c / freq
    k = 2 * numpy.pi / lambda0;
    if not allsky:
        lmext = fov(freq)
    else:
        lmext = 1.0
    nrpix = 101
    l, m = numpy.linspace(-lmext,lmext,nrpix), numpy.linspace(-lmext,lmext,
                                                                          nrpix)
    ll, mm = numpy.meshgrid(l,m)
    #skymap = numpy.zeros((len(l), len(m)))
    bf = numpy.exp(-1.j*k*(numpy.einsum('ij,k->ijk', ll, posU)
                        +numpy.einsum('ij,k->ijk', mm, posV)))
    bfbf = numpy.einsum('ijk,ijl->ijkl',bf,numpy.conj(bf))
    skyimageXX = numpy.einsum('ijkl,kl->ij', bfbf, xstpol[0,0,...].squeeze())
    skyimageXY = numpy.einsum('ijkl,kl->ij', bfbf, xstpol[0,1,...].squeeze())
    skyimageYX = numpy.einsum('ijkl,kl->ij', bfbf, xstpol[1,0,...].squeeze())
    skyimageYY = numpy.einsum('ijkl,kl->ij', bfbf, xstpol[1,1,...].squeeze())
    # TODO: make sure stokes are IAU coordinates.
    skyimageSI = numpy.real(skyimageXX+skyimageYY)
    skyimageSQ = numpy.real(skyimageXX-skyimageYY)
    skyimageSU = numpy.real(skyimageXY+skyimageYX)
    skyimageSV = numpy.imag(skyimageXY-skyimageYX)
    return (skyimageSI, skyimageSQ, skyimageSU, skyimageSV), ll, mm


def cvcimage(cvcpath, cubeslice, req_calsrc, docalibrate = True):
    """Image a CVC file"""
    cvcpath = cvcpath.rstrip('/')
    if '_acc' in cvcpath:
        cvctype = 'acc'
    else:
        cvctype = 'xst'
    cvcobj = dataIO.CVCfiles(cvcpath)
    obsfolderinfo = cvcobj.getobsfolderinfo()
    if type(obsfolderinfo) is list and type(obsfolderinfo[0]) is dataIO.ObsInfo:
        obsinfos = obsfolderinfo
        obsinfo = obsinfos[0]
        t0, rcumode = obsinfo.starttime, obsinfo.beamctl_cmd['rcumode']
    else:
        t0, rcumode = obsfolderinfo['datetime'], obsfolderinfo['rcumode']
    if cvctype == 'acc':
        cvcdata_unc = cvcobj.getdata(0)  # FIXME allow imaging other than index 0
        ts = cvcobj.samptimes[0]
        calsrc = obsfolderinfo['calsrc']
        stnid = obsfolderinfo['stnid']
        sb, nz = stationcontrol.freq2sb(float(cubeslice))
        cubeslice = sb
        t = ts[cubeslice]
    else:
        cvcdata_unc = cvcobj.getdata(0)
        starttime = obsinfo.starttime
        stnid = obsinfo.stnid
        beamctl_cmd = obsinfo.beamctl_cmd
        sb = int(obsinfo.rspctl_cmd['xcsubband'])
        t = t0 + datetime.timedelta(seconds=float(cubeslice))
        cubeslice = int(cubeslice)
    pntstr = ilisa.observations.modeparms.stdPointings(req_calsrc)
    pointing = pntstr.split(',')
    bandarr = stationcontrol.rcumode2antset(rcumode).split("_")[0]
    stnPos, stnRot, antpos, stnIntilePos \
                            = antennafieldlib.getArrayBandParams(stnid, bandarr)
    try:
        obsinfo.septonconf
    except AttributeError:
        obsinfo.septonconf = None
    if obsinfo.septonconf is not None:
        elmap = ilisa.observations.modeparms.str2elementMap2(obsinfo.septonconf)
        for tile, elem in enumerate(elmap):
            antpos[tile] = antpos[tile] + stnIntilePos[elem]
        allsky = True
    else:
        allsky = False
    freqs = stationcontrol.rcumode2sbfreqs(rcumode)
    # stn2Dcoord = stnRot.T * antpos.T
    # Apply calibration
    if cvctype == 'acc':
        cvcdata, caltabhead = calibrationtables.calibrateACC(cvcdata_unc,
                                                rcumode, stnid, t0, docalibrate)
    else:
        cvcdata, caltabhead = calibrationtables.calibrateXST(cvcdata_unc, sb,
                                                rcumode, stnid, t0, docalibrate)
    cvpol = dataIO.cvc2cvpol(cvcdata)
    cvpu, UVWxyz = phaseref_xstpol(cvpol[:,:,cubeslice,...].squeeze(),
                                   t, freqs[sb], stnPos, antpos, pointing)
    skyimages, ll, mm = xst2skyim_stn2Dcoord(cvpu, UVWxyz.T, freqs[sb], True, allsky)
    plotskyimage(ll, mm, skyimages, t, freqs[sb], stnid)

# Conversion between datatypes
def xst2bst(xst, obstime, freq, stnPos, antpos, pointing):
    """Convert xst data to bst data"""
    xstpu, UVWxyz = phaseref_xstpol(xst, obstime, freq, stnPos, antpos,
                                                                       pointing)
    bstXX = numpy.sum(xstpu[0,0,...].squeeze(), axis=(0,1))
    bstXY = numpy.sum(xstpu[0,1,...].squeeze(), axis=(0,1))
    #bstYX = numpy.sum(xstpu[1,0,...].squeeze(), axis=(0,1))
    bstYY = numpy.sum(xstpu[1,1,...].squeeze(), axis=(0,1))
    return bstXX, bstXY, bstYY #, bstYX

def accpol2bst(accpol, sbobstimes, freqs, stnPos, antpos, pointing):
    """Convert spectral cube of visibilities (ACC) to polarized brightness
    towards pointing direction. The output is a set of 2 real and 1 complex
    powers (this is like beamlet statics, bst, data but also has the complex,
    cross-hand polarization components, which the bst does not have)."""
    # Phase up ACC towards pointing direction
    accpu = phaseref_accpol(accpol, sbobstimes, freqs, stnPos, antpos,
                                                                       pointing)
    # Sum up phased up ACCs per polarization component over all baselines
    bstXX = numpy.sum(numpy.real(accpu[0,0,...].squeeze()), axis=(1,2))
    bstXY = numpy.sum(accpu[0,1,...].squeeze(), axis=(1,2))
    # bstYX is redunant: it is conjugate of bstXY.
    #bstYX = numpy.sum(accpu[1,0,...].squeeze(), axis=(1,2))
    bstYY = numpy.sum(numpy.real(accpu[1,1,...].squeeze()), axis=(1,2))
    return bstXX, bstXY, bstYY #, bstYX


def plotskyimage(ll, mm, skyimages, t, freq, stnid):
    """Generic plot of images."""
    plt.subplot(2,2,1)
    #norm=colors.LogNorm()
    norm=None
    plt.pcolormesh(ll, mm, skyimages[0], norm=norm)
    plt.gca().invert_xaxis()
    plt.xlabel('E<-W direction cosine')
    plt.ylabel('S->N direction cosine')
    plt.colorbar()
    plt.title('Stokes I')
    plt.subplot(2,2,2)
    plt.pcolormesh(ll, mm, skyimages[1]/skyimages[0])
    plt.gca().invert_xaxis()
    plt.colorbar()
    plt.title('Stokes q')
    plt.subplot(2,2,3)
    plt.pcolormesh(ll, mm, skyimages[2]/skyimages[0])
    plt.gca().invert_xaxis()
    plt.colorbar()
    plt.title('Stokes u')
    plt.subplot(2,2,4)
    plt.pcolormesh(ll, mm, skyimages[3]/skyimages[0]) #, vmin=-1, vmax=1)
    plt.gca().invert_xaxis()
    plt.colorbar()
    plt.title('Stokes v')
    plt.suptitle('Allsky: {} @ {} MHz, UT={}'.format(stnid, freq/1e6, t))
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    plt.show()
