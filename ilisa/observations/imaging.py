"""Provides support direct imaging of LOFAR stand-alone data such as ACC and
XST."""
from __future__ import print_function
import sys
import casacore.measures
import casacore.quanta.quantity
import numpy
from scipy.constants import speed_of_light
import matplotlib.pyplot as plt
import ilisa.antennameta.antennafieldlib as antennafieldlib
import ilisa.antennameta.calibrationtables as calibrationtables
import ilisa.observations.directions
import ilisa.observations.modeparms as modeparms
import ilisa.observations.dataIO as dataIO
try:
    from dreambeam.polarimetry import convertxy2stokes
    canuse_stokes = True
except ImportError:
    canuse_stokes = False
try:
    from dreambeam.rime.scenarios import primarybeampat
    canuse_dreambeam = True
except ImportError:
    canuse_dreambeam = False

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


def phaseref_accpol(accpol, sbobstimes, freqs, stnPos, antpos, pointing):
    """Phase up a spectral cube of visibilities (ACC) to pointing direction."""
    accphasedup = numpy.zeros(accpol.shape, dtype=complex)
    (pntRA, pntDEC, pntref) = pointing
    pos_ITRF_X = str(stnPos[0, 0])+'m'
    pos_ITRF_Y = str(stnPos[1, 0])+'m'
    pos_ITRF_Z = str(stnPos[2, 0])+'m'
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
        bl.append(obsme.baseline("ITRF",
          *[str(comp)+'m' for comp in numpy.asarray(antpos[antnr, :]).squeeze()])
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
        phasefactors = numpy.exp(-2.0j*numpy.pi*UVWxyz[:, 2]/lambda0)
        PP = numpy.einsum('i,j->ij', phasefactors, numpy.conj(phasefactors))
        accphasedup[0, 0, sb, :, :] = PP*accpol[0, 0, sb, :, :]
        accphasedup[0, 1, sb, :, :] = PP*accpol[0, 1, sb, :, :]
        accphasedup[1, 0, sb, :, :] = PP*accpol[1, 0, sb, :, :]
        accphasedup[1, 1, sb, :, :] = PP*accpol[1, 1, sb, :, :]
    return accphasedup


def xst2skyim_stn2Dcoord(xstpol, stn2Dcoord, freq, include_autocorr=True, allsky=False,
                         polrep_req='Stokes'):
    """Image XSTpol data.

    Parameters
    ----------

    xstpol : array
        The crosslet statistics data.
    stn2Dcoord: array
        The 2D array configuration matrix.
    freq : float
        The frequency of the the data in Hz.
    include_autocorr : bool
        Whether or not to include the autocorrelations.
    allsky : bool
        Should an allsky image be produced?
    polrep_req : str
        Requested type of representation for the polarimetric data.
        Can be 'Stokes' or 'XY'. (If dreamBeam package not accessible only 'XY' is
        possible)

    Returns
    -------

    polrep : str
        The polarization representation of the image tuple.
        Can be 'Stokes' or 'XY'.
    (skyimag_0, skyimag_1, skyimag_2, skyimag_3): tuple
        Polarimetric image maps.
        If polrep='Stokes' (and dreamBeam package accessible) then the elements
        correspond to Stokes I,Q,U,V.
        If  polrep='XY' then the elements correspond to XX,XY,YX,YY.
    ll : array
        The l direction cosine of the image.
    mm : array
        The m direction cosine of the image.

    """
    if not include_autocorr:
        for indi in range(2):
            for indj in range(2):
                numpy.fill_diagonal(xstpol[indi,indj,:,:],0.0)
    posU, posV = stn2Dcoord[0,:].squeeze(), stn2Dcoord[1,:].squeeze()
    lambda0 = c / freq
    k = 2 * numpy.pi / lambda0
    if not allsky:
        lmext = fov(freq)
    else:
        lmext = 1.0
    nrpix = 101
    l, m = numpy.linspace(-lmext, lmext, nrpix), numpy.linspace(-lmext, lmext, nrpix)
    ll, mm = numpy.meshgrid(l,m)
    bf = numpy.exp(-1.j*k*(numpy.einsum('ij,k->ijk', ll, posU)
                           + numpy.einsum('ij,k->ijk', mm, posV)))
    bfbf = numpy.einsum('ijk,ijl->ijkl', bf, numpy.conj(bf))
    skyimag_xx = numpy.einsum('ijkl,kl->ij', bfbf, xstpol[0, 0, ...].squeeze())
    skyimag_xy = numpy.einsum('ijkl,kl->ij', bfbf, xstpol[0, 1, ...].squeeze())
    skyimag_yx = numpy.einsum('ijkl,kl->ij', bfbf, xstpol[1, 0, ...].squeeze())
    skyimag_yy = numpy.einsum('ijkl,kl->ij', bfbf, xstpol[1, 1, ...].squeeze())
    if polrep_req == 'Stokes' and canuse_stokes:
        skyimag_si, skyimag_sq, skyimag_su, skyimag_sv = convertxy2stokes(
            skyimag_xx, skyimag_xy, skyimag_yx, skyimag_yy)
        polrep = 'Stokes'
        return  polrep, (skyimag_si, skyimag_sq, skyimag_su, skyimag_sv), ll, mm
    else:
        polrep = 'XY'
        return polrep, (skyimag_xx, skyimag_xy, skyimag_yx, skyimag_yy), ll, mm


def nearfield_grd_image(vis_S0, stn2Dcoord, freq, include_autocorr=False):
    """Make a nearfield image along the ground from Stokes I visibility.
    (Useful for RFI)."""
    if not include_autocorr:
        numpy.fill_diagonal(vis_S0[:,:],0.0)
    stn2Dcoord = numpy.squeeze(numpy.asarray( stn2Dcoord))
    posU, posV = stn2Dcoord[:,0].squeeze(), stn2Dcoord[:,1].squeeze()
    lambda0 = c / freq
    k = 2 * numpy.pi / lambda0
    r_ext = 500.0
    print(r_ext)
    nrpix = 101
    x, y = numpy.linspace(-r_ext, r_ext, nrpix), numpy.linspace(-r_ext, r_ext, nrpix)
    xx, yy = numpy.meshgrid(x,y)
    xx1 = xx[...,numpy.newaxis]
    yy1 = yy[...,numpy.newaxis]
    rvec = numpy.array([xx1 - posU, yy1 - posV])
    r = numpy.linalg.norm(rvec, axis=0)
    #plt.imshow(r[:,:,50])
    #plt.colorbar()
    #plt.show()
    bf = numpy.exp(-1.j*k*r)
    bfbf = numpy.einsum('ijk,ijl->ijkl', bf, numpy.conj(bf))
    nfhimage = numpy.einsum('ijkl,kl->ij', bfbf, vis_S0)
    blankimage = numpy.zeros(nfhimage.shape)
    return (nfhimage, blankimage, blankimage, blankimage), xx, yy


def cvcimage(cvcobj, filestep, cubeslice, req_calsrc=None, docalibrate=True,
             pbcor=False):
    """Image a CVC object.

    Parameters
    ----------

    cvcobj : CVCfiles()
        The covariance cube files object containing visibility cubes.
    filestep : int
        The file index to select.
    cubeslice : int
        The cube index in file.
    req_calsrc : str
        The requested sky source.
    docalibrate : bool
        Perform calibration or not.
    pbcor : bool
        Perform primary beam correction or not.

    Returns
    -------

    ll : array
        The l-direction cosine map.
    mm : array
        The m-direction cosine map.
    images : tuple
        Tuple of polarized image maps.
    polrep : str
        Polarization representation of the images tuple.
    t : datetime
        Observation time.
    freq : float
        Observation frequency.
    stnid : str
        Station id.
    phaseref : tuple
        Direction of phase reference used for imaging.

    """
    t = cvcobj.samptimeset[filestep][cubeslice]
    freq = cvcobj.freqset[filestep][cubeslice]
    cvcdata_unc = cvcobj.getdata(filestep)

    stnid = cvcobj.stnsesinfo.get_stnid()
    bandarr = cvcobj.stnsesinfo.get_bandarr()
    band = cvcobj.stnsesinfo.get_band()
    rcumode = cvcobj.stnsesinfo.get_rcumode()
    pointingstr = cvcobj.stnsesinfo.get_pointingstr()

    # Get/Compute ant positions
    stnPos, stnRot, antpos, stnIntilePos \
                            = antennafieldlib.getArrayBandParams(stnid, bandarr)
    septon = cvcobj.stnsesinfo.is_septon()
    if septon:
        elmap = cvcobj.stnsesinfo.get_septon_elmap()
        for tile, elem in enumerate(elmap):
            antpos[tile] = antpos[tile] + stnIntilePos[elem]

    # stn2Dcoord = stnRot.T * antpos.T
    # Apply calibration
    datatype = cvcobj.stnsesinfo.get_datatype()
    if datatype == 'acc':
        cvcdata, caltabhead = calibrationtables.calibrateACC(cvcdata_unc, rcumode, stnid,
                                                             t, docalibrate)
    else:
        sb, nz = modeparms.freq2sb(freq)
        cvcdata, caltabhead = calibrationtables.calibrateXST(cvcdata_unc, sb, rcumode,
                                                             stnid, t, docalibrate)
    cvpol = dataIO.cvc2cvpol(cvcdata)

    # Determine if allsky FoV
    if band == '10_90' or band == '30_90' or septon:
        allsky = True
    else:
        allsky = False

    # Determine phaseref
    if req_calsrc is not None:
        pntstr = ilisa.observations.directions.std_pointings(req_calsrc)
    elif allsky:
        pntstr = ilisa.observations.directions.std_pointings('Z')
    else:
        pntstr = pointingstr
    phaseref = pntstr.split(',')

    skyimage = True
    if skyimage:
        # Phase up visibilities
        cvpu, UVWxyz = phaseref_xstpol(cvpol[:,:,cubeslice,...].squeeze(),
                                       t, freq, stnPos, antpos, phaseref)

        # Make image on phased up visibilities
        polrep, images, ll, mm = xst2skyim_stn2Dcoord(
            cvpu, UVWxyz.T, freq, include_autocorr=False, allsky=allsky, polrep_req='XY')
        if pbcor and canuse_dreambeam:
            # Get dreambeam jones:
            pointing = (float(phaseref[0]), float(phaseref[1]), 'STN')
            jonesfld, stnbasis, j2000basis = primarybeampat(
                'LOFAR', stnid, bandarr, 'Hamaker', freq, pointing=pointing, obstime=t,
                lmgrid=(ll,mm))
            ijones = numpy.linalg.inv(jonesfld)
            bri_ant = numpy.array([[images[0], images[1]],
                                   [images[2], images[3]]])

            bri_ant = numpy.moveaxis(numpy.moveaxis(bri_ant, 0, -1), 0, -1)

            ijonesH = numpy.conj(numpy.swapaxes(ijones, -1, -2))
            bri_xy_iau = numpy.matmul(numpy.matmul(ijones, bri_ant), ijonesH)
            images = (bri_xy_iau[:, :, 0, 0], bri_xy_iau[:, :, 0, 1],
                      bri_xy_iau[:, :, 1, 0], bri_xy_iau[:, :, 1, 1])
        if canuse_stokes and polrep == 'XY':
            images = convertxy2stokes(images[0], images[1], images[2], images[3])
            polrep = 'Stokes'
    else:
        vis_S0 = cvpol[0,0,cubeslice,...].squeeze() + cvpol[0,0,cubeslice,...].squeeze()
        nfhimages, ll, mm = nearfield_grd_image(vis_S0, antpos, freq,
                                                include_autocorr=True)
        polrep = 'S0'
        images = numpy.real(nfhimages)

    return ll, mm, images, polrep, t, freq, stnid, phaseref


# Conversion between datatypes
def xst2bst(xst, obstime, freq, stnPos, antpos, pointing):
    """Convert xst data to bst data"""
    xstpu, UVWxyz = phaseref_xstpol(xst, obstime, freq, stnPos, antpos, pointing)
    bstXX = numpy.sum(xstpu[0,0,...].squeeze(), axis=(0,1))
    bstXY = numpy.sum(xstpu[0,1,...].squeeze(), axis=(0,1))
    #bstYX = numpy.sum(xstpu[1,0,...].squeeze(), axis=(0,1))
    bstYY = numpy.sum(xstpu[1,1,...].squeeze(), axis=(0,1))
    return bstXX, bstXY, bstYY #, bstYX


def accpol2bst(accpol, sbobstimes, freqs, stnPos, antpos, pointing, use_autocorr = False):
    """Convert a polarized spectral cube of visibilities (ACC order by two X,Y indices)
    to polarized brightness towards pointing direction. The output is a set of 2 real and
    1 complex powers (this is like beamlet statics, bst, data but also has the complex,
    cross-hand polarization components, which the bst does not have)."""
    nrelems = accpol.shape[-1]
    nrbaselinestot = nrelems**2  # i.e. total number of baselines incl. autocorr and
                                 # conjugate baselines.
    if not use_autocorr:
        # Do not use the autocorrelations (for all pol combos i.e. for XX, YY, XY and YX)
        for idx in range(nrelems):
            accpol[...,idx,idx]= 0.0
        nrbaselinestot -= nrelems
    # Phase up ACC towards pointing direction
    accpu = phaseref_accpol(accpol, sbobstimes, freqs, stnPos, antpos, pointing)
    # Sum up phased up ACCs per pol component over all baselines (Previously average)
    # Note that this sum is also over conjugate baselines, so factor 2 more
    bstXX = numpy.sum(numpy.real(accpu[0,0,...].squeeze()), axis=(1,2))
    bstXY = numpy.sum(accpu[0,1,...].squeeze(), axis=(1,2))
    # bstYX is redundant: it is conjugate of bstXY.
    # bstYX = numpy.sum(accpu[1,0,...].squeeze(), axis=(1,2))
    bstYY = numpy.sum(numpy.real(accpu[1,1,...].squeeze()), axis=(1,2))
    return bstXX, bstXY, bstYY  # , bstYX


def plotskyimage(ll, mm, skyimages, polrep, t, freq, stnid, phaseref, integration,
                 pbcor=None):
    """Generic plot of images of Stokes components from sky map."""

    # Compute extent
    dl = ll[0, 1] - ll[0, 0]
    lmin = ll[0, 0] - dl / 2.0
    lmax = ll[-1, -1] + dl / 2.0
    dm = mm[1, 0] - mm[0, 0]
    mmin = mm[0, 0] - dm / 2.0
    mmax = mm[-1, -1] + dm / 2.0
    domaincheck = False
    if domaincheck:
        fig_dt, (ax_dt_ll, ax_dt_mm) = plt.subplots(2)
        im_ll = ax_dt_ll.imshow(ll)
        ax_dt_ll.set_title('ll')
        fig_dt.colorbar(im_ll, ax=ax_dt_ll)
        im_mm = ax_dt_mm.imshow(mm)
        ax_dt_mm.set_title('mm')
        fig_dt.colorbar(im_mm, ax=ax_dt_mm)
        plt.show()
    hrzrgn = (numpy.sqrt(ll**2+mm**2)>.99)

    #norm=colors.LogNorm()
    norm=None
    xlabel = 'Easting (dir. cos) []'
    ylabel = 'Northing (dir. cos) []'

    def plotcomp(compmap, compname, pos):
        # Avoid horizon:
        compmap = numpy.ma.masked_where(hrzrgn, compmap)
        plt.subplot(2, 2, pos)
        if pos == 1:
            vmax = None
            vmin = None
        else:
            vmax = 1
            vmin = -1
        plt.imshow(compmap, origin='lower', extent=[lmin, lmax, mmin, mmax],
                   interpolation='none') #, vmax=vmax, vmin=vmin)
        plt.gca().invert_xaxis()
        if pos == 3 or pos == 4:
            plt.xlabel(xlabel)
        if pos == 1 or pos == 3:
            plt.ylabel(ylabel)
        plt.colorbar()
        if polrep == 'Stokes':
            polrepstr = 'Stokes'
        else:
            polrepstr = 'Component'
        plt.title('{} {}'.format(polrepstr, compname))

    if polrep == 'Stokes':
        # Stokes I
        sI = skyimages[0]  # numpy.ma.masked_where(skyimages[0]<0.0, skyimages[0])
        plotcomp(sI, 'I', 1)

        # Stokes V
        sV = skyimages[3]
        sv = sV / sI
        plotcomp(sV, 'V', 3)

        # Stokes Q
        sQ = skyimages[1]
        sq = sQ / sI
        plotcomp(sQ, 'Q', 4)

        # Stokes U
        sU = skyimages[2]
        su = sU / sI
        plotcomp(sU, 'U', 2)
    elif polrep == 'XY':
        plotcomp(numpy.real(skyimages[0]), 'XX*', 1)
        plotcomp(numpy.real(skyimages[1]), 'Re(XY*)', 2)
        plotcomp(numpy.imag(skyimages[2]), 'Im(YX*)', 3)
        plotcomp(numpy.real(skyimages[3]), 'YY*', 4)
    plt.suptitle('Sky Image: PhaseRef={} @ {} MHz\nStation {}, int={}s, UT={}, PbCor={}'\
        .format(','.join(phaseref), freq/1e6, stnid, integration, t, pbcor))
    plt.tight_layout(rect=[0, 0.0, 1, 0.95])
    plt.show()
