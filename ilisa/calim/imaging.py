"""
Provides support direct imaging of LOFAR stand-alone data such as ACC and XST.
"""
#from __future__ import print_function

import sys
import os
import argparse

import numpy
import matplotlib.pyplot as plt
from scipy.constants import speed_of_light

import casacore.measures
import casacore.quanta.quantity
import ilisa.antennameta.antennafieldlib as antennafieldlib
from ilisa.monitorcontrol import data_io as dataIO
from ilisa.monitorcontrol.directions import _req_calsrc_proc, pointing_tuple2str,\
                                          directionterm2tuple

try:
    import dreambeam
    CANUSE_DREAMBEAM = True
except ImportError:
    CANUSE_DREAMBEAM = False
if CANUSE_DREAMBEAM:
    from dreambeam.polarimetry import convertxy2stokes, cov_lin2cir
    from dreambeam.rime.scenarios import primarybeampat

c = speed_of_light


def airydisk_radius(freq, d):
    """
    Radius of Airy disk in direction-cosine units

    Parameters
    ----------
    freq : float
        Frequency of wave in Hertz
    d : float
        Aperture diameter in meters
    
    Returns
    -------
    sintheta : float
        Radius of first null in Airy disk in direction-cosine units 
        (i.e. radius is equal to sin(theta)), if lambda
        is less than diameter. If lambda > d/1.22, then 1.0 is returned,
        which corresponds to full hemisphere field-of-view.
    """
    lambda0 = c/freq
    sintheta = 1.22*lambda0/d
    if sintheta > 1.0:
        sintheta = 1.0
    return sintheta


def calc_uvw(obstime, phaseref, stn_pos, stn_antpos):
    """
    Calculate UVW coords from datetime, phaseref and station & antenna positions

    Parameters
    ----------
    obstime : datetime
        Datetime of observation.
    phaseref : tuple
        Phase reference direction given by (lon, lat, ref),
        where lon, lat are floats for longitude, latitude in radians and ref is
        the reference frame str.
    stn_pos : array_like
        Vector of ITRF X,Y,Z coordinates in meters.
    stn_antpos : array_like
        Matrix over antenna versus ITRF X,Y,Z postions relative stn_pos.
    
    Returns
    -------
    uvw_xyz : array_like
        UVW coordinates in meters.
    """
    (pntRA, pntDEC, pntref) = phaseref
    pos_ITRF_X = str(stn_pos[0, 0])+'m'
    pos_ITRF_Y = str(stn_pos[1, 0])+'m'
    pos_ITRF_Z = str(stn_pos[2, 0])+'m'
    # Set up casacore measures object.
    obsme = casacore.measures.measures()
    where = obsme.position("ITRF", pos_ITRF_X, pos_ITRF_Y, pos_ITRF_Z)
    what = obsme.direction(pntref, str(pntRA)+"rad", str(pntDEC)+"rad")
    obsme.doframe(where)
    obsme.doframe(what)
    # Set up baselines of the array
    nrant = stn_antpos.shape[0]
    bls = []
    for antnr in range(nrant):
        bl = obsme.baseline("ITRF",
                            *[str(comp)+'m' for comp in 
                              numpy.asarray(stn_antpos[antnr, :]).squeeze()])
        bls.append(bl)
    uvw_xyz = numpy.zeros((nrant,3))
    # for obstime
    when = obsme.epoch("UTC", obstime.isoformat('T'))
    obsme.doframe(when)
    for antnr in range(nrant):
        uvw_xyz[antnr,:] = numpy.asarray(
                        obsme.to_uvw(bls[antnr])["xyz"].get_value('m'))
    return uvw_xyz


def phaseref_xstpol(xstpol, UVWxyz, freq):
    """
    Phase up polarized visibilities stack to U,V-align them at frequency 
    """
    lambda0 = c/freq
    phasefactors = numpy.exp(-2.0j*numpy.pi*UVWxyz[:,2]/lambda0)
    PP = numpy.einsum('i,k->ik', phasefactors, numpy.conj(phasefactors))
    xstpupol = numpy.array(
           [[PP*xstpol[0, 0, ...].squeeze(), PP*xstpol[0, 1, ...].squeeze()],
            [PP*xstpol[1, 0, ...].squeeze(), PP*xstpol[1, 1, ...].squeeze()]])
    return xstpupol


def phaseref_accpol(accpol, sbobstimes, freqs, stnPos, antpos, pointing):
    """
    Phase up spectral cube of visibilities (ACC) to pointing direction
    """
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
          *[str(comp)+'m' for comp in numpy.asarray(antpos[antnr, :]).squeeze()
           ])
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


def phasedup_vis(vis, srcname, t, freq, polrep, stn_pos, stn_antpos):
    """
    Phase up visibiliies
    """
    # Phase center on src
    dir_src = directionterm2tuple(srcname)
    uvw_src = calc_uvw(t, dir_src, stn_pos, stn_antpos)
    vis_pu = phaseref_xstpol(vis, uvw_src, freq)
    return vis_pu


def beamformed_image(xstpol, stn2Dcoord, freq, use_autocorr=True,
                     lmsize=2.0, polrep='linear', fluxperbeam=True):
    """
    Beamformed image XSTpol data.

    Parameters
    ----------
    xstpol : array
        The crosslet statistics data. Should have format:
        xstpol[polport1, polport2, elemnr1, elemnr2] where polport1 and
        polport2 are the two polarization ports, e.g. X and Y, and elemnr1 and
        elemnr2 are two elements of the interferometer array configuration.
    stn2Dcoord : array
        The 2D array configuration matrix.
    freq : float
        The frequency of the the data in Hz.
    use_autocorr : bool
        Whether or not to include the autocorrelations.
    lmsize : float
        Size of image in (lm) direction-cosine units. Default 2.0 means allsky.
    polrep : str
        Requested type of representation for the polarimetric data.
        Can be 'linear' (default), 'circular', or 'stokes'.
        (If dreamBeam package not accessible only 'linear' is possible)
    fluxperbeam : bool
        If True, then the flux is per beam, else the flux is per sterradian.

    Returns
    -------
    (skyimag_0, skyimag_1, skyimag_2, skyimag_3) : tuple
        Polarimetric image maps.
        If polrep='Stokes' (and dreamBeam package accessible) then the elements
        correspond to Stokes I,Q,U,V.
        If  polrep='XY' then the elements correspond to XX,XY,YX,YY.
    ll : array
        The l direction cosine of the image.
    mm : array
        The m direction cosine of the image.

    """
    if not use_autocorr:
        # Set Autocorrelations to zero:
        xstpol = numpy.copy(xstpol)  # Copy since diagonal will be rewritten
        for indi in range(2):
            for indj in range(2):
                numpy.fill_diagonal(xstpol[indi, indj, :, :], 0.0)
    posU, posV = stn2Dcoord[0, :].squeeze(), stn2Dcoord[1, :].squeeze()
    lambda0 = c / freq
    k = 2 * numpy.pi / lambda0
    lmext = lmsize/2.0
    nrpix = 101
    l, m = numpy.linspace(-lmext, lmext, nrpix), numpy.linspace(-lmext, lmext,
                                                                nrpix)
    ll, mm = numpy.meshgrid(l, m)
    bf = numpy.exp(-1.j*k*(numpy.einsum('ij,k->ijk', ll, posU)
                           + numpy.einsum('ij,k->ijk', mm, posV)))
    bfbf = numpy.einsum('ijk,ijl->ijkl', bf, numpy.conj(bf))
    skyimag_xx = numpy.einsum('ijkl,kl->ij', bfbf, xstpol[0, 0, ...].squeeze())
    skyimag_xy = numpy.einsum('ijkl,kl->ij', bfbf, xstpol[0, 1, ...].squeeze())
    skyimag_yx = numpy.einsum('ijkl,kl->ij', bfbf, xstpol[1, 0, ...].squeeze())
    skyimag_yy = numpy.einsum('ijkl,kl->ij', bfbf, xstpol[1, 1, ...].squeeze())
    if not fluxperbeam:
        ll2mm2 = ll**2+mm**2
        beyond_horizon = ll2mm2 >= 1.0
        nn = numpy.sqrt(1-ll2mm2.astype('complex'))
        # Weight values beyond horizon to zero
        nn[beyond_horizon] = 0.0
        skyimag_xx = skyimag_xx * nn
        skyimag_xy = skyimag_xy * nn
        skyimag_yx = skyimag_yx * nn
        skyimag_yy = skyimag_yy * nn
    (skyimag_0, skyimag_1, skyimag_2, skyimag_3) = (None, None, None, None)
    if not CANUSE_DREAMBEAM or polrep == 'linear':
        (skyimag_0, skyimag_1, skyimag_2, skyimag_3) =\
            (skyimag_xx, skyimag_xy, skyimag_yx, skyimag_yy)
    elif polrep == 'stokes':
        skyimag_si, skyimag_sq, skyimag_su, skyimag_sv = convertxy2stokes(
                            skyimag_xx, skyimag_xy, skyimag_yx, skyimag_yy)
        (skyimag_0, skyimag_1, skyimag_2, skyimag_3) =\
                            (skyimag_si, skyimag_sq, skyimag_su, skyimag_sv)
    elif polrep == 'circular':
        skyimag_circ = cov_lin2cir([[skyimag_xx, skyimag_xy],
                                    [skyimag_yx, skyimag_yy]])
        skyimag_0 = skyimag_circ[0][0]
        skyimag_1 = skyimag_circ[0][1]
        skyimag_2 = skyimag_circ[1][0]
        skyimag_3 = skyimag_circ[1][1]
    return (skyimag_0, skyimag_1, skyimag_2, skyimag_3), ll, mm


def nearfield_grd_image(cvcobj, filestep, cubeslice, use_autocorr=False):
    """
    Make a nearfield image along the ground from Stokes I visibility.
    (Useful for RFI).
    """
    freq = cvcobj.freqset[filestep][cubeslice]
    cvcpol_lin = dataIO.cvc2polrep(cvcobj[filestep], crlpolrep='lin')
    vis_S0 = cvcpol_lin[0, 0, cubeslice, ...] + cvcpol_lin[1, 1, cubeslice, ...]
    stn_antpos = cvcobj.stn_antpos
    if not use_autocorr:
        numpy.fill_diagonal(vis_S0[: ,:], 0.0)
    stn2dcoord = numpy.asarray(stn_antpos @ cvcobj.stn_rot)
    pos_u, pos_v = stn2dcoord[: ,0].squeeze(), stn2dcoord[: ,1].squeeze()
    lambda0 = c / freq
    k = 2 * numpy.pi / lambda0
    r_ext = 100.0
    nrpix = 2*101
    x = numpy.linspace(-r_ext, r_ext, nrpix)
    y = numpy.linspace(-r_ext, r_ext, nrpix)
    xx, yy = numpy.meshgrid(x,y)
    xx1 = xx[..., numpy.newaxis]
    yy1 = yy[..., numpy.newaxis]
    rvec = numpy.array([xx1 - pos_u, yy1 - pos_v])
    r = numpy.linalg.norm(rvec, axis=0)
    bf = numpy.exp(-1.j*k*r)
    bfbf = numpy.einsum('ijk,ijl->ijkl', bf, numpy.conj(bf))
    nfhimage = numpy.einsum('ijkl,kl->ij', bfbf, vis_S0)
    blankimage = numpy.zeros(nfhimage.shape)
    nfhimage = numpy.real(nfhimage)
    nfimages = (nfhimage, blankimage, blankimage, blankimage)
    return xx, yy, nfimages


def cvc_image(cvcobj, filestep, cubeslice, req_calsrc=None, pbcor=False,
              fluxperbeam=True, polrep='stokes'):
    """
    Image CVC object using beamformed synthesis

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
    pbcor : bool
        Perform primary beam correction or not.
    polrep : str
        Polarization representation to use for image.

    Returns
    -------
    ll : array
        The l-direction cosine map.
    mm : array
        The m-direction cosine map.
    images : tuple
        Tuple of polarized image maps.
    t : datetime
        Observation time.
    freq : float
        Observation frequency.
    phaseref : tuple
        Direction of phase reference used for imaging.
    """
    t = cvcobj.samptimeset[filestep][cubeslice]
    freq = cvcobj.freqset[filestep][cubeslice]

    pointingstr = cvcobj.scanrecinfo.get_pointingstr()
    stn_pos = cvcobj.stn_pos
    stn_antpos = cvcobj.stn_antpos

    cvcpol_lin = dataIO.cvc2polrep(cvcobj[filestep], crlpolrep='lin')

    allsky = cvcobj.scanrecinfo.get_allsky()
    phaseref = _req_calsrc_proc(req_calsrc, allsky, pointingstr)

    # Select a visibility snapshot
    cvpol_lin = cvcpol_lin[:, :, cubeslice, ...].squeeze()

    # Calculate UVW coords
    UVWxyz = calc_uvw(t, phaseref, stn_pos, stn_antpos)

    # Phase up visibilities
    cvpu_lin = phaseref_xstpol(cvpol_lin, UVWxyz, freq)

    # Determine FoV and image lm size
    bandarr = cvcobj.scanrecinfo.get_bandarr()
    lmsize = 2.0
    if not allsky:
        d = antennafieldlib.ELEMENT_DIAMETER[bandarr]
        fov = 2*airydisk_radius(freq, d)
        lmsize = 1.0*fov

    # Make image on phased up visibilities
    imgs_lin, ll, mm = beamformed_image(
        cvpu_lin, UVWxyz.T, freq, use_autocorr=False, lmsize=lmsize,
        polrep='linear', fluxperbeam=fluxperbeam)

    # Potentially apply primary beam correction 
    if pbcor and CANUSE_DREAMBEAM:
        # Get dreambeam jones:
        pointing = (float(phaseref[0]), float(phaseref[1]), 'STN')
        stnid = cvcobj.scanrecinfo.get_stnid()
        jonesfld, _stnbasis, _j2000basis = primarybeampat(
            'LOFAR', stnid, bandarr, 'Hamaker', freq, pointing=pointing,
            obstime=t, lmgrid=(ll, mm))
        ijones = numpy.linalg.inv(jonesfld)
        bri_ant = numpy.array([[imgs_lin[0], imgs_lin[1]],
                               [imgs_lin[2], imgs_lin[3]]])

        bri_ant = numpy.moveaxis(numpy.moveaxis(bri_ant, 0, -1), 0, -1)

        ijonesH = numpy.conj(numpy.swapaxes(ijones, -1, -2))
        bri_xy_iau = numpy.matmul(numpy.matmul(ijones, bri_ant), ijonesH)
        imgs_lin = (bri_xy_iau[:, :, 0, 0], bri_xy_iau[:, :, 0, 1],
                    bri_xy_iau[:, :, 1, 0], bri_xy_iau[:, :, 1, 1])
    
    # Convert to requested polarization representation
    if polrep == 'stokes' and CANUSE_DREAMBEAM:
        images = convertxy2stokes(imgs_lin[0], imgs_lin[1], imgs_lin[2],
                                  imgs_lin[3])
    elif polrep == 'circular' and CANUSE_DREAMBEAM:
        imgpolmat_lin = numpy.array([[imgs_lin[0], imgs_lin[1]],
                                     [imgs_lin[2], imgs_lin[3]]])
        imgpolmat = cov_lin2cir(imgpolmat_lin)
        images = (imgpolmat[0][0], imgpolmat[0][1], imgpolmat[1][0],
                  imgpolmat[1][1])
    else:
        # polrep == 'linear'
        images = imgs_lin

    return images, ll, mm, phaseref


def rm_redundant_bls(cvc, rmconjbl=True, use_autocorr=False):
    """Remove redundant baselines from covariance matrices.
    Assumes baseline indices are in last components."""
    nrelems = cvc.shape[-1]
    # Total number of baselines incl. autocorr and conjugate baselines.
    nrbaselinestot = nrelems**2
    if rmconjbl:
        # Remove conjugate baselines
        for idx_i in range(1, nrelems):
            for idx_j in range(idx_i):
                cvc[..., idx_i, idx_j] = 0.0
        nrbaselinestot -= nrelems*(nrelems-1)/2
    if not use_autocorr:
        # Do not use the autocorrelations (for all pol combos
        # i.e. for XX, YY, XY and YX)
        for idx in range(nrelems):
            cvc[..., idx, idx] = 0.0
        nrbaselinestot -= nrelems
    return cvc, nrbaselinestot


# Conversion between datatypes
def xst2bst(xst, obstime, pointing, stn_pos, stn_antpos, freq):
    """Convert xst data to bst data"""
    xst, nrbaselinestot = rm_redundant_bls(xst)
    UVWxyz = calc_uvw(obstime, pointing, stn_pos, stn_antpos)
    xstpu = phaseref_xstpol(xst, UVWxyz, freq)
    bstXX = numpy.sum(xstpu[0, 0, ...].squeeze(), axis=(0, 1))/nrbaselinestot
    bstXY = numpy.sum(xstpu[0, 1, ...].squeeze(), axis=(0, 1))/nrbaselinestot
    bstYY = numpy.sum(xstpu[1, 1, ...].squeeze(), axis=(0, 1))/nrbaselinestot
    return bstXX, bstXY, bstYY


def accpol2bst(accpol, sbobstimes, freqs, stn_pos, stn_antpos, pointing,
               use_autocorr=False):
    """Convert a polarized spectral cube of visibilities (ACC order by two X,Y
    indices) to polarized brightness towards pointing direction. The output is
    a set of 2 real and 1 complex powers (this is like beamlet statics, bst,
    data but also has the complex, cross-hand polarization components, which
    the bst does not have)."""
    accpol, nrbaselinestot = rm_redundant_bls(accpol)
    # Phase up ACC towards pointing direction
    accpu = phaseref_accpol(accpol, sbobstimes, freqs, stn_pos, stn_antpos,
                            pointing)
    # Sum up phased up ACCs per pol component over all baselines (Previously
    #  average)
    # Note that this sum is also over conjugate baselines, so factor 2 more
    bstXX = numpy.sum(numpy.real(accpu[0, 0, ...].squeeze()), axis=(1, 2)
                      )/nrbaselinestot
    bstXY = numpy.sum(accpu[0, 1, ...].squeeze(), axis=(1, 2))/nrbaselinestot
    bstYY = numpy.sum(numpy.real(accpu[1, 1, ...].squeeze()), axis=(1, 2)
                      )/nrbaselinestot
    return bstXX, bstXY, bstYY


def plotskyimage(ll, mm, skyimages, polrep, t, freq, stnid, integration,
                 phaseref=None, calibrated=None, pbcor=None, maskhrz=True,
                 fluxperbeam=True):
    """
    Generic plot of images of Stokes components from sky map.
    
    Parameters
    ----------
    ll : array_like
        Direction cosine, x-aligned image coordinate.
    mm : array_like
        Direction cosine, y-aligned image coordinate.
    skyimages : array_like
        Polarized sky images.
    polrep : str
        Polarization representation for skyimages.
        Can be 'stokes', 'linear', 'circular' or 'S0'.
    t : datetime
        UT date-time of image.
    freq : float
        Frequency of image in Hz.
    stnid : str
        ID of station from which image was taken.
    phaseref : (lon, lat, ref)
        Phase reference of image given as a tuple.
    integration : float
        Image exposure in seconds.
    calibrated : boolean
        Whether visibilities were calibrated or not.
    pbcor : boolean
        Primary beam correction applied or not, default None.
    maskhrz : boolean
        Mask horizon, default True.
    fluxperbeam : boolean
        Normalize flux values to be in units of flux per beam.
        Default True. If False, flux is in units of flux per steradian (s.r.).
    """

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
    hrzrgn = (ll**2+mm**2>=1.0)
    pbeamfld = (ll**2+mm**2<0.75**2)  # Primary beam field (used for color
                                    # scaling, thus choosen slightly larger
                                    # than 0.5)

    xlabel = 'Easting (dir. cos) []'
    ylabel = 'Northing (dir. cos) []'
    if fluxperbeam:
        fluxnrmlabel = 'beam'
    else:
        fluxnrmlabel = 's.r.'

    def plotcomp(compmap, compname, pos):
        if maskhrz:
            # Avoid horizon:
            compmap = numpy.ma.masked_where(hrzrgn, compmap)
        plt.subplot(2, 2, pos+1)
        vmax = numpy.amax(compmap[pbeamfld])
        vmin = numpy.amin(compmap[pbeamfld])
        plt.imshow(compmap, origin='lower', extent=[lmin, lmax, mmin, mmax],
                   interpolation='none', cmap=plt.get_cmap("jet"),
                   vmax=vmax,vmin=vmin)
        plt.gca().invert_xaxis()
        if pos == 2 or pos == 3:
            plt.xlabel(xlabel)
        else:
            plt.gca().get_xaxis().set_visible(False)
        if pos == 0 or pos == 2:
            plt.ylabel(ylabel)
        else:
            plt.gca().get_yaxis().set_visible(False)
        plt.colorbar(label="flux/"+fluxnrmlabel)
        if polrep == 'stokes':
            polrepstr = 'Stokes'
        else:
            polrepstr = 'Component'
        plt.title('{} {}'.format(polrepstr, compname))

    if polrep == 'stokes':
        # Stokes I
        sI = skyimages[0]

        # Stokes Q
        sQ = skyimages[1]
        sq = sQ / sI

        # Stokes U
        sU = skyimages[2]
        su = sU / sI

        # Stokes V
        sV = skyimages[3]
        sv = sV / sI

        s0, s0lbl = sI, 'I'
        if numpy.amin(sI) > 0.0:
            s1, s1lbl = sq, 'q'
            s2, s2lbl = su, 'u'
            s3, s3lbl = sv, 'v'
            print('Normalized')
        else:
            # If we have negative Stokes I values
            # then relative Stokes are not meaningful
            s1, s1lbl = sQ, 'Q'
            s2, s2lbl = sU, 'U'
            s3, s3lbl = sV, 'V'
        plotcomp(s0, s0lbl, 0)
        plotcomp(s1, s1lbl, 1)
        plotcomp(s2, s2lbl, 2)
        plotcomp(s3, s3lbl, 3)

    elif polrep == 'linear':
        plotcomp(numpy.real(skyimages[0]), 'XX*', 0)
        plotcomp(numpy.real(skyimages[1]), 'Re(XY*)', 1)
        plotcomp(numpy.imag(skyimages[2]), 'Im(YX*)', 2)
        plotcomp(numpy.real(skyimages[3]), 'YY*', 3)
    elif polrep == 'S0':
        #print(skyimages.shape)
        plt.imshow(skyimages[0], origin='lower',
                   extent=[lmin, lmax, mmin, mmax], interpolation='none',
                   cmap=plt.get_cmap("jet")) #, vmax=vmax, vmin=vmin)
        plt.gca().invert_xaxis()
    if calibrated:
        caltag = 'Cal'
    else:
        caltag = 'Raw'
    if not phaseref:
        phaseref = ('', '', '')
    plt.suptitle(
        """Sky Image: PhaseRef={} @ {} MHz,
        Station {}, int={}s, UT={}, PbCor={}, {}
        """.format(pointing_tuple2str(phaseref), freq/1e6, stnid, integration,
                   t, pbcor, caltag), fontsize=8)
    plt.show()


def pntsrc_hmsph(*pntsrcs, imsize=101):
    """Generate point sources on a hemisphere.
    """
    lmext = 1.0
    (l, m) = (numpy.linspace(-lmext, lmext, imsize),
              numpy.linspace(-lmext, lmext, imsize))
    ll, mm = numpy.meshgrid(l, m)
    img_S0 = numpy.zeros(ll.shape)
    img_S1 = numpy.zeros(ll.shape)
    img_S2 = numpy.zeros(ll.shape)
    img_S3 = numpy.zeros(ll.shape)
    for pntsrc in pntsrcs:
        lidx = numpy.argmin(numpy.abs(l-pntsrc[0]))
        midx = numpy.argmin(numpy.abs(m-pntsrc[1]))
        img_S0[lidx, midx] = pntsrc[2]
    return ll, mm, (img_S0, img_S1, img_S2, img_S3)


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
        intgs = len(cvcobj.samptimeset[fileidx])
        for tidx in range(args.sampnr, intgs):
            t = cvcobj.samptimeset[fileidx][tidx]
            freq = cvcobj.freqset[fileidx][tidx]
            skyimages, ll, mm, phaseref = \
                cvc_image(cvcobj, fileidx, tidx, args.phaseref,
                                  polrep=polrep, pbcor=args.correctpb,
                                  fluxperbeam=fluxperbeam)
            plotskyimage(ll, mm, skyimages, polrep, t, freq, stnid,
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
        intgs = len(cvcobj.samptimeset[fileidx])
        for tidx in range(args.sampnr, intgs):
            xx, yy, nfimages = nearfield_grd_image(cvcobj, fileidx,
                                                           tidx)
            t = cvcobj.samptimeset[fileidx][tidx]
            freq = cvcobj.freqset[fileidx][tidx]
            plotskyimage(xx, yy, nfimages, polrep, t, freq, stnid,
                                 integration, maskhrz=False)


def main_cli():
    """\
    CLI to image LOFAR station data
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-n', '--filenr', type=int, default=0)
    parser.add_argument('-s', '--sampnr', type=int, default=0)

    subparsers = parser.add_subparsers(help='sub-command help')

    parser_image = subparsers.add_parser('bf', help='beamform image')
    parser_image.set_defaults(func=image)
    parser_image.add_argument('dataff', help="acc or xst filefolder")
    parser_image.add_argument('-p', '--phaseref', type=str, default=None)
    parser_image.add_argument('-c', '--correctpb',
                              help="Correct for primary beam",
                              action="store_true")
    parser_image.add_argument('-f', '--fluxpersterradian',
                              help="Normalize flux per sterradian",
                              action="store_true")

    parser_image = subparsers.add_parser('nf', help='nearfield image')
    parser_image.set_defaults(func=nfimage)
    parser_image.add_argument('dataff', help="acc or xst filefolder")

    args = parser.parse_args()
    args.dataff = os.path.normpath(args.dataff)
    args.func(args)

if __name__ == "__main__":
    main_cli()
