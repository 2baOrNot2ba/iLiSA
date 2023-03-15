"""
Module dealing with aspects of radio telescope beams

Beams here is in the most general sense, and could be synthesized or primary.
"""
import argparse

import numpy
import pkg_resources
from matplotlib import pyplot as plt

import ilisa.calim.flagging
from ilisa.calim.flagging import Flags
from ilisa.antennameta import antennafieldlib as antennafieldlib
from ilisa.operations.modeparms import rcumode2sbfreqs
from . import SPEED_OF_LIGHT as c
from . import visibilities as vsb
import ilisa.calim.imaging


def beamformed_pattern(stn2Dcoord, freq, vis_flags):
    """
    Compute beamformed beam flux pattern

    Parameters
    ----------
    stn2Dcoord: array
        Positions of array elements.
    freq: float
        Frequency.
    flagged_vis: dict of bool matrices
        Keyed set of visibility flag-matrices.

    Returns
    -------
    ll, mm: array
        Direction cosine grids.
    skyimages: tuple
        Polarized sky images.
    """
    nrants = stn2Dcoord.shape[0]
    vis_xx = numpy.ones((nrants, nrants))
    vis_yy = vis_xx
    vis_xy = numpy.zeros_like(vis_xx)
    vis_yx = vis_xy
    vis_pol = numpy.array([[vis_xx, vis_xy],[vis_yx, vis_yy]])
    vis_flags.vis = vis_pol
    vis_pol = vis_flags.apply_vispol_flags()
    skyimages, ll, mm \
        = ilisa.calim.imaging.beamformed_image(vis_pol, stn2Dcoord.T, freq, nrpix=201)
    return ll, mm, skyimages


def beam_pat_shape(ll, mm, images_bf_pat):
    """
    Assess beamformed pattern

    Parameters
    ----------
    ll, mm: array
        Direction cosine grids.
    images_bf_pat: tuple
        Polarized images of beamformed pattern.

    Returns
    -------
    major_diam: float
        Major diameter of beam ellipse
    minor_diam: float
        Minor diameter of beam FWHM.
    elltilt: float
        Angle of tilt of beam ellipse in radians.
    fov_area: float
        FoV area in direction cosine squared units.
    """
    def maskbb(crdgrd, msk):
        msk_crdgrd = msk * crdgrd
        crd_max = numpy.amax(msk_crdgrd)
        crd_min = numpy.amin(msk_crdgrd)
        return (crd_max+crd_min)/2, crd_max-crd_min

    image_0 = numpy.real(images_bf_pat[0])
    beamcenter = numpy.unravel_index(numpy.argmax(image_0), image_0.shape)
    fluxmax = image_0[beamcenter]
    halfmax = fluxmax/2.0
    beammask = numpy.zeros_like(image_0)
    beammask[image_0>halfmax]=1.0
    ll_me, mm_me = numpy.mean(beammask*ll), numpy.mean(beammask*mm)
    ll_cntr, mm_cntr = ll-ll_me, mm-mm_me
    ll2mm2 = numpy.sqrt(ll_cntr**2+mm_cntr**2)
    semimajor_idx = numpy.unravel_index(numpy.argmax(beammask*ll2mm2),
                                        image_0.shape)
    semimajor_vec = numpy.array([ll_cntr[semimajor_idx], mm_cntr[semimajor_idx]]
                                )
    if semimajor_vec[0]<0.0:
        semimajor_vec = -semimajor_vec
    semiminor_uvec = numpy.array([-semimajor_vec[1], semimajor_vec[0]])/ \
                     (numpy.linalg.norm(semimajor_vec))
    semiminor_llmm = semiminor_uvec[0]*ll_cntr+semiminor_uvec[1]*mm_cntr
    smin_cntr, smin_ext = maskbb(semiminor_llmm, beammask)
    major_diam = 2*numpy.linalg.norm(semimajor_vec)
    minor_diam = smin_ext
    elltilt = numpy.arctan2(semimajor_vec[1], semimajor_vec[0])
    fov_area = area_beamell(major_diam, minor_diam)
    return major_diam, minor_diam, elltilt, fov_area


def get_beam_shape_parms(stnid, antset, freq, flagged_vis, _use_lookuptab=None):
    """
    Return station beam field-of-view size

    Parameters
    ----------
    stnid: str
        Station ID.
    antset: str
        Band ID.
    freq: float
        Frequency.
    flagged_vis: Flags
        Visibility flags
    _use_lookuptab : bool
        Use lookup table instead of calculating FoV size.

    Returns
    -------
    major_diam: float
        Major diameter of beam ellipse
    minor_diam: float
        Minor diameter of beam FWHM.
    elltilt: float
        Angle of tilt of beam ellipse in radians.
    fov_area: float
        FoV size in direction cosine.
    """
    if  _use_lookuptab is None:
        # Decide whether to use lookup table or not (since it is not set)
        if flagged_vis is None:
            # No flags set so lookup table works
            _use_lookuptab = True
        else:
            bl_flags = flagged_vis.bl_mask
            # Check if vis_flag is not just autocorrelation flags
            ac_mask = numpy.ones_like(bl_flags, dtype=bool)
            numpy.fill_diagonal(ac_mask, False)
            if numpy.any(numpy.logical_and(bl_flags, ac_mask)):
                # Lookup table won't work since nonzero baselines are flagged
                _use_lookuptab = False
                print("Will compute beampattern")
            else:
                # Just autocorrelations flagged, so looktable applies
                _use_lookuptab = True
    if not _use_lookuptab:
        stn_pos, stn_rot, stn_antpos, stn_intilepos \
            = antennafieldlib.get_antset_params(stnid, antset)
        antpos_uv = vsb.rot2uv(stn_antpos, stn_rot)
        ll, mm, bfps = beamformed_pattern(antpos_uv, freq, flagged_vis)
        major_diam, minor_diam, elltilt, fov_area = beam_pat_shape(ll, mm, bfps)
    else:
        beamshape_path = pkg_resources.resource_filename(__name__,
                                'beamshape_' + stnid + '_' + antset + '.npy')
        try:
            beamshapes = numpy.load(beamshape_path)
        except FileNotFoundError:
            major_diam, minor_diam, elltilt, fov_area = 0., 0., 0., None
        else:
            freqidx = numpy.argmin(numpy.abs(beamshapes[:,0] - freq))
            freq_cntr, major_diam, minor_diam, elltilt, fov_area \
                = beamshapes[freqidx, :]
    return major_diam, minor_diam, elltilt, fov_area


def area_beamell(major_diam, minor_diam):
    """
    Compute beam ellipse area

    Parameters
    ----------
    major_diam: float
        Major diameter in units of direction cosine.
    minor_diam: float
        Minor diameter in units of direction cosine.

    Returns
    -------
    fov_sz_dlm: float
        FoV size in direction cosine squared.
    """
    fov_sz_dlm = major_diam*minor_diam/4*numpy.pi
    return fov_sz_dlm


def nrpixels_hint(minor_diam, lm_extent, pixsperminor=5):
    """
    Suggest the number of pixels to use for an image

    minor_diam: float
        Beam ellipse minor diameter in direction-cosine units.
    lm_extent: int
        Extent (length of edge) of image in direction-cosine units.
    pixsperminor: int
        Number of pixels desired to fit along beam ellipse minor axis.

    Returns
    -------
    nrpixels: int
        Number of pixels to use along image axis.
    """
    if minor_diam == 0.0:
        # Invalid so send back invalid nrpixels for caller to deal with.
        return 0
    nrpixels = int(numpy.ceil(pixsperminor*(lm_extent/minor_diam)))
    return nrpixels


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


def main_cli():
    """
    Beampattern CLI

    Plots or prints out the synthesized beam of a specified station, band-array,
    and frequency. It also precomputes beampattern parameters to save in a
    lookup table .npy file.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-p', '--plot', action="store_true",
                        help="Plot the results")
    parser.add_argument('-s' , '--stnid', type=str, default='')
    parser.add_argument('-a', '--antset', type=str, default='')
    parser.add_argument('-f', '--freq', type=float, default=0.0)
    parser.add_argument('-b', '--blflags', type=str, default='[]',
                        help="Baseline flag select")
    args = parser.parse_args()
    if not args.stnid:
        raise(RuntimeError, 'Choose stnid')
    if not args.antset:
        args.antset = 'LBA'  # TODO: HBA
    if not args.freq:
        if 'LBA' in args.antset:
            rcumode = 3
        freqs = rcumode2sbfreqs(rcumode)
    else:
        freqs = [args.freq]
    stn_pos, stn_rot, stn_antpos, stn_intilepos \
        = antennafieldlib.get_antset_params(args.stnid, args.antset)
    antpos_uv = vsb.rot2uv(stn_antpos, stn_rot)

    beamshapes = []
    print('Freq Major Minor Tilt FoV')
    for freq in freqs:
        flagged_bls = Flags(nrelems=antpos_uv.shape[0]).set_blflagargs(args.blflags)
        ll, mm, bfps = beamformed_pattern(antpos_uv, freq, flagged_bls)
        madi, midi, tlt, fov_area = beam_pat_shape(ll, mm, bfps)
        beamshape = freq, madi, midi, numpy.rad2deg(tlt), fov_area
        print(*beamshape)
        print("Pixels over hemisphere:", nrpixels_hint(midi, 2))
        beamshapes.append(beamshape)
        if args.plot and len(freqs)==1:
            ilisa.calim.imaging.plotskyimage(ll, mm, bfps, 'linear', 0, freq,
                                             args.stnid, 0)
            plt.show()
    if len(freqs) > 1:
        outfilename = 'beamshape_' + args.stnid + '_' + args.antset + '.npy'
        numpy.save(outfilename, beamshapes)


if __name__ == "__main__":
    main_cli()
