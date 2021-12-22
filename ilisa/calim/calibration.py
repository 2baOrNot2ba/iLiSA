import shutil
import numpy
from scipy.constants import speed_of_light as c

from casacore.measures import measures
import ilisa.calim
from ilisa.antennameta import calibrationtables as calibrationtables
from ilisa.operations import data_io as dataIO, modeparms as modeparms


def applycaltab_cvc(cvcunc, caltab, sb=None):
    """Apply a caltable to CVC data.

    Note
    ----
    Formula used is:
        G_ijk = g_ij*g^*_ik (no sum over i)
        V'_ijk = G_s0jk*V_ijk  (no sum over j,k & s0 is explicitly given)

    (Note that since function was designed for simplicity, it determines
    whether the CVC is an ACC (with all subbands) or XST (only one given
    subband) based on if no subband is given and its first index has size 512.
    There is a very small chance for the user to make a mistake by not
    setting 'sb' and the data happens to be a 512 samples XST.)

    Parameter
    ----------
    cvcunc : array [:,nr_rcus_0,nr_rcus_1]
        Uncalibrated CVC array, where 1st index is time. nr_rcus is
        usually 192.
    sb : int
        Subband in which the XST data was taken.
    caltab: array [512,nr_rcus]
        Calibration table array to apply to xstunc.

    Returns
    -------
    cvccal: array [:,nr_rcus,nr_rcus]
        Calibrated XST array.
    """
    nrsbs = cvcunc.shape[0]
    if not sb and nrsbs != 512:
        # Cannot assume its an ACC
        raise ValueError("Must give sb for XST data.")
    gg = numpy.einsum('ij,ik->ijk', caltab, numpy.conj(caltab))
    if not sb and nrsbs == 512:
        # Assume it's an ACC
        g_apply = gg
    else:
        # It's an XST
        g_apply = gg[sb, :, :]
    cvccal = g_apply*cvcunc
    return cvccal


def applycal_cvcfolder(cvcpath, caltabpath):
    """Apply a calibration table file to a CVC folder.

    This creates a copy of the folder pointed to by cvcpath renamed with
    a '_cal_' before the ldat suffix. Then it applies the calibration table
    contained in the caltab file pointed to by caltabpath, to the CVC dataset
    and copies the caltab file used in the calibrated CVC folder.

    Parameters
    ----------
    cvcpath: str
        Path to CVC folder
    caltabpath: str
        Path to caltab file
    """
    try:
        caltab, header = calibrationtables.readcaltab(caltabpath)
    except:
        raise
    ldat_type = dataIO.datafolder_type(cvcpath)
    if ldat_type != "acc" and ldat_type != "xst":
        raise ValueError("Not CVC data.")
    # Copy CVC folder within parent folder and add the tag "_cal_" in the name
    # right before ldat_type suffix:
    spltpath = cvcpath.split("_")
    cvccalpath = "_".join(spltpath[:-1]) + "_cal_" + spltpath[-1]
    shutil.copytree(cvcpath, cvccalpath)
    # Read in cvcobj:
    cvcobj_cal = dataIO.CVCfiles(cvccalpath)
    nrfiles = cvcobj_cal.getnrfiles()
    # Loop over files in CVC folder:
    for filestep in range(nrfiles):
        if ldat_type == "xst":
            freq = cvcobj_cal.freqset[filestep][0]  # Freq const. over xst file
            sb, nz = modeparms.freq2sb(freq)
        else:
            sb = None  # Because this signals ACC data
        # Get actual covariance cubes:
        cvcdata_unc = cvcobj_cal[filestep]
        # Apply calibration
        cvcdata = ilisa.calim.calibration.applycaltab_cvc(cvcdata_unc, caltab,
                                                          sb)
        # Replace uncalibrated data with calibrated:
        cvcobj_cal[filestep] = cvcdata
    # Note in ScanRecInfo about calibrating this dataset:
    cvcobj_cal.scanrecinfo.set_postcalibration(caltabpath, cvccalpath)


def apply_polgains(gainspol, vispol):
    """
    Apply polarized biscalar gains to polarized visibility

    Parameters
    ----------
    gainspol : array_like
        Polarized biscalar gains. g_pq has shape (2,N) where N is number of
        dual-polarized elements.
    vispol : array_like
        Polarized visibility array. Its shape is (2,2,N,N)

    Returns
    -------
    vispol_gapp : array_like
        Gain applied polarized visibility. Formula is g_ik*vis_ijkl*conj(g_jl),
        where none of the indices are summed over.
        It will have same shape as vis_pq.
    
    Examples
    --------
    >>> from ilisa.calim.calibration import apply_polgains
    >>> import numpy
    >>> n=3
    >>> vispqtrue = numpy.ones((2,2,n,n))
    >>> gpq=numpy.ones((2,n))

    Double the P pol channel and zeros the Q channel:
    >>> gpq[0,:]=2.0
    >>> gpq[1,:]=0.0
    Apply this to visibiities
    >>> vis=apply_polgains(gpq, vispqtrue)
    >>> vis
    array([[[[4., 4., 4.],
             [4., 4., 4.],
             [4., 4., 4.]],

            [[0., 0., 0.],
             [0., 0., 0.],
             [0., 0., 0.]]],


           [[[0., 0., 0.],
             [0., 0., 0.],
             [0., 0., 0.]],

            [[0., 0., 0.],
             [0., 0., 0.],
             [0., 0., 0.]]]])
    """
    #vispq = numpy.copy(vispq)
    #    (2,N)^(2,N) swap 1,2 => (2,2,N,N) 
    gg = numpy.tensordot(gainspol, numpy.conj(gainspol),0).swapaxes(1,2)
    vispol_gapp = gg * vispol
    return vispol_gapp


def stefcal(r, m, niter=100):
    """
    Compute solution to the quadratic matrix equation using the StefCal
    algorithm [stefs]_.

    Returns a solution, up to a phase factor, to the equation
    .. math:: g r g^H = m
    where g is a column vector and r, m are Hermitian matrices. 

    Parameters
    ----------
    r: array [nr_rcus, nr_rcus]
        The measured covariance matrix.
    m: array [nr_rcus, nr_rcus]
        The model covariance matrix.
    niter: int, optional
        Number of iterations. Default equal to 100.
    
    Returns
    -------
    g: array [nr_rcus]
        Complex gain solution vector, unnormalized.
    
    Examples
    --------
    >>> from ilisa.calim.calibration import stefcal
    >>> import numpy as np
    >>> dim = 3
    >>> gtrue = np.random.randn(dim)+1j*np.random.randn(dim)
    >>> m = np.ones((dim, dim))
    >>> r = gtrue[:, np.newaxis]*m*np.conj(gtrue)
    >>> g = stefcal(r,m)

    Normalize w.r.t to first component and check result:

    >>> np.allclose(g/g[0],gtrue/gtrue[0])
    True

    Notes
    -----
    As pointed out in [Bhatnagar]_, StefCal is essentially equivalent to the
    'antsol' algorithm used in CASA and APES.
    
    References
    ----------
    .. [stefs]: S. Salvini, & S. J. Wijnholds,
       "Fast gain calibration in radio astronomy using alternating direction
       implicit methods: Analysis and applications", A&A 571 A97, 2014.
       DOI: 10.1051/0004-6361/201424487
    .. [Bhatnagar]: S. Bhatnagar,
       "StefCal vs. Classical antsol: A critique", 2013.
       URL: http://www.aoc.nrao.edu/~sbhatnag/misc/stefcal.pdf
    """
    dim = r.shape[0]
    g_prev = numpy.ones((dim,), dtype=complex)
    g_curr = numpy.zeros((dim,), dtype=complex)
    for indi in range(2, niter):
        for indp in range(dim):
            z = g_prev*m[:, indp]
            g_curr[indp] = (numpy.inner(numpy.conj(r[:, indp]), z)
                            /numpy.inner(numpy.conj(z), z))
        if indi % 2 == 0:
            g_curr = (g_curr+g_prev)/2.0
        g_prev = g_curr
    g = g_curr
    return g


def gain_cal_bs_lin(vis_pol_src):
    """
    Apply stefcal to XX and YY visibilities

    Parameters
    ----------
    vis_pol_obs : (2, 2, N, N) array
        Observsered (measured) polarized visibilities
    
    Returns
    -------
    g_bs_lin : (2, N) array
        Gain solutions for X and Y channels
    """
    vis_mod_xx = numpy.ones(vis_pol_src[0, 0, ...].shape)
    #vis_mod_xy = numpy.zeros(vis_pol_obs[0, 1, ...].shape)
    #vis_mod_yx = numpy.zeros(vis_pol_obs[1, 0, ...].shape)
    vis_mod_yy = numpy.ones(vis_pol_src[1, 1, ...].shape)
    g_xx = stefcal(vis_pol_src[0, 0, ...], vis_mod_xx)
    #g_xy = stefcal(vis_pol_obs[0, 1, ...], vis_mod_xy)
    #g_yx = stefcal(vis_pol_obs[1, 0, ...], vis_mod_yx)
    g_yy = stefcal(vis_pol_src[1, 1, ...], vis_mod_yy)
    g_bs_lin = numpy.array([g_xx, g_yy])
    return g_bs_lin


def vcz(ll, mm, skyimage, freq, ant_pos, imag_is_fd=False):
    """\
    Compute visibility from image via the van Cittert-Zernicke relation

    Parameters
    ----------
    ll : array_like
        Direction cosine, x-aligned image coordinate.
    mm : array_like
        Direction cosine, y-aligned image coordinate.
    skyimage : array_like
        Total flux image of sky.
    uv: array_like
        U,V vectors corresponding to 2D array configuration.
    imag_is_fd : bool
        Input image is flux distribution map (fd) rather than flux density
        distribution (fdd) map.

    Returns
    -------
    vis : array_like
        Visibility corresponding to input skyimage via vCZ relation.
    """
    pos_x = ant_pos[:, 0].squeeze()
    pos_y = ant_pos[:, 1].squeeze()
    pos_z = ant_pos[:, 2].squeeze()

    nr_ants = pos_x.shape[0]
    vis = numpy.zeros((nr_ants, nr_ants), dtype=complex)
    phas = 0.0
    k = 2 * numpy.pi * freq / c
    lm_r2 = (ll**2+mm**2).astype(numpy.complex)
    nn = numpy.sqrt(1-lm_r2)
    nn[lm_r2 >= 1.0] = 1.0
    numpy.imag(nn)
    if imag_is_fd:
        fdd = skyimage / nn
    else:
        fdd = skyimage
    for ant_i in range(nr_ants):
        for ant_j in range(ant_i, nr_ants):
            u = pos_x[ant_i] - pos_x[ant_j]
            v = pos_y[ant_i] - pos_y[ant_j]
            w = pos_z[ant_i] - pos_z[ant_j]
            vis[ant_i, ant_j] = numpy.sum(
                fdd * numpy.exp(+1.0j * k * (ll * u + mm * v + (nn-1) * w)))
    do_conj = True
    if do_conj:
        for ant_i in range(nr_ants):
            for ant_j in range(0, ant_i):
                vis[ant_i, ant_j] = numpy.conj(vis[ant_j, ant_i])
    return vis

from ilisa.calim import imaging
from ilisa.calim import skymodels # import gdskymodel
import matplotlib.pyplot as plt
def gsmcal(dataff, filenr, sampnr, fluxpersterradian):
    ccm = measures()
    gs_model = 'LFSM'
    imsize = 200
    fluxperbeam = True
    l = numpy.linspace(-1, 1, imsize)
    m = numpy.linspace(-1, 1, imsize)
    ll, mm = numpy.meshgrid(l, m)
    #normcolor = colors.LogNorm()
    # The code below is almost cut-n-pasted from imaging.image()
    polrep = 'stokes'
    lofar_datatype = dataIO.datafolder_type(dataff)
    fluxperbeam = not fluxpersterradian
    if lofar_datatype != 'acc' and lofar_datatype != 'xst':
        raise RuntimeError("Datafolder '{}'\n not ACC or XST type data."
                           .format(dataff))
    cvcobj = dataIO.CVCfiles(dataff)
    calibrated = False
    if cvcobj.scanrecinfo.calibrationfile:
        calibrated = True
    stnid = cvcobj.scanrecinfo.get_stnid()
    lon, lat, h = imaging.ITRF2lonlat(cvcobj.stn_pos[0, 0],
                                      cvcobj.stn_pos[1, 0],
                                      cvcobj.stn_pos[2, 0])
    stn_pos_x, stn_pos_y, stn_pos_z = cvcobj.stn_pos[0, 0], cvcobj.stn_pos[1, 0], \
                                      cvcobj.stn_pos[2, 0]
    ccm.doframe(ccm.position('ITRF', str(stn_pos_x) + 'm', str(stn_pos_y) + 'm',
                             str(stn_pos_z) + 'm'))
    for fileidx in range(filenr, cvcobj.getnrfiles()):
        integration = cvcobj.scanrecinfo.get_integration()
        intgs = len(cvcobj.samptimeset[fileidx])
        for tidx in range(sampnr, intgs):
            t = cvcobj.samptimeset[fileidx][tidx]
            freq = cvcobj.freqset[fileidx][tidx]
            img = skymodels.globaldiffuseskymodel(t, (lon, lat, h), freq,
                                                  gs_model=gs_model,
                                                  imsize=imsize)
            ccm.doframe(ccm.epoch('UTC', t.isoformat('T')))
            phaseref_ccm = ccm.measure(ccm.direction('AZEL', '0.0rad', str(numpy.deg2rad(90))+'rad'), 'J2000')
            phaseref = (phaseref_ccm['m0']['value'],  phaseref_ccm['m1']['value'], phaseref_ccm['refer'])
            uvw_sl = imaging.calc_uvw(t, phaseref, cvcobj.stn_pos, cvcobj.stn_antpos)
            vis_mod = vcz(ll, mm, img, freq, uvw_sl, imag_is_fd=not(fluxperbeam))
            cvcpol_lin = dataIO.cvc2polrep(cvcobj[fileidx], crlpolrep='lin')
            cvpol_lin = cvcpol_lin[:, :, tidx, ...].squeeze()
            cvpol_x = cvpol_lin[0,0,...].squeeze()
            cvpol_y = cvpol_lin[1,1,...].squeeze()
            vis_meas_xx = cvpol_x
            vis_meas_yy = cvpol_y
            g_xx = stefcal(vis_meas_xx, vis_mod/2.0)
            g_yy = stefcal(vis_meas_yy, vis_mod/2.0)
            inv_g_xx = 1/g_xx
            inv_g_yy = 1/g_yy
            vis_cal_xx = inv_g_xx[:,numpy.newaxis]*vis_meas_xx*numpy.conj(inv_g_xx)
            vis_cal_yy = inv_g_yy[:,numpy.newaxis]*vis_meas_yy*numpy.conj(inv_g_yy)
            vis_cal_I = vis_cal_xx + vis_cal_yy
            vis_resid_I = vis_cal_I - vis_mod
            xstpol = numpy.array([[vis_cal_I, numpy.zeros_like(vis_mod)],
                                  [numpy.zeros_like(vis_mod), vis_resid_I]])
            skyimages, _l, _m = imaging.beamformed_image(xstpol, uvw_sl.T, freq,
                                                         use_autocorr=True,
                                                         lmsize=2.0,
                                                         nrpix=imsize,
                                                         polrep='linear',
                                                         fluxperbeam=fluxperbeam
                                                         )
            imaging.plotskyimage(_l, _m, skyimages, 'linear', t, freq, stnid,
                                 integration, phaseref, calibrated, pbcor=False,
                                 maskhrz=False, fluxperbeam=fluxperbeam)


def gsmcal_cli():
    parser = argparse.ArgumentParser()
    parser.add_argument('dataff',
                        help="""Path to CVC folder""")
    args = parser.parse_args()
    args.fluxpersterradian = False
    args.filenr = 180
    args.sampnr = 0
    gains = gsmcal(args.dataff, args.filenr, args.sampnr,
                   args.fluxpersterradian)


import argparse


def main_cli():
    """
    Apply a calibration file to ACC or XST data folder.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('cvcpath',
                        help="""Path to CVC folder""")
    parser.add_argument('caltabpath', help="""Path to caltab file""")
    args = parser.parse_args()

    applycal_cvcfolder(args.cvcpath, args.caltabpath)


if __name__ == "__main__":
    #main_cli()
    gsmcal_cli()
