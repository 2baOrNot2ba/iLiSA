import sys
import os
import shutil

import numpy
import numpy as np
from numpy.linalg import norm

from ilisa.calim.geodesy import ITRF2lonlat
from ilisa.antennameta import calibrationtables as calibrationtables
from ilisa.operations import data_io as data_io, modeparms as modeparms
from ilisa.operations.directions import _req_calsrc_proc
from . import visibilities, skymodels, beam, imaging


def reldiffnorm(x, y):
    """
    Calculate the relative difference norm between arguments

    Parameters
    ----------
    x, y: array
        Two arrays

    Returns
    -------
    reldif: float
        Relative difference norm
    """
    reldif = 2*norm(x - y)/norm(x + y)
    return reldif


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
    """
    Apply a calibration table file to a CVC folder

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
        caltab, header = calibrationtables.read_caltabfile(caltabpath)
    except:
        raise
    ldat_type = data_io.datafolder_type(cvcpath)
    if ldat_type != "acc" and ldat_type != "xst":
        raise ValueError("Not CVC data.")
    # Copy CVC folder within parent folder and add the tag "_cal_" in the name
    # right before ldat_type suffix:
    spltpath = cvcpath.split("_")
    cvccalpath = "_".join(spltpath[:-1]) + "_cal_" + spltpath[-1]
    shutil.copytree(cvcpath, cvccalpath)
    # Read in cvcobj:
    cvcobj_cal = data_io.CVCfiles(cvccalpath)
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
        cvcdata = applycaltab_cvc(cvcdata_unc, caltab, sb)
        # Replace uncalibrated data with calibrated:
        cvcobj_cal[filestep] = cvcdata
    # Note in ScanRecInfo about calibrating this dataset:
    cvcobj_cal.scanrecinfo.set_postcalibration(caltabpath, cvccalpath)


def apply_gains_noises(vispol, gainspol, noises=None, variant='legacy',
                       ac_vis=False):
    """
    Apply polarized biscalar gains to polarized visibility

    Parameters
    ----------
    vispol : array_like
        Polarized visibility array. Its shape is (...,2,2,N,N), for pol-indexed,
        full visibility matrix, or (...,2,N). for pol-indexed autocorrelations.
    gainspol : array_like
        Polarized biscalar gains. g_pq has shape (T,2,N) where N is number of
        dual-polarized elements.
    noises: array_like
        Antenna noise powers and has shape (T,2,N). This vector represents
        the diagonal of the noise power matrix.
    variant : {'legacy', 'inv'}
        Which gain solution variant to use. 'legacy' uses formula
            V^{est} = (1/g) * (r - n) * (1/g)^H
        while 'inv' uses
            V^{est} = g * r * g^H - n
    ac_vis : bool, default=True
        True if `vispol` array a pol-indexed autocorrelation, if False then
        `vis_pol` is a normal, full visibility.
    Returns
    -------
    vispol_gapp : array_like
        Gain applied polarized visibility. Formula is g_ik*vis_ijkl*conj(g_jl),
        where none of the indices are summed over.
        It will have same shape as vis_pq.
    
    Examples
    --------
    >>> from ilisa.calim.calibration import apply_gains_noises
    >>> import numpy
    >>> n=3
    >>> vispqtrue = numpy.ones((2,2,n,n))
    >>> gpq=numpy.ones((2,n))

    Double the P pol channel and zeros the Q channel:
    >>> gpq[0,:]=2.0
    >>> gpq[1,:]=0.0
    Apply this to visibiities
    >>> vis=apply_gains_noises(gpq, vispqtrue)
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
    if not ac_vis:
        # Full visibility
        if variant=='legacy':
            gainspol = 1.0 / gainspol
        #    (T,2,N),conj(T,2,N) => (T,2,1,N,1)*conj(T,1,2,1,N) => (T,2,2,N,N)
        gg = (gainspol[...,    :, None,    :, None] * numpy.conj(
              gainspol[..., None,    :, None,    :]))
        vispol_gapp = gg * vispol
        if noises is not None:
            # Update vis estimate `V^{est}` by removing noise estimate `n` from
            # measured vis `r`, i.e.
            #    V^{est}=(1/g)*(r-n)*(1/g)^H  'legacy'
            #    V^{est}=g*r*g^H-n            'inv'
            if variant=='legacy':
                noises *= numpy.abs(gainspol)**2
            # noises.shape == (T,2,N) so apply diag on axis 1 AND -1:
            vispol_gapp -= numpy.apply_along_axis(numpy.diag, 1,
                numpy.apply_along_axis(numpy.diag, -1, noises))
    else:
        # Auto-correlation visibility
        gg = numpy.abs(gainspol) ** 2
        if variant == 'legacy':
            vispol_gapp = 1/gg * (vispol - noises)
        else:  # 'inv'
            vispol_gapp = gg * vispol - noises
    return vispol_gapp

def apply_polgains_cvcfolder(dataff, variant='legacy'):
    """
    Apply polarimetric gains to CVC data

    Parameters
    ----------
    dataff_mod: str
        Name of file-folder
    variant: str
        Which gain solution variant to use. Can be either: 'legacy' or 'inv'.

    Returns
    -------
    cvcobj_cal: CVCfiles
        The calibrated CVCfiles object
    """
    dataff_raw, dataff_mod, dataff_cal = data_io.dataff_raw_model_cal(dataff)
    # Copy raw into cal
    shutil.copytree(dataff_raw, dataff_cal)
    # Read in cvcobj:
    cvcobj_cal = data_io.CVCfiles(dataff_cal)
    # Read in gain solutions
    gainsolfile = 'gain_solutions.npy'
    noisesolfile = 'noises_raw.npy'
    if variant=='inv':
        gainsolfile = 'gain_inv_solutions.npy'
        noisesolfile = 'noises_model.npy'
    gainsolpath = os.path.join(dataff_mod, gainsolfile)
    gainsols = numpy.load(gainsolpath)
    noisesols = numpy.load(os.path.join(dataff_mod, noisesolfile))
    nrfiles = cvcobj_cal.getnrfiles()
    # Loop over files in CVC folder:
    for filestep in range(nrfiles):
        # Get actual covariance cubes:
        cvpol_unc = visibilities.cov_flat2polidx(cvcobj_cal[filestep])
        # Apply gains to uncal vis:
        cvcdata_cal = apply_gains_noises(cvpol_unc, gainsols[filestep],
                                         noisesols[filestep],
                                         variant=variant)
        # Replace uncalibrated data with calibrated:
        cvcobj_cal[filestep] = visibilities.cov_polidx2flat(cvcdata_cal)
    # Note in ScanRecInfo about calibrating this dataset:
    cvcobj_cal.scanrecinfo.set_postcalibration(gainsolpath, dataff_cal)
    return cvcobj_cal


def scale_lin(r, m):
    r_f = numpy.ravel(r)
    m_f = numpy.ravel(m)
    o_f = numpy.ravel(numpy.eye(r.shape[0]))
    a = numpy.column_stack((r_f, o_f))
    gc = numpy.matmul(numpy.linalg.pinv(a), m_f)
    g = gc[0]
    c = gc[1]
    return g, c


def stefcal(r, m, niter=100, incl_autocor=True):
    """
    Compute solution to the quadratic matrix equation using the StefCal
    algorithm [stefs]_.

    Returns a solution, up to a phase factor, to the equation
    .. math:: g m g^H = r
    where g is a column vector and r, m are Hermitian matrices.

    This function is used for calibration and can be used in two main ways:
    Alt I: Take r to be the measured covariance (visibility) matrix,
           and m to be the model covariance matrix. In this case the function
           returned gain solutions g can be applied to the measured visibility
           matrix r to get calibrated matrix by using the formula
           .. math:: r_cal = (1/g) r (1/g)^H
    Alt II: Take m to be the measured covariance matrix and r as the desired
            cov. matrix, and for which the solutions can be applied as
            .. math:: m_cal = g m g^H
            This latter method tends to work best.

    Parameters
    ----------
    r: array [nr_rcus, nr_rcus]
        The measured covariance matrix.
    m: array [nr_rcus, nr_rcus]
        The model covariance matrix.
    niter: int, optional
        Number of iterations. Default equal to 100.
    incl_autocor: bool
        Whether to include autocorrelations in the calculations.
    
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
    For stefcal to work, m matrix needs to have nondegenerate singular values.

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
    if not incl_autocor:
        # Disassociate input r, m, so autocor can be removed only locally
        r = r.copy()
        m = m.copy()
        numpy.fill_diagonal(r, 0.0)
        numpy.fill_diagonal(m, 0.0)
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

def wals(r, m, variant='legacy', nitr=100, err_tol=1e-3, mask=None):
    """
    Weighted alternating least-square solver

    Solves argmin_{g,n}(||g*r*g^H-m-n||).

    Parameters
    ----------
    r: array
        Covariance matrix
    m: array
        Covariance matrix
    variant: str
        What variant of WALS to use: 'legacy' for legacy variant, 'inv' for
        inverse variant.
    nitr: int
        Max number of iterations
    err_tol: float
        Relative error tolerance

    Returns
    -------
    g: vec
        Gains solution
    n: vec
        Noise diagonal

    Rasies
    ------
    ValueError
        If algorithm does not converge.
    """
    def wals_legacy(r, m, nitr, err_tol):
        """WALS legacy version

        Solves: argmin_{g,n}(||r-g*m*g^H-n||)
        """
        g_prev = numpy.ones(r.shape[0], dtype=complex)
        n_prev = numpy.zeros(r.shape[0], dtype=float)
        for itr in range(nitr):
            incl_autocor = True
            if itr == 0:
                incl_autocor = True
            # Use stefcal using Alt I
            g = stefcal(r-numpy.diag(n_prev), m, incl_autocor=incl_autocor)
            gg = (g[:, None] * numpy.conj(g[None, :]))
            n = numpy.real(numpy.diag(r-gg*m).copy())
            #n_med = numpy.median(n)
            #n[numpy.argmax(n)]=n_med
            #n[numpy.argmin(n)]=n_med
            serr = numpy.abs(  numpy.vdot(numpy.linalg.pinv(g_prev[None,:]), g)
                             + numpy.vdot(numpy.linalg.pinv(n_prev[None,:]), n)
                             - 1)
            if serr < err_tol:
                pass
            g_prev = g
            n_prev = n
            print('walsleg', itr, serr, numpy.amax(n)/numpy.amin(n),
                  numpy.amax(numpy.abs(g))/numpy.amin(numpy.abs(g)))
        if itr+1 == nitr:
            pass
            #raise ValueError('Has not converged')
        return g, n

    def wals_inv(r, m, nitr, err_tol):
        """WALS inverse version
        Solves: argmin_{g,n}(||g*r*g^H-m-n||)
        """
        g_prev = numpy.ones(r.shape[0], dtype=complex)
        n_prev = numpy.zeros(r.shape[0], dtype=float)
        for itr in range(nitr):
            incl_autocor = True
            if itr == 0:
                incl_autocor = True
            # Use stefcal using Alt II
            g = stefcal(m+numpy.diag(n_prev), r, incl_autocor=incl_autocor)
            gg = (g[:, None] * numpy.conj(g[None, :]))
            n = numpy.real(numpy.diag(gg*r - m)).copy()
            #n_med = numpy.median(n)
            #n[numpy.argmax(n)]=n_med
            #n[numpy.argmin(n)]=n_med
            #err = norm(gg * r - m - numpy.diag(n)) / (norm(gg * r)+norm(m)+norm(n))
            #rerr = reldiffnorm(n, n_prev)+reldiffnorm(g, g_prev)
            serr = numpy.abs(  numpy.vdot(numpy.linalg.pinv(g_prev[None,:]), g)
                             + numpy.vdot(numpy.linalg.pinv(n_prev[None,:]), n)
                             - 1)
            if serr < err_tol:
                pass
            g_prev = g
            n_prev = n
            print('walsinv', itr, serr, numpy.amax(n)/numpy.amin(n),
                  numpy.amax(numpy.abs(g))/numpy.amin(numpy.abs(g)))
        if itr + 1 == nitr:
            pass
            #raise ValueError('Has not converged')
        return g, n

    if mask is not None:
        r = numpy.ma.array(r, mask=mask).filled(0.0)
        m = numpy.ma.array(m, mask=mask).filled(0.0)
    if variant == 'inv':
        return wals_inv(r, m, nitr=nitr, err_tol=err_tol)
    return wals_legacy(r, m, nitr=nitr, err_tol=err_tol)

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


def gainsolve(cvcobj_uncal, cvcobj_model, wals_variant='legacy', nitr=100):
    """
    Solve for gains based on uncalibrated and model data

    Parameters
    ----------
    dataff: str
        Path to CVC filefolder.
    gs_model:  str
        Name Global skymodel. Choose 'LFSM', 'GSM'.

    Returns
    -------
    gainsolutions: array
        The gain solutions.
    """
    gainsolutions = []
    noisesolutions = []
    nrfiles = cvcobj_uncal.getnrfiles()
    for fileidx in range(nrfiles):
        print('Solving for file {}/{}'.format(fileidx, nrfiles))
        intgs = len(cvcobj_uncal.samptimeset[fileidx])
        gainsol_t = []
        noise_t = []
        cvcpol_uncal = visibilities.cov_flat2polidx(cvcobj_uncal[fileidx])
        cvcpol_model = visibilities.cov_flat2polidx(cvcobj_model[fileidx])
        t0 = cvcobj_uncal.samptimeset[fileidx][0]
        refdir = _req_calsrc_proc(None, cvcobj_uncal.scanrecinfo.get_allsky(),
                                  cvcobj_uncal.scanrecinfo.get_pointingstr())
        #uvw_xyz = visibilities.calc_uvw(t0, refdir, cvcobj_uncal.stn_pos,
        #                                cvcobj_uncal.stn_antpos)
        mask = None
        # Put back autocorrelations
        #numpy.fill_diagonal(mask, False)
        for tidx in range(intgs):
            t = cvcobj_uncal.samptimeset[fileidx][tidx]
            freq = cvcobj_uncal.freqset[fileidx][tidx]
            vis_uncal = cvcpol_uncal[tidx]
            vis_model = cvcpol_model[tidx]
            print('loop', fileidx, tidx, 'xx')
            g, c = scale_lin(vis_uncal[0, 0, ...], vis_model[0, 0, ...])
            cd = c*numpy.eye(vis_uncal[0, 0, ...].shape[-1])
            # print('scal',g,c)
            g_xx, n_xx = wals(g*vis_uncal[0, 0, ...]+cd, vis_model[0, 0, ...],
                              variant=wals_variant,  nitr=nitr, err_tol=1e-2,
                              mask=mask)
            print('loop', fileidx, tidx, 'yy')
            g, c = scale_lin(vis_uncal[0, 0, ...], vis_model[0, 0, ...])
            cd = c*numpy.eye(vis_uncal[0, 0, ...].shape[-1])
            #print('scal',g,c)
            g_yy, n_yy = wals(g*vis_uncal[1, 1, ...]+cd, vis_model[1, 1, ...],
                              variant=wals_variant,  nitr=nitr, err_tol=1e-2,
                              mask=mask)
            # Use Stefcal with Alt II, so measured as 2nd arg
            #g_xx = stefcal(vis_model[0, 0, ...], vis_uncal[0, 0, ...], incl_autocor=incl_autocor)
            #g_yy = stefcal(vis_model[1, 1, ...], vis_uncal[1, 1, ...], incl_autocor=incl_autocor)
            gainsol_t.append([g_xx, g_yy])
            noise_t.append([n_xx, n_yy])
        gainsolutions.append(gainsol_t)
        noisesolutions.append(noise_t)
    gainsolutions = numpy.asarray(gainsolutions)
    noisesolutions = numpy.asarray(noisesolutions)
    return gainsolutions, noisesolutions


def autocorr_gain_solve(vis_uncal, hdsm_file=None):
    """\
    Solve for gains using autocorrelation over sidereal day

    This function solves r = g*m + n, with r the measured autocorrelations,
    m the model sky powers, g the power gains and n the noise power.
    Note that to apply the calibration to measured power data, one needs to
    invert this formula, so r_cal = 1/g*r-n/g.

    Parameters
    ----------
    vis_uncal : VisDataset
        The uncalibrated visibilities to be used to find gains and noise.
    hsdm_file: str
        Hemispheric diffuse model file name.

    Returns
    -------
    powgains : array
        Power gain solutions.
    pownoise : array
        Noise power solutions.
    """
    if vis_uncal.attrs['datatype'] != 'sst':
        # Extract auto-correlations from XST
        raise NotImplementedError(vis_uncal.attrs['datatype'])
    acc = vis_uncal.get_values(squash=True)
    #acc = vis_uncal.values.reshape(len(vis_uncal), *vis_uncal.shape()[2:])
    nrant = acc.shape[-1]
    nrsb = acc.shape[1]
    if not hdsm_file:
        abs_positions = vis_uncal.attrs['positions']
        centroid, rel_pos  = visibilities.layout_abs2rel(abs_positions)
        geopos = ITRF2lonlat(*centroid)
    else:
        geopos = hdsm_file
    hdsm = skymodels.HemiDiffuseSkyModel(geopos)
    del_ts = vis_uncal.delta_time.flatten()
    ll, mm = imaging.lmgrid(hdsm.imsize)
    bmjones = beam.horizontaldipoles_jones(ll, mm, rotzen=np.deg2rad(-45.0))
    powgains = np.zeros((nrsb, 2, nrant))
    pownoise = np.zeros_like(powgains)
    lc_mod = np.zeros((len(del_ts), nrsb, 2))
    relerr = np.zeros((nrsb, 2, nrant))
    # For each s and pol, solve eq.:
    #     w_isp .* [lcmod_isp, 1_isp] * [G_spn, N_spn] = w_isp .* lcmes_ispn
    w = np.diag(vis_uncal.flag2weights(['delta_time']).squeeze().flatten())
    for sb in range(nrsb):
        print('On subband ', sb, '/', nrsb)
        freq = vis_uncal.frequencies[sb]
        lsts, lc_mod_sb_XX, lc_mod_sb_YY, dat00 = \
            skymodels.lightcurve_bl(freq, hdsm, del_ts,
                                    vis_uncal.attrs['start_datetime'], bmjones)
        lc_mod[:, sb, 0] = lc_mod_sb_XX
        lc_mod[:, sb, 1] = lc_mod_sb_YY
        for pol in range(2):
            # Solve lin eq a*x=b for XX and YY polarizations, where
            #   a=[lc_mod(UT), 1] and b=acc(UT,elem).
            # (Previously UT index was sorted by LST, but not needed)
            a = np.vstack([lc_mod[:, sb, pol], np.ones(lc_mod.shape[0])]).T
            b = acc[:, sb, pol, :].squeeze()
            lsq_res = np.linalg.lstsq(w @ a, w @ b)
            powgains[sb, pol, :] = lsq_res[0][0, :]
            pownoise[sb, pol, :] = lsq_res[0][1, :]
            # Calculate rel.err. of solutions as
            lc_cal = (b - pownoise[sb, pol, :]) * (1 / powgains[sb, pol, :])
            relerr[sb, pol, :] = np.linalg.norm(
                lc_mod[:, sb, pol][..., np.newaxis] - lc_cal, axis=0
            ) / np.linalg.norm(lc_cal)
    # Axes have form [sampnr, sbnr, polnr, antnr]
    gains = numpy.emath.sqrt(powgains)
    return gains, pownoise, lc_mod, relerr

import argparse


def solvegains_cli():
    """
    Compute gain solutions for CVC files

    Parameters
    ----------
    variant: str
        Variant of WALS to use. Can be 'legacy' or 'inv'.

    Returns
    -------
    gainsolutions : array
        Gain solutions
    noises: array
        Noise solutions
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--legacy_variant', action="store_true",
                        help="""If raised, use legacy cal variant, rather than\
                                inverse 'inv' variant.""")
    parser.add_argument('dataff',
                        help="""Path to CVC folder""")
    args = parser.parse_args()
    variant = 'inv'
    if args.legacy_variant:
        variant = 'legacy'
    dataff_raw, dataff_mod, _dataff_cal = \
        data_io.dataff_raw_model_cal(args.dataff)
    if not os.path.isdir(dataff_mod):
        raise IsADirectoryError('No model file-folder {} found.'.format(dataff_mod))
    # Direct output to '_mod' data file-folder
    cvcobj_raw = data_io.CVCfiles(dataff_raw)
    cvcobj_model = data_io.CVCfiles(dataff_mod)
    gainsolutions, noises = gainsolve(cvcobj_raw, cvcobj_model,
                                      wals_variant=variant)
    gsol_filename = 'gain_solutions'
    noise_filename = 'noises_raw'
    if variant == 'inv':
        gsol_filename = 'gain_inv_solutions'
        noise_filename = 'noises_model'
    gsol_file = os.path.join(dataff_mod, gsol_filename)
    noise_file = os.path.join(dataff_mod, noise_filename)
    numpy.save(gsol_file, gainsolutions)
    numpy.save(noise_file, noises)
    return gainsolutions, noises


def applygains_cli():
    """
    Apply a calibration file to ACC or XST data folder.

    Returns
    -------
    cvcobj_cal : CVCfiles
        Calibrated CVCfiles.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-l', '--legacy_variant', action="store_true",
                        help="""If raised, use legacy cal variant, rather than\
                                inverse 'inv' variant.""")
    parser.add_argument('cvcpath', help="Path to CVC folder")
    parser.add_argument('caltabpath', nargs='?', default='',
                        help="Path to caltab file")

    args = parser.parse_args()
    variant = 'inv'
    if args.legacy_variant:
        variant = 'legacy'
    if args.caltabpath:
        applycal_cvcfolder(args.cvcpath, args.caltabpath)
    else:
        cvcobj_cal = apply_polgains_cvcfolder(args.cvcpath, variant=variant)
    return cvcobj_cal


def main_cli():
    cmd = sys.argv.pop(1)
    ran_sol = False
    ran_app = False
    if cmd.startswith('sol'):
        print("Solving for gains...")
        try:
            solvegains_cli()
        except IsADirectoryError as err:
            print(err)
            print('Please generate model first.')
            sys.exit(1)
        ran_sol = True
    if cmd.endswith('app'):
        print("Applying solutions...")
        applygains_cli()
        ran_app = True
    if not ran_sol and not ran_app:
        print("Nothing chosen.\n"
              "Choose: 'sol' to solve, 'app' to apply solutions.")


if __name__ == "__main__":
    main_cli()
