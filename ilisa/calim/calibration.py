import sys
import os
import shutil

import numpy
import matplotlib.pyplot as plt

from ilisa.antennameta import calibrationtables as calibrationtables
from ilisa.operations import data_io as data_io, modeparms as modeparms
from . import visibilities


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


def apply_polgains(vispol, gainspol):
    """
    Apply polarized biscalar gains to polarized visibility

    Parameters
    ----------
    vispol : array_like
        Polarized visibility array. Its shape is (2,2,N,N)
    gainspol : array_like
        Polarized biscalar gains. g_pq has shape (2,N) where N is number of
        dual-polarized elements.

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
    #    (T,2,N),conj(T,2,N) => (T,2,1,N,1)*conj(T,1,2,1,N) => (T,2,2,N,N)
    gg = (gainspol[...,    :, None,    :, None] * numpy.conj(
          gainspol[..., None,    :, None,    :]))
    vispol_gapp = gg * vispol
    return vispol_gapp


def apply_polgains_cvcfolder(dataff, gainsolfile='gainsolutions.npy'):
    """
    Apply polarimetric gains to CVC data

    Parameters
    ----------
    dataff_mod: str
        Name of file-folder
    gainsolfile: str
        Name of gain solutions in model vis directory

    Returns
    -------
    cvcobj_cal: CVCfiles
        The calibrated CVCfiles object
    """
    ldat_type = data_io.datafolder_type(dataff)
    if ldat_type != "acc" and ldat_type != "xst":
        raise ValueError("Not CVC data.")
    dataff = os.path.normpath(dataff)
    dataff_dir, dataff = os.path.split(dataff)
    obsinfo_raw = data_io.filefolder2obsinfo(dataff)
    # Determine 'raw' and 'mod' dataffs,
    obsinfo_raw.pop('cal', None)
    obsinfo_raw.pop('model', None)
    dataff_raw = data_io.obsinfo2filefolder(obsinfo_raw)
    dataff_raw = os.path.join(dataff_dir, dataff_raw)
    # Model
    obsinfo_mod = dict(obsinfo_raw)
    obsinfo_mod['model'] = True
    dataff_mod = data_io.obsinfo2filefolder(obsinfo_mod)
    dataff_mod = os.path.join(dataff_dir, dataff_mod)
    # Cal
    obsinfo_cal = dict(obsinfo_raw)
    obsinfo_cal['cal'] = gainsolfile
    dataff_cal = data_io.obsinfo2filefolder(obsinfo_cal)
    dataff_cal = os.path.join(dataff_dir, dataff_cal)
    # Copy raw into cal
    shutil.copytree(dataff_raw, dataff_cal)
    # Read in cvcobj:
    cvcobj_cal = data_io.CVCfiles(dataff_cal)
    # Read in gain solutions
    gainsolpath = os.path.join(dataff_mod, gainsolfile)
    gainsols = numpy.load(gainsolpath)
    nrfiles = cvcobj_cal.getnrfiles()
    # Loop over files in CVC folder:
    for filestep in range(nrfiles):
        # Get actual covariance cubes:
        cvpol_unc = visibilities.cov_flat2polidx(cvcobj_cal[filestep])
        # Apply calibration
        _g = 1/numpy.conjugate(gainsols[filestep])
        cvcdata_cal = apply_polgains(cvpol_unc, _g)
        # Replace uncalibrated data with calibrated:
        cvcobj_cal[filestep] = visibilities.cov_polidx2flat(cvcdata_cal)
    # Note in ScanRecInfo about calibrating this dataset:
    cvcobj_cal.scanrecinfo.set_postcalibration(gainsolpath, dataff_cal)
    return cvcobj_cal


def stefcal(r, m, niter=100):
    """
    Compute solution to the quadratic matrix equation using the StefCal
    algorithm [stefs]_.

    Returns a solution, up to a phase factor, to the equation
    .. math:: g m g^H = r
    where g is a column vector and r, m are Hermitian matrices.
    To apply gain solutions g (as returned from this function) in order
    calibrate a measurement r, one should therefore use the formula
    .. math:: r_cal = (1/g) r (1/g)^H

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


def gainsolve(dataff, gs_model='LFSM'):
    """
    Solve for gains based on uncalibrated and model data

    Parameters
    ----------
    dataff: str
        Path to CVC filefolder.
    gs_model:  str
        Name Global skymodel. Choose 'LFSM', 'GSM'.
    """
    dataff = os.path.normpath(dataff)
    dataff_dir, dataff = os.path.split(dataff)
    lofar_datatype = data_io.datafolder_type(dataff)
    if lofar_datatype != 'acc' and lofar_datatype != 'xst':
        raise RuntimeError("Datafolder '{}'\n not ACC or XST type data."
                           .format(dataff))
    obsinfo_raw = data_io.filefolder2obsinfo(dataff)
    # Determine 'raw' and 'mod' dataffs,
    obsinfo_raw.pop('cal', None)
    obsinfo_raw.pop('model', None)
    dataff_raw = data_io.obsinfo2filefolder(obsinfo_raw)
    dataff_raw = os.path.join(dataff_dir, dataff_raw)
    obsinfo_mod = dict(obsinfo_raw)
    obsinfo_mod['model'] = gs_model
    dataff_mod = data_io.obsinfo2filefolder(obsinfo_mod)
    dataff_mod = os.path.join(dataff_dir, dataff_mod)
    # Direct output to '_mod' data file-folder
    cvcobj_uncal = data_io.CVCfiles(dataff_raw)
    cvcobj_model = data_io.CVCfiles(dataff_mod)
    gainsolutions = []
    nrfiles = cvcobj_uncal.getnrfiles()
    for fileidx in range(nrfiles):
        print('Solving for file {}/{}'.format(fileidx, nrfiles))
        intgs = len(cvcobj_uncal.samptimeset[fileidx])
        gainsol_t = []
        cvcpol_uncal = visibilities.cov_flat2polidx(cvcobj_uncal[fileidx])
        cvcpol_model = visibilities.cov_flat2polidx(cvcobj_model[fileidx])
        for tidx in range(intgs):
            t = cvcobj_uncal.samptimeset[fileidx][tidx]
            freq = cvcobj_uncal.freqset[fileidx][tidx]
            vis_uncal = cvcpol_uncal[tidx]
            vis_model = cvcpol_model[tidx]
            g_xx = stefcal(vis_uncal[0, 0, ...], vis_model[0, 0, ...])
            g_yy = stefcal(vis_uncal[1, 1, ...], vis_model[1, 1, ...])
            gainsol_t.append([g_xx, g_yy])
        gainsolutions.append(gainsol_t)
    gainsolutions = numpy.asarray(gainsolutions)
    gsol_file = os.path.join(dataff_mod, 'gainsolutions')
    numpy.save(gsol_file, gainsolutions)


import argparse


def solvegains_cli():
    """
    Compute gain solutions for CVC files

    Returns
    -------
    gains : array
        Gain solutions
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('dataff',
                        help="""Path to CVC folder""")
    args = parser.parse_args()
    gains = gainsolve(args.dataff)
    return gains


def applygains_cli():
    """
    Apply a calibration file to ACC or XST data folder.

    Returns
    -------
    cvcobj_cal : CVCfiles
        Calibrated CVCfiles.
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('cvcpath', help="Path to CVC folder")
    parser.add_argument('caltabpath', nargs='?', default='',
                        help="Path to caltab file")
    args = parser.parse_args()
    if args.caltabpath:
        applycal_cvcfolder(args.cvcpath, args.caltabpath)
    else:
        cvcobj_cal = apply_polgains_cvcfolder(args.cvcpath)
    return cvcobj_cal


if __name__ == "__main__":
    cmd = sys.argv.pop(1)
    ran_sol = False
    ran_app = False
    if cmd.startswith('sol'):
        print("Solving for gains...")
        solvegains_cli()
        ran_sol = True
    if cmd.endswith('app'):
        print("Applying solutions...")
        applygains_cli()
        ran_app = True
    if not ran_sol and not ran_app:
        print("Nothing was done.")
