import shutil

import numpy

import ilisa.calim
from ilisa.antennameta import calibrationtables as calibrationtables
from ilisa.observations import dataIO as dataIO, modeparms as modeparms


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
        cvcdata_unc = cvcobj_cal.covcube_fb(filestep, crlpolrep=None)
        # Apply calibration
        cvcdata = ilisa.calim.calibration.applycaltab_cvc(cvcdata_unc, caltab,
                                                          sb)
        # Replace uncalibrated data with calibrated:
        cvcobj_cal.dataset[filestep] = cvcdata
    cvcobj_cal.save_ldat()
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
