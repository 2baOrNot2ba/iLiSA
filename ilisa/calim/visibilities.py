import numpy
try:
    import dreambeam
    canuse_dreambeam = True
except ImportError:
    canuse_dreambeam = False
if canuse_dreambeam:
    from dreambeam.polarimetry import cov_lin2cir, convertxy2stokes


def cov_polidx2flat(cvcpol, parity_ord=True):
    """
    Convert polarization indexed array covariance matrix to flat array

    Parameters
    ----------
    cvcpol : array_like
        Polarization index array covariance cube to be converted.
        Shape will be (2, 2, ..., N, N).
    parity_ord : Boolean
        If True (default) then the polarization component is determined from
        the index's parity: even index maps to component 0, odd to 1.
        If False the baseline indices are split into a first and second half
        and mapped to pol component 0,1 respectively.

    Returns
    -------
    cvc : array_like
        Converted flat covariance cube. Shape should be (..., 2*N, 2*N),
        where N is the number of dual-polarized elements.

    Examples
    --------
    Minimalist example:
    >>> cov_polidx2flat(numpy.arange(4*4).reshape((2,2,2,2)))
    array([[ 0.,  4.,  1.,  5.],
           [ 8., 12.,  9., 13.],
           [ 2.,  6.,  3.,  7.],
           [10., 14., 11., 15.]])

    See Also
    --------
    cov_flat2polidx : Inverse function.
    """
    pp = cvcpol[0, 0, ...]
    qq = cvcpol[1, 1, ...]
    pq = cvcpol[0, 1, ...]
    qp = cvcpol[1, 0, ...]
    if parity_ord:
        N = cvcpol.shape[-1]
        restshape = cvcpol.shape[2:-2]
        flatshape = restshape + (2 * N, 2 * N)
        cvcflat = numpy.zeros(flatshape, dtype=numpy.complex)
        cvcflat[..., ::2, ::2] = pp
        cvcflat[..., 1::2, 1::2] = qq
        cvcflat[..., ::2, 1::2] = pq
        cvcflat[..., 1::2, ::2] = qp
    else:
        cvcflat = numpy.block([[pp, pq], [qp, qq]])
    return cvcflat


def cov_flat2polidx(cvc, parity_ord=True):
    """
    Convert flat array covariance matrix (visibilities) to polarization
    component indexed array covariance matrix.

    Parameters
    ----------
    cvc : array_like
        Covariance cube to be converted. Shape should be (..., 2*N, 2*N),
        where N is the number of dual-polarized elements.
    parity_ord : Boolean
        If True (default) then the polarization component is determined from
        the index's parity: even index maps to component 0, odd to 1.
        If False the baseline indices are split into a first and second half
        and mapped to pol component 0,1 respectively.

    Returns
    -------
    cvcpol : array_like
        Polarization index array covariance cube. Shape will be
        (2, 2, ..., N, N).

    Notes
    -----
    This function is agnostic to whether the polarization basis is
    linear or circular. Also the ordering is conserved from the
    flat ordering, so the mapping of index 0,1 to say L,R or X,Y
    components is determined from the original ordering.

    Examples
    --------
    Minimal example:
    >>> cov_flat2polidx(numpy.arange(4).reshape((2,2)))[1,0,...]
    array([[2]])

    See Also
    --------
    cov_polidx2flat : Inverse function.
    """
    if parity_ord:
        pp = cvc[..., ::2, ::2]
        qq = cvc[..., 1::2, 1::2]
        pq = cvc[..., ::2, 1::2]
        qp = cvc[..., 1::2, ::2]
    else:
        # First-half, second-half order
        n = cvc.shape[-1]//2
        pp = cvc[..., :n, :n]
        qq = cvc[..., n:, n:]
        pq = cvc[..., :n, n:]
        qp = cvc[..., n:, :n]
    cvpol = numpy.array([[pp, pq], [qp, qq]])
    return cvpol


def cvc2polrep(cvc, crlpolrep='lin'):
    """
    Convert a flat-indexed polarized covariance cube `cvc` into an array
    indexed according to correlated polarization representation
    `crlpolrep`.

    Parameters
    ----------
    cvc : (T,2*N,2*N) array (usually T=512 and 2*N=196)
        The Covariance Cube array produced by an International LOFAR
        station when it is in calibration mode. It is the covariance
        matrices of the 196 rcus (98 X-polarized & 98 Y-polarized
        interleaved) over 512subbands.
    crlpolrep : str
        Correlated polarization representation can be 'lin', 'cir', 'sto' where:
            'lin' is linear pol rep
            'cir' is circular pol rep
            'sto' is Stokes pol rep.

    Returns
    -------
    cvpol : array
        Polarized visibilities. Shape depends on `crlpolrep`:
            * (2,2,T,N,N) for 'lin' or 'cir'. Basis is indexed with
              0,1 corresponding to X,Y or L,R resp.
            * (4,T,N,N) for 'sto'. Basis is indexed with
              0,1,2,3 corresponding to Stokes I,Q,U,V resp.

    Notes
    -----
    Requires dreamBeam.
    """
    cvcpolidx = cov_flat2polidx(cvc)
    if crlpolrep == 'lin':
        cvpol = cvcpolidx
    elif crlpolrep == 'cir':
        cvpol = cov_lin2cir(cvcpolidx)
    elif crlpolrep == 'sto':
        (S0, S1, S2, S3) =\
            convertxy2stokes(cvcpolidx[0][0], cvcpolidx[0][1], cvcpolidx[1][0],
                             cvcpolidx[1][1])
        cvpol = numpy.array([S0, S1, S2, S3])
    else:
        raise ValueError("No such correlated pol rep: '{}'".format(crlpolrep))
    return cvpol
