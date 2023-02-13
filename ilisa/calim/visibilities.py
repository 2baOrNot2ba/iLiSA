import casacore.measures
import numpy
import numpy as np
import numpy.ma as ma

try:
    import dreambeam
    canuse_dreambeam = True
except ImportError:
    canuse_dreambeam = False
if canuse_dreambeam:
    from dreambeam.polarimetry import cov_lin2cir, convertxy2stokes


def apply_vispol_flags(vispol, flags):
    """
    Flag visibilities

    Parameters
    ----------
    vispol: array
        Polarized visibility matrix.
    flag: dict
        Dict of flags.

    Returns
    -------
    vispol_flagged: array
        Flagged polarized visibilities.
    """
    flagbls = flags.get('bls')
    flagpols = flags.get('pol')
    if flagbls is None and flagpols is None:
        vispol_flagged = ma.array(vispol)
        return vispol_flagged
    if flagpols is None:
        polshape = vispol.shape[:-2]
        flagpols = numpy.zeros(polshape, dtype=bool)
    else:
        # flagbls is None, so set to unselected
        blshape = vispol.shape[-2:]
        flagbls = numpy.zeros(blshape, dtype=bool)
    flagpolmat = numpy.expand_dims(flagpols, axis=(-1, -2))
    flagblspolmat = flagbls + flagpolmat
    vispol_flagged = ma.array(vispol, mask=flagblspolmat)
    return vispol_flagged


def select_cov_mask(selections, nr_cov_el):
    """
    From a covariance selection specification, make a mask matrix

    This function can be used for covariances such as baselines or
    polarizations.

    Parameters
    ----------
    selections: list
        Covariance selection specification.
        Integeters in list select corresponding elements.
        Tuples of length 1 select auto-correlation elements.
        Tuples of length 2, where 2nd item is None, select all covariance
        between all elements in 1st item slot.
        Tuples of length 2, select all covariances with one element from 1st
        list and the other element from 2nd list.

    nr_cov_el: int
        Number of covariance elements.

    Returns
    -------
    maskmat: array
        Matrix to use as mask with correlation matrix to effect selections.

    Examples
    --------
    >>> import ilisa.calim.visibilities as vsb
    Select antenna 2 and 3:
    >>> vsb.select_cov_mask([2,3], 4)
    array([[False, False,  True,  True],
       [False, False,  True,  True],
       [ True,  True,  True,  True],
       [ True,  True,  True,  True]])
    Select auto-correlation of antenna 1:
    >>> vsb.select_cov_mask([(1,)], 4)
    array([[False, False, False, False],
       [False,  True, False, False],
       [False, False, False, False],
       [False, False, False, False]])
    Select baseline 0-3:
    >>> vsb.select_cov_mask([(0,3)], 4)
    array([[False, False, False,  True],
       [False, False, False, False],
       [False, False, False, False],
       [ True, False, False, False]])
    Select all auto-correlations:
    >>> vsb.select_cov_mask([(None,)], 4)
    array([[ True, False, False, False],
           [False,  True, False, False],
           [False, False,  True, False],
           [False, False, False,  True]])

    """
    maskmat = np.zeros((nr_cov_el, nr_cov_el), dtype=bool)
    for sel in selections:
        if type(sel) is int:
            # Antenna select
            maskmat[sel, :] = True
            maskmat[:, sel] = True
        elif type(sel) is tuple:
            if len(sel) == 1:
                if sel[0] is None:
                    sel =  (range(nr_cov_el),)
                # Autocorrelations
                maskmat[sel[0], sel[0]] = True
            elif len(sel) == 2:
                if sel[1] is None:
                    if sel[0] is not None:
                        _sel = tuple(numpy.meshgrid(sel[0], sel[0]))
                        maskmat[_sel] = True
                    else:
                        maskmat[:, :] = True
                else:
                    maskmat[sel[0], sel[1]] = True
                    maskmat[sel[1], sel[0]] = True
    return maskmat


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
    For one sample and 2 elements, an example where the two polarization
    indices are embedded along with elements as in
    >>> v0=numpy.arange(1*2*2*2*2).reshape((1,2,2,2,2))
    array([[[[[ 0,  1],
          [ 2,  3]],
         [[ 4,  5],
          [ 6,  7]]],
        [[[ 8,  9],
          [10, 11]],
         [[12, 13],
          [14, 15]]]]])

    then if odd pols are index 0 and even pols are index 1, then
    >>> cov_polidx2flat(v0)
    array([[[ 0,  4,  1,  5],
        [ 8, 12,  9, 13],
        [ 2,  6,  3,  7],
        [10, 14, 11, 15]]])

    while if one the other hand the lower half is pol 0 and upper half pol 1:
    >>> cov_polidx2flat(v0, parity_ord=False)
    array([[[ 0,  1,  4,  5],
            [ 2,  3,  6,  7],
            [ 8,  9, 12, 13],
            [10, 11, 14, 15]]])

    See Also
    --------
    cov_flat2polidx : Inverse function.
    """
    pp = cvcpol[..., 0, 0, :, :]
    qq = cvcpol[..., 1, 1, :, :]
    pq = cvcpol[..., 0, 1, :, :]
    qp = cvcpol[..., 1, 0, :, :]
    if parity_ord:
        N = cvcpol.shape[-1]
        restshape = cvcpol.shape[:-4]
        flatshape = restshape + (2 * N, 2 * N)
        cvcflat = numpy.empty(flatshape, dtype=cvcpol.dtype)
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
        (..., 2, 2, N, N).

    Notes
    -----
    This function is agnostic to whether the polarization basis is
    linear or circular. Also the ordering is conserved from the
    flat ordering, so the mapping of index 0,1 to say L,R or X,Y
    components is determined from the original ordering.

    Examples
    --------
    For one sample and 2 elements, an example where the two polarization
    indices are embedded along with elements as in
    >>> v0=numpy.arange(1*2*2*2*2).reshape((1,2*2,2*2))
    array([[[ 0,  1,  2,  3],
        [ 4,  5,  6,  7],
        [ 8,  9, 10, 11],
        [12, 13, 14, 15]]])

    then if odd pols are index 0 and even pols are index 1, then
    >>> cov_flat2polidx(v0)
    array([[[[[ 0,  2],
          [ 8, 10]],
         [[ 1,  3],
          [ 9, 11]]],
        [[[ 4,  6],
          [12, 14]],
         [[ 5,  7],
          [13, 15]]]]])
    while if one the other hand the lower half is pol 0 and upper half pol 1:
    >>> cov_flat2polidx(v0, parity_ord=False)
    array([[[[[ 0,  1],
          [ 4,  5]],
         [[ 2,  3],
          [ 6,  7]]],
        [[[ 8,  9],
          [12, 13]],
         [[10, 11],
          [14, 15]]]]])

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
    # Move pol-axes from leftmost axises rightward just before element indexes:
    cvpol = numpy.moveaxis(cvpol, 1, -3)
    cvpol = numpy.moveaxis(cvpol, 0, -4)
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
        interleaved) over 512 subbands.
    crlpolrep : str
        Correlated polarization representation can be 'lin', 'cir', 'sto' where:
            'lin' is linear pol rep
            'cir' is circular pol rep
            'sto' is Stokes pol rep.

    Returns
    -------
    cvpol : array
        Polarized visibilities. Shape depends on `crlpolrep`:
            * (T,2,2,N,N) for 'lin' or 'cir'. Basis is indexed with
              0,1 corresponding to X,Y or L,R resp.
            * (T,4,N,N) for 'sto'. Basis is indexed with
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
        cvpol = numpy.moveaxis(cvpol, 1, -3)
        cvpol = numpy.moveaxis(cvpol, 0, -4)
    elif crlpolrep == 'sto':
        (S0, S1, S2, S3) =\
            convertxy2stokes(cvcpolidx[0][0], cvcpolidx[0][1], cvcpolidx[1][0],
                             cvcpolidx[1][1])
        cvpol = numpy.array([S0, S1, S2, S3])
        cvpol = numpy.moveaxis(cvpol, 0, -3)
    else:
        raise ValueError("No such correlated pol rep: '{}'".format(crlpolrep))
    return cvpol


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


def rot2uv(stn_antpos, stn_rot):
    """
    Calculate 2D UV coords from antenna positions and rotation

    Parameters
    ----------
    stn_antpos : array_like
        Matrix over antenna versus ITRF X,Y,Z postions relative stn_pos.
    stn_rot : array_like
        Vector of ITRF X,Y,Z coordinates in meters.

    Returns
    -------
    uv_xy : array_like
        UV coordinates in meters.
    """
    uv_xy = np.matmul(stn_antpos, stn_rot)[:,:3]
    return uv_xy


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
