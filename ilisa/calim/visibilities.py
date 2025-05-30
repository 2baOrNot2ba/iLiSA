import datetime
import sys

import casacore.measures
import numpy as np

from . import SPEED_OF_LIGHT as c
from ..operations.directions import directionterm2tuple

try:
    import dreambeam
    canuse_dreambeam = True
except ImportError:
    canuse_dreambeam = False
if canuse_dreambeam:
    from dreambeam.polarimetry import cov_lin2cir, convertxy2stokes


class VisDatasetFile:
    """"""
    coord_names = ('delta_time', 'frequencies')
    _axes_all = ('filenr',         # 0
                 'dt_perfile',     # 1
                 coord_names[1],   # 2
                 'polarization1',  # 3
                 'polarization2',  # 4
                 'antenna1',       # 5
                 'antenna2'        # 6
                 )
    _multi_axes = {coord_names[0]: (0,1),
                   'rcu1': (3,5),
                   'rcu2': (4,5)}
    _axes_gen = {'full': _axes_all,
                 'autocorr': (*_axes_all[0:4], _axes_all[5])}

    def __init__(self, *data_files, coords=None, attrs=None, vis_type='full'):
        self.values = data_files
        self.delta_time = coords[VisDatasetFile.coord_names[0]]
        self.frequencies = coords[VisDatasetFile.coord_names[1]]
        self.attrs = attrs
        dt_shape = list(self.delta_time.shape)
        if len(dt_shape) > 1:
            self.nrfile = dt_shape.pop(0)
        self.nrsmpsfile = dt_shape.pop(0)
        self.vis_type = vis_type
        self._dim_basis = VisDatasetFile._axes_gen[self.vis_type]
        _dt_axes = VisDatasetFile._multi_axes[VisDatasetFile.coord_names[0]]
        # Create virtual axes with delta_time axes replacing filenr, dt_perfile
        __dim_basis = list(self._dim_basis)
        for _axis in reversed(_dt_axes):
            __dim_basis.pop(_axis)
        __dim_basis.insert(sorted(_dt_axes)[0], VisDatasetFile.coord_names[0])
        self.dims = tuple(__dim_basis)

    @classmethod
    def load_from_np(cls, npfile):
        ds_np = np.load(npfile, allow_pickle=True)
        arr_names = tuple(filter(lambda s: s.startswith('arr_'), ds_np.files))
        data_arrs = (ds_np[en] for en in arr_names)
        coords = {nm: ds_np[nm] for nm in VisDatasetFile.coord_names}
        attrs = {k: ds_np[k] for k in ds_np.files
                 if k not in arr_names+VisDatasetFile.coord_names}
        vis_type = 'full'
        if attrs['telescope'] == 'LOFAR' and attrs['datatype'] == 'sst':
            vis_type = 'autocorr'
        return cls(*data_arrs, coords=coords, attrs=attrs, vis_type=vis_type)

    def save_to_np(self, filename):
            np.savez_compressed(filename, *self.values, **self._get_coords(),
                                **self.attrs)

    def __len__(self):
        """Returns number of temporal samples"""
        return self.nrfile * self.nrsmpsfile

    def __deepcopy__(self):
        new_vis = VisDatasetFile(*self.values, coords=self._get_coords(),
                             attrs=self.attrs, vis_type=self.vis_type)
        return new_vis

    def _get_coords(self):
        coords = {self.coord_names[0]: self.delta_time,
                  self.coord_names[1]: self.frequencies}
        return coords

    def shape(self, squashed=False):
        shape = (self.nrfile,) + self.values[0].shape
        if squashed:
            shape = (self.__len__(), *self.values[0].shape[1:])
        return shape

    def get_values(self, squash=True, wflagdims=None):
        values = np.asarray(self.values).reshape(self.shape(squash))
        if wflagdims is not None:
            wf = self.flag2weights(wflagdims, squash)
            values *= wf
        return values

    def refile(self, nrfile=1):
        new_vd = self.__deepcopy__()
        new_vd.values = tuple(new_vd.values.reshape(nrfile, -1, *self.shape()[2:]))
        new_vd.nrfile = nrfile
        new_vd.nrsmpsfile = self.__len__() // nrfile
        # Coords
        new_vd.delta_time = new_vd.delta_time.reshape(new_vd.nrfile, new_vd.nrsmpsfile)
        if new_vd.frequencies.shape > 1:
            new_vd.frequencies = new_vd.frequencies.reshape(new_vd.nrfile, new_vd.nrsmpsfile)
        return new_vd

    def _flagondim(self, dim=None, indices=None):
        """Flag indices on a dimension

        Parameters
        ----------
        dim : str
            Name of the dimension to add flags too.
        indices : list of ints
            List of indices along dimension `dim` to flag.
        """
        if self.attrs.get('flagaxes') is None:
            self.attrs['flagaxes'] = {}
        if dim is not None:
            if indices is not None:
                self.attrs['flagaxes'][dim] = indices
            else:
                del self.attrs['flagaxes'][dim]
        else:
            del self.attrs['flagaxes']

    def flagfromfile(self, filename):
        """
        Use flags specified in text file

        Parameters
        ----------
        filename: str
            Name of the flag file.
        """
        flagaxes = {}
        with open(filename, 'r') as f:
            for l in f.readlines():
                if l.startswith('#'): continue
                axis, indices_str = l.rstrip().split(' ', 1)
                indices = eval('np.r_'+indices_str)
                flagaxes[axis] = indices
        for _axis in flagaxes:
            self._flagondim(_axis, flagaxes[_axis])

    def flag2weights(self, dims=None, squashfile=True):
        """Get Weights as per Flags and Axes"""
        flagaxes = self.attrs.get('flagaxes', np.array({})).item()
        wvec_axes = [np.ones(1)] * len(self._dim_basis)
        if dims == '*':
            dims = self.dims
        for dim in dims:
            flagindices = flagaxes.get(dim, [])
            if dim in VisDatasetFile._multi_axes:
                merge_axis_0, merge_axis_1 = VisDatasetFile._multi_axes[dim]
                # Map across file axis
                wvv = np.tensordot(np.ones(self.shape(False)[merge_axis_0]),
                                   np.ones(self.shape(False)[merge_axis_1]),
                                   axes=0)
                if flagindices != []:
                    idxs = np.unravel_index(flagindices, wvv.shape)
                    wvv[idxs] = 0.0
                wvec_axes[0] = wvv
                wvec_axes[1] = None
            else:
                _axis = self._dim_basis.index(dim)
                nrcmp = self.shape(False)[_axis]
                wv = np.ones(nrcmp)
                if flagindices != []:
                    wv[flagindices] = 0.0
                wvec_axes[_axis] = wv
        wvec_axes = list(filter(lambda x: x is not None, wvec_axes))
        # Compute einsum index expression
        p = ord('i')
        idxxpr = []
        for b in (np.ndim(wv) for wv in wvec_axes):
            l = ''
            for bi in range(b):
                c = chr(p)
                l += c
                p += 1
            idxxpr.append(l)
        idxxpr = ','.join(idxxpr)
        w = np.einsum(idxxpr, *wvec_axes)
        if squashfile:
            w = w.reshape(w.shape[0]*w.shape[1], *w.shape[2:])
        return w


def fiducial_visibility(nrelems=2):
    off_diag_val = 0.5
    vis = np.ones((nrelems, nrelems), dtype=complex)
    vis = off_diag_val * vis
    np.fill_diagonal(vis, 1.0)
    return vis


def point_source_vis2d(uv_lam, l0=0.0, m0=0.0, amp=1.0):
    # For point src at l0,m0
    # Vis(U,V)=exp(i*(U*l0+V*m0))
    arrayvec = np.exp(-1j*2*np.pi*(uv_lam[:, 0]*l0+uv_lam[:, 1]*m0))
    n0 = np.sqrt(1-l0**2-m0**2)
    vis = amp*np.einsum('i,k->ik', arrayvec, np.conj(arrayvec))/n0
    return vis


def baseline_distances(uvw_xyz):
    dist = np.linalg.norm(uvw_xyz[None, :, :] - uvw_xyz[:, None, :], axis=-1)
    return dist


def select_baselines_by_dist(uvw_xyz, freq, crit):
    dist_uvw = baseline_distances(uvw_xyz) /(c / freq)
    sel_mat = eval('dist_uvw'+crit)
    return sel_mat


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
    >>> v0=np.arange(1*2*2*2*2).reshape((1,2,2,2,2))
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

    while if on the other hand the lower half is pol 0 and upper half pol 1:
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
        cvcflat = np.empty(flatshape, dtype=cvcpol.dtype)
        cvcflat[..., ::2, ::2] = pp
        cvcflat[..., 1::2, 1::2] = qq
        cvcflat[..., ::2, 1::2] = pq
        cvcflat[..., 1::2, ::2] = qp
    else:
        cvcflat = np.block([[pp, pq], [qp, qq]])
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
    >>> v0=np.arange(1*2*2*2*2).reshape((1,2*2,2*2))
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
    while if on the other hand the lower half is pol 0 and upper half pol 1:
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
    cvpol = np.array([[pp, pq], [qp, qq]])
    # Move pol-axes from leftmost axises rightward just before element indexes:
    cvpol = np.moveaxis(cvpol, 1, -3)
    cvpol = np.moveaxis(cvpol, 0, -4)
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
        cvpol = np.moveaxis(cvpol, 1, -3)
        cvpol = np.moveaxis(cvpol, 0, -4)
    elif crlpolrep == 'sto':
        (S0, S1, S2, S3) =\
            convertxy2stokes(cvcpolidx[0][0], cvcpolidx[0][1], cvcpolidx[1][0],
                             cvcpolidx[1][1])
        cvpol = np.array([S0, S1, S2, S3])
        cvpol = np.moveaxis(cvpol, 0, -3)
    else:
        raise ValueError("No such correlated pol rep: '{}'".format(crlpolrep))
    return cvpol


def pos2uv_flat(antpos):
    """
    Convert 3D antenna positions (approx flat config) to 2D, UV baseline

    Parameters
    ----------
    antpos: (nrants, 3) array
        3D positions of antennas

    Returns
    -------
    uv : (nrants, nrants, 3) array
        Flat 3D baselines of configuration
    """
    # antpos is of the form [nants, 3]
    nants = antpos.shape[0]
    # Compute baselines in XYZ
    ant_pos_rep = np.repeat(antpos[:, :], nants, axis=0
                          ).reshape((nants, nants, 3))
    # Compute square baseline tensor
    _blsq = (ant_pos_rep - np.transpose(ant_pos_rep, (1, 0, 2)))
    # Vectorize square baseline matrix
    _idx = np.triu_indices(ant_pos_rep.shape[0])
    baselines = _blsq[_idx]

    # Compute local UV coord sys by first estimating normal to
    # xyz (uvw) array, using the fact that its normal vec is a
    # null space vector.
    _u_svd, _d_svd, _vt_svd = np.linalg.svd(baselines)
    nrmvec = -_vt_svd[2, :] / np.linalg.norm(_vt_svd[2, :])
    lon_nrm = np.arctan2(nrmvec[1], nrmvec[0])
    lat_nrm = np.arcsin(nrmvec[2])
    # Transform by rotations xyz to local UV crd sys, which has
    # normal along its z-axis and long-axis along x-axis:
    # First rotate around z so nrmvec x is in (+x,+z) quadrant
    # (this means Easting is normal to longitude 0 meridian plane)
    _rz = np.array([[np.cos(lon_nrm), np.sin(lon_nrm), 0.],
                    [-np.sin(lon_nrm), np.cos(lon_nrm), 0.],
                    [0., 0., 1.]])
    # Second rotate around y until normal vec is along z (NCP)
    _tht = (np.pi / 2 - lat_nrm)
    _ry = np.array([[np.cos(_tht), 0., -np.sin(_tht)],
                    [0., 1., 0.],
                    [+np.sin(_tht), 0., np.cos(_tht)]])
    # Third rotate around z so Easting is along (final) x
    _r_xy90 = np.array([[0, 1, 0], [-1, 0, 0], [0, 0, 1]])
    rot_mat = _r_xy90 @ _ry @ _rz
    uv = baselines @ rot_mat.T
    return uv


def calc_uvw(obstime, phaseref, stn_pos, stn_antpos):
    """
    Calculate UVW coords from datetime, phaseref and station & antenna positions

    Parameters
    ----------
    obstime : datetime
        Datetime of observation. Type can be datetime.datetime or
        np.datetime64.
    phaseref : tuple
        Phase reference direction given by (lon, lat, ref),
        where lon, lat are floats for longitude, latitude in radians and ref is
        the reference frame str.
    stn_pos : array_like
        Vector, with components X,Y,Z in axis 0, giving station's ITRF position
        in meters.
    stn_antpos : array_like
        Matrix over antenna versus X,Y,Z  of ITRF positions relative stn_pos.

    Returns
    -------
    uvw_xyz : array_like
        UVW coordinates in meters.
    """
    (pntRA, pntDEC, pntref) = phaseref
    stn_pos = np.asarray(stn_pos)
    pos_ITRF_X = str(stn_pos[0].squeeze())+'m'
    pos_ITRF_Y = str(stn_pos[1].squeeze())+'m'
    pos_ITRF_Z = str(stn_pos[2].squeeze())+'m'
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
                              np.asarray(stn_antpos[antnr, :]).squeeze()])
        bls.append(bl)
    uvw_xyz = np.zeros((nrant,3))
    # Set obstime
    if type(obstime) is np.datetime64:
        obstime = datetime.datetime.fromtimestamp(
            (obstime - np.datetime64('1970-01-01T00:00:00Z'))
                  / np.timedelta64(1, 's'), datetime.timezone.utc)
    when = obsme.epoch("UTC", obstime.isoformat('T'))
    obsme.doframe(when)
    for antnr in range(nrant):
        uvw_xyz[antnr,:] = np.asarray(
                        obsme.to_uvw(bls[antnr])["xyz"].get_value('m'))
    return uvw_xyz


def layout_abs2rel(layout_abs):
    """\
    Convert absolute layout positions to positions relative a centroid

    Parameters
    ----------
    layout_abs : array
        Array of absolute 3D positions of layout

    Returns
    -------
    centroid : array
            Absolute vector to layout centroid
    layout_rel : array
        Array of relative 3D positions of layout
    """
    # Compute barycenter of abs 3D positions using mean:
    centroid = np.mean(layout_abs, axis=0)
    layout_rel = layout_abs - centroid
    return centroid, layout_rel


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


def phaseref_xstpol(xstpol, UVWxyz, freq):
    """
    Phase up polarized visibilities stack to U,V-align them at frequency
    """
    lambda0 = sys.float_info.max
    if freq != 0.0:
        lambda0 = c / freq
    # Phase-factor corresponds to exp(-i*2*pi*W/lambda),
    # this is so that B(l,m)*exp(i*2*pi*(U*l+V*m+W*n)/lambda)*phasefactor=
    #   B(l,m)*exp(i*2*pi*(U*l+V*m+W*(n-1))/lambda) which
    phasefactors = np.exp(-2.0j*np.pi*UVWxyz[:,2]/lambda0)
    PP = np.einsum('i,k->ik', phasefactors, np.conj(phasefactors))
    #xstpupol = np.array(
    #       [[PP*xstpol[0, 0, ...].squeeze(), PP*xstpol[0, 1, ...].squeeze()],
    #        [PP*xstpol[1, 0, ...].squeeze(), PP*xstpol[1, 1, ...].squeeze()]])
    xstpupol = PP*xstpol
    return xstpupol


def phasedup_vis(vis, srcname, t, freq, stn_pos, stn_antpos):
    """
    Phase up visibiliies
    """
    # Phase center on src
    dir_src = directionterm2tuple(srcname)
    uvw_src = calc_uvw(t, dir_src, stn_pos, stn_antpos)
    vis_pu = phaseref_xstpol(vis, uvw_src, freq)
    return vis_pu
