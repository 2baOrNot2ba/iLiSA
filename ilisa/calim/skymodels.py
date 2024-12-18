import sys
import argparse
import os.path
import shutil
import warnings
import datetime

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import pylab
import healpy as hp
from scipy.constants import speed_of_light as c

from pygdsm.pygsm import GlobalSkyModel
from pygdsm.pygsm2016 import GlobalSkyModel2016
from pygdsm.lfsm import LowFrequencySkyModel
from pygdsm.base_observer import BaseObserver

from casacore.measures import measures

from ilisa.antennameta.export import ITRF2lonlat
import ilisa.operations.modeparms as modeparms
import ilisa.operations.data_io as data_io
from .imaging import plotskyimage, imggrid_res, fiducial_image
from .visibilities import cov_polidx2flat, calc_uvw, rot2uv,\
    point_source_vis2d, layout_abs2rel
from .beam import dualdipole45_cov_patt
from .im_process import n_from_lm


class CommonGSMObs(BaseObserver):
    def __init__(self, gsm_instance):
        """Initialize the common Observer object with a GSM instance
        Calls ephem.Observer.__init__ function and adds on gsm
        """
        super(BaseObserver, self).__init__()
        self.observed_sky = None
        self.gsm = gsm_instance
        # Inline replacement of self._setup() follows.
        # Generate mapping from pix <-> angles
        gen_freq_Hz = 100e6
        fu = gsm_instance.freq_unit
        gen_freq = 100
        if fu == 'Hz':
            gen_freq = 100e6
        elif fu == 'GHz':
            gen_freq = 0.1
        self.gsm.generate(gen_freq)
        self._n_pix  = hp.get_map_size(self.gsm.generated_map_data)
        self._n_side = hp.npix2nside(self._n_pix)
        self._theta, self._phi = hp.pix2ang(self._n_side, np.arange(self._n_pix))


def globaldiffuseskymodel(dattim, geopos, freq, gs_model='LFSM', imsize=200):
    """\
    Generate hemisphere sky model image based on GSM, epoch & geolocation

    Parameters
    ----------
    dattim: datetime
        Epoch for model sky
    geopos: tuple
        Longitude, latitude, elevation tuple of geographic position for model
    freq: float
        Frequency for which model should be generated
    gs_model: str
        ID of global sky model
    imsize: int
        Number of pixels along one dimension of image

    Returns
    -------
    img: array_like
        The generated sky model image
    """
    if freq < 10e6:
        warnings.warn('Freq =< 10 MHz, will use model for 10.1 MHz instead.')
        freq = 10.1e6
    (longitude, latitude, elevation) = geopos
    freq_unit = 'Hz'  # ('Hz', 'MHz', 'GHz')
    if gs_model =='LFSM':
        gsm = LowFrequencySkyModel(freq_unit=freq_unit)
    elif gs_model == 'GSM' or gs_model == 'GSM2008':
        gsm = GlobalSkyModel(freq_unit=freq_unit,
                             basemap='haslam',  # 'haslam', 'wmap' or '5deg'
                             interpolation='pchip'  # 'cubic' or 'pchip'
                             )
    else:
        gsm = GlobalSkyModel2016(freq_unit=freq_unit,
                                 data_unit='MJysr',  # ('TCMB', 'MJysr', 'TRJ')
                                 resolution='hi',  # ('hi', 'lo')
                                 theta_rot=0, phi_rot=0)
    # NOTE: CommonGSMObserver() is my addition to PyGDSM
    gsm_obs = CommonGSMObs(gsm)
    gsm_obs.lon = str(longitude)
    gsm_obs.lat = str(latitude)
    gsm_obs.elev = elevation
    gsm_obs.date = dattim
    try:
        gsm_map = gsm_obs.generate(freq)
    except RuntimeError as e:
        raise ValueError(e)
    f = pylab.figure(None, figsize=None)
    extent = (0.0, 0.0, 1.0, 1.0)
    ax = hp.projaxes.HpxOrthographicAxes(f, extent)
    img_ma = ax.projmap(gsm_map, xsize=imsize, half_sky=True)
    img = np.ma.getdata(img_ma)
    img[img == -np.inf] = 0.0
    img = np.fliplr(img)  # Sky-model has a flip along East-West, so flipback
    pylab.close()
    return img


def plot_gsm_for_obsdata(cvcobj, filenr=0, sampnr=0, gs_model='LFSM',
                         imsize=200):
    stnid = cvcobj.scanrecinfo.get_stnid()
    for fileidx in range(filenr, cvcobj.getnrfiles()):
        integration = cvcobj.scanrecinfo.get_integration()
        intgs = len(cvcobj.samptimeset[fileidx])
        for tidx in range(sampnr, intgs):
            t = cvcobj.samptimeset[fileidx][tidx]
            freq = cvcobj.freqset[fileidx][tidx]
            lon, lat, h = ITRF2lonlat(cvcobj.stn_pos[0, 0],
                                      cvcobj.stn_pos[1, 0],
                                      cvcobj.stn_pos[2, 0])
            plot_gsm(t, (lon, lat, h), freq, gs_model, imsize, stnid)


def plot_gsm(dattim, loc_lonlat, freq, gs_model='LFSM', imsize=200,
             stnid='Unknown'):
    lon, lat, h = loc_lonlat[0], loc_lonlat[1], loc_lonlat[2]
    try:
        skyimg_model = globaldiffuseskymodel(dattim, (lon, lat, h), freq,
                                             gs_model=gs_model, imsize=imsize)
    except ValueError:
        warnings.warn("Skipping GSM plot since frequency invalid")
    l, m = np.linspace(-1, 1, imsize), np.linspace(-1, 1, imsize)
    ll, mm = np.meshgrid(l, m)
    img_zero = np.zeros_like(skyimg_model, dtype=float)
    modality = 'model:'+gs_model
    _phaseref_ = (0, np.pi/2,'AZEL')
    integration = None
    correctpb = True
    fluxperbeam = False
    plotskyimage(ll, mm, (skyimg_model, img_zero, img_zero, img_zero),
                 'stokes', dattim, freq, stnid, integration, _phaseref_,
                 modality, pbcor=correctpb, maskhrz=False,
                 fluxperbeam=fluxperbeam, plot_title='Model image')
    plt.show()


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
    freq: float
        Center frequency in Hz
    ant_pos: array_like
        Position vectors of array elements
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
    vis = np.zeros((nr_ants, nr_ants), dtype=complex)
    k = 2 * np.pi * freq / c
    nn = n_from_lm(ll, mm)
    if imag_is_fd:
        fdd = skyimage / nn
    else:
        fdd = skyimage
    dll, dmm = imggrid_res(ll, mm)
    for ant_i in range(nr_ants):
        for ant_j in range(ant_i, nr_ants):
            u = pos_x[ant_i] - pos_x[ant_j]
            v = pos_y[ant_i] - pos_y[ant_j]
            w = pos_z[ant_i] - pos_z[ant_j]
            vis[ant_i, ant_j] = np.sum(
                fdd * np.exp(+1.0j * k * (ll * u + mm * v + (nn-1) * w))
                *dll*dmm)
    do_conj = True
    if do_conj:
        for ant_i in range(nr_ants):
            for ant_j in range(0, ant_i):
                vis[ant_i, ant_j] = np.conj(vis[ant_j, ant_i])
    return vis


def skymodel_visibility(t, stn_pos, freq, stn_antpos, gs_model, ant_model=''):
    """
    Create sky-model visibilities

    Parameters
    ----------
    t: datetime
        Epoch date-time
    stn_pos: array
        Station position on Earth as 3D ITRF vector.
    freq: float
        Frequency of observation in Hz.
    stn_antpos: array [3, N]
        Positions, 3D cartesian, of antennas in array.
    gs_model: str
        Global skymodel name.
    ant_model:
        Name of antenna model. '' for no antenna, or 'dual_dipole45' for dual
        dipoles tilted 45 deg.

    Returns
    -------
    skymod_vis: array
        Complex array of visibilities
    """
    #nr_ants = stn_pos.shape[1]
    imsize = 100  # ToDo: improve this expression?
    l = np.linspace(-1, 1, imsize)
    m = np.linspace(-1, 1, imsize)
    ll, mm = np.meshgrid(l, m)
    stn_pos_x, stn_pos_y, stn_pos_z \
        = stn_pos[0, 0], stn_pos[1, 0], stn_pos[2, 0]
    ccm = measures()
    ccm.doframe(ccm.position('ITRF', str(stn_pos_x) + 'm', str(stn_pos_y) + 'm',
                             str(stn_pos_z) + 'm'))
    lon, lat, h = ITRF2lonlat(stn_pos_x, stn_pos_y, stn_pos_z)
    img_S0 = globaldiffuseskymodel(t, (lon, lat, h), freq,
                                gs_model=gs_model,
                                imsize=imsize)
    ccm.doframe(ccm.epoch('UTC', t.isoformat('T')))
    phaseref_ccm = ccm.measure(
        ccm.direction('AZEL', '0.0rad', str(np.deg2rad(90)) + 'rad'), 'J2000')
    phaseref = (phaseref_ccm['m0']['value'], phaseref_ccm['m1']['value'],
                phaseref_ccm['refer'])
    uvw_sl = calc_uvw(t, phaseref, stn_pos, stn_antpos)
    if ant_model == 'dual_dipole45':
        (cov_xx, cov_xy, cov_yx, cov_yy) = dualdipole45_cov_patt(ll, mm)
        imag_xx = cov_xx * img_S0 / 2.0
        imag_xy = cov_xy * img_S0 / 2.0
        imag_yx = cov_yx * img_S0 / 2.0
        imag_yy = cov_yy * img_S0 / 2.0
        vis_xx = vcz(ll, mm, imag_xx, freq, uvw_sl, imag_is_fd=False)
        vis_xy = vcz(ll, mm, imag_xy, freq, uvw_sl, imag_is_fd=False)
        vis_yx = vcz(ll, mm, imag_yx, freq, uvw_sl, imag_is_fd=False)
        vis_yy = vcz(ll, mm, imag_yy, freq, uvw_sl, imag_is_fd=False)
    else:
        imag_xx = img_S0 / 2.0
        imag_yy = imag_xx
        vis_xx = vcz(ll, mm, imag_xx, freq, uvw_sl, imag_is_fd=False)
        vis_yy = vis_xx
        vis_xy = np.zeros_like(vis_xx)
        vis_yx = np.zeros_like(vis_xx)
    return vis_xx, vis_xy, vis_yx, vis_yy


def datetime64_to_datetime(dattim64):
    """\
    Convert datetime64 to python datetime

    Parameters
    ----------
    dattim64 : numpy.datetime64
        A numpy datetime64 date-time.

    Returns
    -------
    dattim : datetime.datetime
        A python datetime date-time.
    """
    dattim = datetime.datetime.utcfromtimestamp(
        (dattim64-np.datetime64('1970-01-01T00:00:00'))/np.timedelta64(1, 's'))
    return dattim


def create_vis_model(visdata, gs_model, ant_model=''):
    """
    Create visibility model

    Parameters
    ----------
    visdata: dict
        The visibilities themselves will computed overwriting existing
        values, this being the returned quantity.
        Visibility metadata used for observational setup.
    gs_model: str
        Name of global sky model to use
    ant_model: str
        Name of antenna model to use. Could be '' or 'dual_dipole45'. If empty
        str then no antenna model will be applied.

    Returns
    -------
    visdat: dict
        The model visibilities as a CVCfiles object

    Raises
    ------
    FileExitsError
        If the model filefolder already exists.
    """
    ldat_type = visdata['datatype']
    if ldat_type != 'acc' and ldat_type != 'xst':
        raise ValueError("ldat_type not set")
    # Look for visibility arrays
    arr_keys = list(filter(lambda k: k.startswith('arr_'), visdata.keys()))
    nrfiles = len(arr_keys)
    stn_pos, antpos_uv = layout_abs2rel(visdata['positions'])
    if visdata.get('stn_rot') is not None:
        antpos_uv = rot2uv(antpos_uv, visdata['stn_rot'].T)
    nrants = antpos_uv.shape[0]
    for filestep in range(nrfiles):
        print("Making model for file {}/{}".format(filestep, nrfiles))
        if nrfiles > 1:
            ts = visdata['delta_secs'][filestep]
        else:
            # There 2 cases for nrfiles==1: need to remove singleton
            ts = visdata['delta_secs'].squeeze()
        cvc_infile = np.empty((len(ts), 2, 2, nrants, nrants), dtype=complex)
        for sampstep in range(len(ts)):
            print("  sample {}/{}".format(sampstep, len(ts)), end='\r')
            if ldat_type == "xst":
                if visdata['frequencies'].ndim == 0:
                    freq = visdata['frequencies']  # Freq const. over xst file
            else:
                if nrfiles > 1:
                    freq = visdata['frequencies'][filestep][sampstep]
                else:
                    freq = visdata['frequencies'][sampstep]
            t = datetime64_to_datetime(visdata['start_datetime']
                                       + np.timedelta64(int(ts[sampstep]),'s'))
            if not gs_model.startswith('fid:'):
                vis_XX_mod, vis_XY_mod, vis_YX_mod, vis_YY_mod \
                    = skymodel_visibility(t, np.asmatrix(stn_pos).T, freq,
                                          np.asmatrix(visdata['positions']),
                                          gs_model, ant_model=ant_model)
            else:
                _fid, fid_typ = gs_model.split(':', 1)
                lam = c/freq
                stn_antpos_lambda = antpos_uv / lam
                if fid_typ == 'map':
                    img_I, ll, mm = fiducial_image(background=0.0)
                    # N.B. freq=c => k=2*np.pi since antpos in lambda units
                    # (nominally freq=c/lam)
                    vis_I_mod = vcz(ll, mm, img_I, c, stn_antpos_lambda,
                                    imag_is_fd=False)
                elif fid_typ == 'vis':
                    vis_I_mod = point_source_vis2d(stn_antpos_lambda, l0=0.0,
                                                   m0=0.0, amp=1.0)
                vis_XX_mod = vis_I_mod / 2.0
                vis_YY_mod = vis_I_mod / 2.0
                vis_XY_mod = np.zeros_like(vis_I_mod)
                vis_YX_mod = np.zeros_like(vis_I_mod)
            cvc_infile[sampstep, ...] = np.asarray(
                                           [[vis_XX_mod, vis_XY_mod],
                                            [vis_YX_mod, vis_YY_mod]])
            # Replace observed data with model data:
        visdata[arr_keys[filestep]] = cvc_infile
        print()
    # Note in ScanRecInfo about model for this dataset:
    visdata['model'] = gs_model + '+' + ant_model
    return visdata


def create_vis_model_ff(cvcpath, gs_model, ant_model=''):
    """
    Create visibility model file-folder based on measured data file-folder

    Parameters
    ----------
    cvcpath: str
        Path to measured visibilities
    gs_model: str
        Name of global sky model to use
    ant_model: str
        Name of antenna model to use. Could be '' or 'dual_dipole45'. If empty
        str then no antenna model will be applied.

    Returns
    -------
    cvcobj_mod: CVCfiles
        The model visibilities as a CVCfiles object

    Raises
    ------
    FileExitsError
        If the model filefolder already exists.
    """
    ldat_type = data_io.datafolder_type(cvcpath)
    if ldat_type != 'acc' and ldat_type != 'xst':
        raise ValueError("Not CVC data.")
    # Copy CVC folder within parent folder and add the tag "_cal_" in the name
    # right before ldat_type suffix:
    spltpath = cvcpath.split("_")
    cvcmodpath = "_".join(spltpath[:-1]) + "_mod_" + spltpath[-1]
    if os.path.exists(cvcmodpath):
        raise FileExistsError('Model folder already exists: '+cvcmodpath)
    shutil.copytree(cvcpath, cvcmodpath)
    # Read in as cvcobj model to be:
    cvcobj_mod = data_io.CVCfiles(cvcmodpath)
    nrfiles = cvcobj_mod.getnrfiles()
    antpos_uv = rot2uv(cvcobj_mod.stn_antpos, cvcobj_mod.stn_rot)
    # Loop over files in CVC folder:
    for filestep in range(nrfiles):
        if ldat_type == "xst":
            freq = cvcobj_mod.freqset[filestep][0]  # Freq const. over xst file
            sb, nz = modeparms.freq2sb(freq)
        else:
            sb = None  # Because this signals ACC data
        print("Making model for file {}/{}".format(filestep, nrfiles))
        freqs = cvcobj_mod.freqset[filestep]
        ts = cvcobj_mod.samptimeset[filestep]
        cvc_infile = np.empty((len(ts), cvcobj_mod.cvcdim1,
                               cvcobj_mod.cvcdim2), dtype=complex)
        for sampstep in range(len(ts)):
            print("  sample {}/{}".format(sampstep, len(ts)), end='\r')
            freq = freqs[sampstep]
            t = ts[sampstep]
            if not gs_model.startswith('fid:'):
                vis_XX_mod, vis_XY_mod, vis_YX_mod, vis_YY_mod \
                    = skymodel_visibility(t, cvcobj_mod.stn_pos, freq,
                                          cvcobj_mod.stn_antpos, gs_model,
                                          ant_model=ant_model)
            else:
                _fid, fid_typ = gs_model.split(':', 1)
                lam = c/freq
                stn_antpos_lambda = antpos_uv / lam
                if fid_typ == 'map':
                    img_I, ll, mm = fiducial_image(background=0.0)
                    # N.B. freq=c => k=2*np.pi since antpos in lambda units
                    # (nominally freq=c/lam)
                    vis_I_mod = vcz(ll, mm, img_I, c, stn_antpos_lambda,
                                    imag_is_fd=False)
                elif fid_typ == 'vis':
                    vis_I_mod = point_source_vis2d(stn_antpos_lambda, l0=0.0,
                                                   m0=0.0, amp=1.0)
                vis_XX_mod = vis_I_mod / 2.0
                vis_YY_mod = vis_I_mod / 2.0
                vis_XY_mod = np.zeros_like(vis_I_mod)
                vis_YX_mod = np.zeros_like(vis_I_mod)
            vis_flat_mod = cov_polidx2flat(np.asarray(
                                           [[vis_XX_mod, vis_XY_mod],
                                            [vis_YX_mod, vis_YY_mod]]))
            cvc_infile[sampstep, ...] = vis_flat_mod
            # Replace observed data with model data:
        cvcobj_mod[filestep] = cvc_infile
        print()
    # Note in ScanRecInfo about model for this dataset:
    cvcobj_mod.scanrecinfo.scanrecpath = cvcmodpath
    cvcobj_mod.scanrecinfo.set_model(gs_model+'+'+ant_model)
    return cvcobj_mod


def main_cli():
    """
    Compute a model based on LOFAR visibility datafile

    Resulting model can either be represented as a visibility datafile or as a
    sky image.
    """
    parser = argparse.ArgumentParser(
        description='Compute a model based on LOFAR datafiles.')
    parser.add_argument('rep', type=str,
                        help="Model representation. Choose: 'vis' or 'img'.")
    parser.add_argument('-n', '--filenr', type=int, default=0)
    parser.add_argument('-s', '--sampnr', type=int, default=0)
    parser.add_argument('-f', '--fluxpersterradian', action="store_true",
                        help="Normalize flux per sterradian")
    parser.add_argument('-g', '--gs_model', type=str, default='LFSM',
                        help="Name of global sky model")
    parser.add_argument('-a', '--ant_model', type=str, default='',
                        help="Name of antenna model")
    parser.add_argument('dataff',
                        help="""Path to CVC folder""")
    args = parser.parse_args()
    cvcobj = data_io.CVCfiles(args.dataff)
    if args.rep == 'vis':
        try:
            create_vis_model_ff(args.dataff, args.gs_model, args.ant_model)
        except FileExistsError as e:
            print(e, file=sys.stderr)
            print('Remove it to create a new vis model!', file=sys.stderr)
    elif args.rep == 'img':
        plot_gsm_for_obsdata(cvcobj, args.filenr, args.sampnr, args.gs_model)


if __name__ == "__main__":
    main_cli()
