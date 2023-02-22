import argparse
import shutil
import warnings

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
from .imaging import plotskyimage, imggrid_res
from .visibilities import cov_polidx2flat, calc_uvw, rot2uv
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
    Generate hemisphere of global diffuse sky model (GDSM) over a position and for
    given datetime and freq.
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
            plot_gsm(t, cvcobj.stn_pos, freq, gs_model, imsize, stnid)


def plot_gsm(dattim, stn_pos, freq, gs_model='LFSM', imsize=200,
             stnid='Unknown'):
    lon, lat, h = ITRF2lonlat(stn_pos[0,0], stn_pos[1,0], stn_pos[2,0])
    try:
        skyimg_model = globaldiffuseskymodel(dattim, (lon, lat, h), freq,
                                             gs_model=gs_model, imsize=imsize)
    except ValueError:
        warnings.warn("Skipping GSM plot since frequency invalid")
    l, m = np.linspace(-1, 1, imsize), np.linspace(-1, 1, imsize)
    ll, mm = np.meshgrid(l, m)
    img_zero = np.zeros_like(skyimg_model, dtype=float)
    calibrated = 'N.A.'
    _phaseref_ = (0,0,'AZEL')
    integration = 0
    correctpb = True
    fluxperbeam = True
    plotskyimage(ll, mm, (skyimg_model, img_zero, img_zero, img_zero),
                 'stokes', dattim, freq, stnid, integration, _phaseref_,
                 calibrated, pbcor=correctpb, maskhrz=False,
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


def fiducial_visibility(stn_antpos_lambda, background=1.0):
    """
    Generate a fiducial visibility

    Parameters
    ----------
    stn_antpos: array
        Positions of antennas in station.
    background: float
        Background flux density level.

    Returns
    -------
    vis: array
        The fiducial visibility.
    """
    imsize = 100  #21, 43
    l = np.linspace(-1, 1, imsize)
    m = np.linspace(-1, 1, imsize)
    ll, mm = np.meshgrid(l, m)
    img = background*np.ones_like(ll)
    make_sq_reg = True
    if make_sq_reg:
        sqregion = np.zeros_like(ll)
        sqregion[(-0.0<ll) & (ll<0.2) & (-0.9<mm) & (mm<-0.7)] = 2.0
        img += sqregion
    # Zero beyond horizon
    img[ll**2+mm**2>=1.0] = 0.0
    # N.B. freq=c/(2*np.pi) => k=1
    vis_mod = vcz(ll, mm, img, c/(2*np.pi), stn_antpos_lambda, imag_is_fd=True)
    return vis_mod


def model_visibility(t, stn_pos, freq, stn_antpos, gs_model):
    """
    Create visibility model

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

    Returns
    -------
    model_visibility: array
        Complex array of visibilities.
    """
    nr_ants = stn_pos.shape[1]
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
    img = globaldiffuseskymodel(t, (lon, lat, h), freq,
                                gs_model=gs_model,
                                imsize=imsize)
    ccm.doframe(ccm.epoch('UTC', t.isoformat('T')))
    phaseref_ccm = ccm.measure(
        ccm.direction('AZEL', '0.0rad', str(np.deg2rad(90)) + 'rad'), 'J2000')
    phaseref = (phaseref_ccm['m0']['value'], phaseref_ccm['m1']['value'],
                phaseref_ccm['refer'])
    uvw_sl = calc_uvw(t, phaseref, stn_pos, stn_antpos)
    vis_mod = vcz(ll, mm, img, freq, uvw_sl, imag_is_fd=False)
    return vis_mod


def create_vis_model_ff(cvcpath, gs_model):
    """
    Create visibility model file-folder based on measured data file-folder

    Parameters
    ----------
    cvcpath: str
        Path to measured visibilities
    gs_model: str
        Name of global sky model to use

    Returns
    -------
    cvcobj_mod: CVCfiles
        The model visibilities as a CVCfiles object
    """
    ldat_type = data_io.datafolder_type(cvcpath)
    if ldat_type != 'acc' and ldat_type != 'xst':
        raise ValueError("Not CVC data.")
    # Copy CVC folder within parent folder and add the tag "_cal_" in the name
    # right before ldat_type suffix:
    spltpath = cvcpath.split("_")
    cvcmodpath = "_".join(spltpath[:-1]) + "_mod_" + spltpath[-1]
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
            freq = freqs[sampstep]
            t = ts[sampstep]
            if not gs_model.startswith('fid:'):
                vis_I_mod = model_visibility(t, cvcobj_mod.stn_pos, freq,
                                             cvcobj_mod.stn_antpos, gs_model)
            else:
                _fid, fid_mod = gs_model.split(':', 1)
                lam = c/freq
                vis_I_mod = fiducial_visibility(antpos_uv/lam*2*np.pi,
                                                background=0.)
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
    # Note in ScanRecInfo about model for this dataset:
    cvcobj_mod.scanrecinfo.scanrecpath = cvcmodpath
    cvcobj_mod.scanrecinfo.set_model(gs_model)
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
    parser.add_argument('dataff',
                        help="""Path to CVC folder""")
    args = parser.parse_args()
    cvcobj = data_io.CVCfiles(args.dataff)
    if args.rep == 'vis':
        create_vis_model_ff(args.dataff, args.gs_model)
    elif args.rep == 'img':
        plot_gsm_for_obsdata(cvcobj, args.filenr, args.sampnr, args.gs_model)


if __name__ == "__main__":
    main_cli()
