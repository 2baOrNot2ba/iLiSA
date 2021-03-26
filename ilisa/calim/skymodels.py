import numpy as np
from matplotlib import pylab
import healpy as hp
from pygdsm.pygsm import GlobalSkyModel
from pygdsm.pygsm2016 import GlobalSkyModel2016
from pygdsm.lfsm import LowFrequencySkyModel
from pygdsm.base_observer import BaseObserver


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
    return img
