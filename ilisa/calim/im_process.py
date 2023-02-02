"""Image process module"""
import numpy

def lm_grid(nrpix=128,l_ext=1.0, m_ext=1.0):
    """
    Create a grid of direction cosine coordinates l and m.

    Parameters
    ----------
    nrpix : int
        Number of pixels.
    l_ext : float
        Extent of l.
    m_ext : float
        Extent of m.

    Returns
    -------
    ll : array_like
        2D grid of ll coordinates.
    mm :
        2D grid of mm coordinates.
    """
    l, m = (numpy.linspace(-l_ext, l_ext, nrpix),
            numpy.linspace(-m_ext, m_ext, nrpix))
    ll, mm = numpy.meshgrid(l, m)
    return ll, mm


def dynamic_range(im):
    """
    Simplistic dynamic range
    """
    max = numpy.max(im)
    mean = numpy.std(im)
    return max/mean


def split_horizon(im, hor_width=0.0, fill=0.0):
    """
    Split image(s) into regions: above horizon, below horizon and horizon.

    The horizon is the circle where l^2+m^2==1.

    Parameters
    ----------
    im : array_like
        Image(s) to be processed. If ndim=3, the 0 axis is assumed to be the
        stack index of the images. The images are sssumed to be 
        over l=-1,...,+1 and m=-1,...,+1. Furthermore the image is assumed to
        be gridded evenly, i.e. the images are square arrays.
    hor_width : float
        Width of horizon in l,m units. Horizon region is centered on the
        horizon and extends hor_width/2 inwards and hor_width/2 outwards.
    fill : float
        Value to fill excluded region with. 
    
    Returns
    -------
    above : array_like
        Region above (inside) the horizon.
    below : array_like
        Region below (outside) the horizon.
    horiz : array_like
        Horizon region.
    
    Examples
    --------
    >>> import numpy
    >>> import ilisa.calim.imageprocess as imageprocess
    >>> im = numpy.ones((4,4))
    >>> above, below, horiz = imageprocess.split_horizon(im, hor_width=0.2)
    >>> above
    [[0. 0. 0. 0.]
     [0. 1. 1. 0.]
     [0. 1. 1. 0.]
     [0. 0. 0. 0.]]
    >>> below
    [[1. 0. 0. 1.]
     [0. 0. 0. 0.]
     [0. 0. 0. 0.]
     [1. 0. 0. 1.]]
    >>> horiz
    [[0. 1. 1. 0.]
     [1. 0. 0. 1.]
     [1. 0. 0. 1.]
     [0. 1. 1. 0.]]
    """
    if im.ndim == 2:
        im = im[numpy.newaxis, :, :]
    elif im.ndim != 3:
        raise TypeError("Only 2 or 3 dimensional images allowed.")
    nrimg = im.shape[0]
    nrpix = im.shape[-1]
    ll, mm = lm_grid(nrpix=nrpix)
    lm_norm = numpy.sqrt(ll**2+mm**2)
    lm_above = lm_norm < 1.0-hor_width/2.0
    lm_below = lm_norm > 1.0+hor_width/2.0
    lm_horiz = numpy.logical_and(numpy.logical_not(lm_above),
                                 numpy.logical_not(lm_below))
    above = numpy.copy(im)
    below = numpy.copy(im)
    horiz = numpy.copy(im)
    for imgnr in range(nrimg):
        numpy.place(above[imgnr, :, :], numpy.logical_not(lm_above), fill)
        numpy.place(below[imgnr, :, :], numpy.logical_not(lm_below), fill)
        numpy.place(horiz[imgnr, :, :], numpy.logical_not(lm_horiz), fill)
    above=above.squeeze()
    below=below.squeeze()
    horiz=horiz.squeeze()
    return above, below, horiz


def n_from_lm(ll,mm):
    """
    Compute n of an l,m grid

    Parameters
    ----------
    ll: array
        l grid values
    mm: array
        m grid values

    Returns
    -------
    nn : array
        Grid values for n=sqrt(1-l^2-m^2). Beyond horizon l^2+m^2=1, n is set
        to unity.
    """
    lm_r2 = (ll**2+mm**2).astype(numpy.complex)
    nn = numpy.sqrt(1-lm_r2)
    nn[lm_r2 >= 1.0] = 1.0
    return nn


def brightness_sr_2_lm(map_sr, ll, mm):
    """
    Convert brightness per sterradian distribution to brightness per l,m

    Brightness distribution is nominally given in units flux per sterradian.
    For the integration in the van Cittert-Zernike relation, it is more
    convenient to use direction cosines l,m, i.e.

        int map_sr dSR = int map_lm / sqrt(1-l^2-m^2) dl dm

    The conversion in this case thus: map_sr dSR = int map_lm / sqrt(1-l^2-m^2)

    Parameters
    ----------
    map_sr: array
        Brightness distribution per sterradian
    ll, mm: array
        l,m grid of brightness map

    Returns
    -------
    map_lm: array
        Brightness distribution per l,m pixel
    """
    nn = n_from_lm(ll, mm)
    map_lm = map_sr / nn
    return map_lm
