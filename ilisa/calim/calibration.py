import numpy


def applycaltab_cvc(cvcunc, caltab, sb=None):
    """Apply a caltable to CVC data.

    Note
    ----
    Formula used is:
        G_ijk = g_ij*g^*_ik (no sum over i)
        V'_ijk = G_s0jk*V_ijk  (no sum over j,k & s0 is explicitly given)

    (Note that since function was designed for simplicity, it determines whether the
    CVC is an ACC (with all subbands) or XST (only one given subband) based on if
    no subband is given and its first index has size 512. There is a very small
    chance for the user to make a mistake by not setting 'sb' and the data happens to
    be a 512 samples XST.)

    Parameter
    ----------
    cvcunc : array [:,192,192]
        Uncalibrated CVC array, where 1st index is time.
    sb : int
        Subband in which the XST data was taken.
    caltab: array [512,192]
        Calibration table array to apply to xstunc.

    Returns
    -------
    cvccal: array [:,192,192]
        Calibrated XST array.
    """
    nrsbs, nrrcus0, nrrcus1 = cvcunc.shape
    if not sb and nrsbs != 512:
        # Cannot assume its an ACC
        raise ValueError("Must give sb for XST data.")
    gg = numpy.einsum('ij,ik->ijk', caltab, numpy.conj(caltab))
    if not sb and nrsbs == 512:
        # Assume it's an ACC
        g_apply = gg
    else:
        # It's an XST
        g_apply = gg[sb,:,:]
    cvccal = g_apply*cvcunc
    return cvccal