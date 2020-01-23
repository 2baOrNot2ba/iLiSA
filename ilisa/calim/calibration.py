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


def cvcfolder_applycal(cvcpath, caltabpath):
    """Apply a calibration table file to a CVC folder.
    This creates a copy of the folder pointed to by cvcpath renamed with a '_cal_'
    before the ldat suffix. Then it applies the calibration table contained in the
    caltab file pointed to by caltabpath, to the CVC dataset and copies the caltab
    file used in the calibrated CVC folder.

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
    # Copy CVC folder within parent folder and add the tag "_cal_" in the name right
    # before ldat_type suffix:
    spltpath = cvcpath.split("_")
    cvccalpath = "_".join(spltpath[:-1]) + "_cal_" + spltpath[-1]
    shutil.copytree(cvcpath, cvccalpath)
    # Read in cvcobj:
    cvcobj_cal = dataIO.CVCfiles(cvccalpath)
    nrfiles = cvcobj_cal.getnrfiles()
    # Loop over files in CVC folder:
    for filestep in range(nrfiles):
        if ldat_type == "xst":
            freq = cvcobj_cal.freqset[filestep][0]  # Freq constant over xst file
            sb, nz = modeparms.freq2sb(freq)
        else:
            sb = None  # Because this signals ACC data
        # Get actual covariance cubes:
        cvcdata_unc = cvcobj_cal.getdata(filestep)
        # Apply calibration
        cvcdata = ilisa.calim.calibration.applycaltab_cvc(cvcdata_unc, caltab, sb)
        # Replace uncalibrated data with calibrated:
        cvcobj_cal.dataset[filestep] = cvcdata
    cvcobj_cal.save_ldat()
    # Note in ScanRecInfo about calibrating this dataset:
    cvcobj_cal.scanrecinfo.set_postcalibration(caltabpath, cvccalpath)