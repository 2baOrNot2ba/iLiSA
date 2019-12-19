#!/usr/bin/python
"""Provides support for handling LOFAR calibration tables.

LOFAR caltables are files on the station (and included in the './share/CalTables/'
directory) that contain estimated corrective complex gains for each subband and each rcu
per band. They also contain metadata about e.g. when they made.

The format of a caltable is an complex array indexed by [subband,rcunr].
To apply them to data, one uses the generic formula:
    V' = g.*V.*g^H
    g = C[subband,:]
where the rcunr index is implicit, C is the full caltable, g is a vector of gains for one
subband and V is a visibility matrix.
"""
import os
import sys
import numpy
import datetime
import matplotlib.pyplot as plt

import ilisa.observations.modeparms

__version__ = '0.1'
CALTABDIRROOT = os.path.join(os.path.dirname(__file__), 'share/CalTables/')

def findcaltabpath(rcumode, stnid, obsdatestr=None):
    """Find appropriate caltab file based on rcumode stnid and observation
    date.
    """
    def adddatestr(cthist, ctfilename, ct_header):
        use_date = 'Calibration'
        if use_date == 'Observation':
            datestr = ct_header['Observation']['Date']
            datestr = datestr[:8]+'T'+datestr[8:]+'00'
        else:
            datestr = ct_header['Calibration']['Date']
            datestr = datestr+'T120000'
        cthist[datestr] = ctfilename
    caltabdirroot = CALTABDIRROOT
    caltabdirstn = os.path.join(caltabdirroot, stnid)
    # In practice rcumode 4 uses rcumode 3's caltab
    # while rcumode 6 uses rcumode 5's caltab.
    # So map accordingly
    if rcumode == '4':
        rcumode = '3'
        print("Warning: using caltab for rcumode 3 instead.")
    if rcumode == '6':
        rcumode = '5'
        print("Warning: using caltab for rcumode 5 instead.")
    caltabfilename = 'CalTable_'+stnid[2:]+'_mode'+rcumode+'.dat'
    caltabarchive = os.path.join(caltabdirstn, 'old_data/')
    if obsdatestr is not None and os.path.isdir(caltabarchive):
        # Need to determine most appropriate caltab to use
        caltabstndatestr = os.listdir(caltabarchive)
        cthist = {}
        for ctdir in caltabstndatestr:
            ctfullpath = caltabarchive+ctdir+'/'+caltabfilename
            try:
                (cttable, ctheader) = readcaltab(ctfullpath)
            except:
                continue
            adddatestr(cthist, ctfullpath, ctheader)
        # Get latest:
        ctfullpath = caltabdirstn+'/data/'+caltabfilename
        print(ctfullpath)
        (caltab_latest, ctheader_latest) = readcaltab(ctfullpath)
        adddatestr(cthist, ctfullpath, ctheader_latest)
        obsdate = datetime.datetime.strptime(obsdatestr, "%Y%m%dT%H%M%S")
        cthistdates = cthist.keys()
        caltabdates = [datetime.datetime.strptime(d,"%Y%m%dT%H%M%S") for d in cthistdates]
        difobscal = [abs(obsdate - d) for d in caltabdates]
        caltabpath = cthist[cthistdates[difobscal.index(min(difobscal))]]
    else:
        caltabpath = os.path.join(caltabdirstn, 'data', caltabfilename)
    return caltabpath


def getcaltab(rcumode, stnid, obsdatestr=None):
    """Return appropriate calibration table content based on rcumode stnid and
    observation date.

    Parameters
    ----------
    rcumode : str, int
        The rcumode. Can be 3, 4, 5, 6, 7, either as string or integer.
    stnid : str
        The station id, e.g. 'SE607'.
    obsdatestr : str
        Date of the observation. If None (default), will use latest.

    Returns
    -------
    caltab : (512, 192) array
        The calibration (gains) table: 0-axis=subband, 1-axis=rcunr
    header : dict
        The header of the calibration table file.
    """
    rcumode = str(rcumode)
    caltabpath = findcaltabpath(rcumode, stnid, obsdatestr)
    (caltab, header) = readcaltab(caltabpath)
    return caltab, header


def readcaltab(caltabfile):
    """Readin a calibration table file by name.
    
    Parameters
    ----------
    caltabfile : str
        The name of the calibration table file.
    
    Returns
    -------
    caltab : (512, 192) array
        The calibration (gains) table: 0-axis=subband, 1-axis=rcunr
    header : dict
        The header of the calibration table file.
    """
    
    try:
        fin = open(caltabfile, 'rb')
        headline = fin.readline().decode('UTF-8').rstrip()
        if headline != 'HeaderStart':
          raise Exception("{} is not a CalTable file.".format(caltabfile))
    except :
        raise Exception("Cannot use {} as CalTable file.".format(caltabfile))
    observation = {}
    calibration = {}
    comment = []
    while True:
        headline = fin.readline().decode('UTF-8').rstrip()
        if headline == 'HeaderStop': break
        (caltabheadmark, caltableheadline) = headline.split('.',1)
        var, val = caltableheadline.split('=', 1)
        var = var.rstrip(' ')
        val = val.lstrip(' ')
        if var == 'Comment':
            comment.append(val)
        else:
            cat, key = var.split('.',1)
            if cat == 'Observation':
                observation[key] = val
            elif cat == 'Calibration':
                calibration[key] = val
    caltab = numpy.fromfile(fin, dtype=('c16', (ilisa.observations.modeparms.nrofrcus,)),
                            count=ilisa.observations.modeparms.TotNrOfsb)
    fin.close()
    header = {'Observation': observation,
              'Calibration': calibration,
              'Comment':     comment}

    return caltab, header


def plotcaltab(caltab, header):
    """Plot a calibration table."""
    plt.subplot(211)
    plt.pcolormesh(numpy.abs(caltab.T))
    plt.xlabel('subband [#]')
    plt.ylabel('RCU [#]')
    plt.title('Gain (Abs)\nRCU mode: '+header['Observation']['Mode'])
    plt.colorbar()
    plt.subplot(212)
    plt.pcolormesh(numpy.rad2deg(numpy.angle(caltab.T)))
    plt.xlabel('subband [#]')
    plt.ylabel('RCU [#]')
    plt.title('Gain (Phase)\nRCU mode: '+header['Observation']['Mode'])
    plt.colorbar()
    plt.show()


def writecaltab(caltabfile, caltab, observation, calibration, comments):
    """Write a calibration table to file. Inverse of readcaltab()."""
    fout = open(caltabfile,'w')
    fout.write('HeaderStart\n')
    for k in observation.keys():
        fout.write('CalTableHeader.Observation.'+k+' = '+observation[k]+'\n')
    for k in calibration.keys():
        fout.write('CalTableHeader.Calibration.'+k+' = '+calibration[k]+'\n')
    for c in comments:
        fout.write('CalTableHeader.Comment = '+c+'\n')
    fout.write('HeaderStop\n')
    caltab.tofile(fout)
    fout.close()


def calibrateACC(accunc, rcumode, stnid, obsdate, docalibrate=True):
    """Calibrate ACC data by using appropriate caltable."""
    obsdatestr = obsdate.strftime("%Y%m%dT%H%M%S")
    if docalibrate:
        caltabtable, caltabheader = getcaltab(rcumode, stnid, obsdatestr)
        acc = calibrateACCwithtab(accunc, caltabtable)
    else:
        caltabtable, caltabheader = None, None
        acc = accunc
    return acc, (caltabtable, caltabheader)


def calibrateACCwithtab(accunc, caltab):
    """Apply a caltable to ACC data.

    Note
    ----
    Formula used is:
        G_ijk = g_ij*g^*_ik (no sum over i)
        V'_ijk = G_ijk*V_ijk  (no sum over i,j,k)

    Parameters
    ----------
    accunc : array [512,192,192]
        Uncalibrated ACC array.
    caltab: array [512,192]
        Calibration table array to apply to accunc.

    Returns
    -------
    acccal: array [512,192,192]
        Calibrated ACC array.
    """
    gg = numpy.einsum('ij,ik->ijk', caltab, numpy.conj(caltab))
    acccal = gg*accunc
    return acccal


def calibrateXST(xstunc, sb, rcumode, stnid, obsdate, docalibrate=True):
    """Calibrate XST data by using appropriate caltable."""
    obsdatestr = obsdate.strftime("%Y%m%dT%H%M%S")
    if docalibrate:
        caltabtable, caltabheader = getcaltab(rcumode, stnid, obsdatestr)
        xst = calibrateXSTwithtab(xstunc, sb, caltabtable)
    else:
        caltabtable, caltabheader = None, None
        xst = xstunc
    return xst, (caltabtable, caltabheader)


def calibrateXSTwithtab(xstunc, sb, caltab):
    """Apply a caltable to XST data.

    Note
    ----
    Formula used is:
        G_ijk = g_ij*g^*_ik (no sum over i)
        V'_ijk = G_s0jk*V_ijk  (no sum over j,k & s0 is explicitly given)

    Parameter
    ----------
    xstunc : array [:,192,192]
        Uncalibrated XST array, where 1st index is time.
    sb : int
        Subband in which the XST data was taken.
    caltab: array [512,192]
        Calibration table array to apply to xstunc.

    Returns
    -------
    xstcal: array [:,192,192]
        Calibrated XST array.
    """
    gg = numpy.einsum('ij,ik->ijk', caltab, numpy.conj(caltab))
    xstcal = gg[sb,:,:]*xstunc
    return xstcal


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


def getelemgainampdel(caltab):
    """Get the amplitudes and delays for the element (rcu) gains in the given
    calibration table."""
    (nrsbs, nrrcus) = caltab.shape
    ampsrcusb = numpy.abs(caltab)
    ampsrcu = numpy.mean(ampsrcusb, axis=0)
    argsrcusb = numpy.unwrap(numpy.angle(caltab))
    fs = numpy.arange(nrsbs)
    o = numpy.ones((nrsbs,))
    beta = numpy.asarray([fs, o]).T
    betainv = numpy.linalg.pinv(beta)
    delayNphasercu = numpy.dot(betainv,argsrcusb)
    reconstruct = True
    if reconstruct:
        caltabr=numpy.outer(numpy.ones((nrsbs,)), ampsrcu) * numpy.exp(1j*(
                        delayNphasercu[1,:]+numpy.outer(fs, delayNphasercu[0,:])))
        assert numpy.allclose(caltab,caltabr)
    return ampsrcu, delayNphasercu


if __name__=="__main__":
    caltab, header = readcaltab(sys.argv[1])
    print(header)
    print(caltab)
    plotcaltab(caltab, header)

