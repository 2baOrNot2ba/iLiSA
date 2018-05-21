#!/usr/bin/python
"""Provides support for handling LOFAR calibration tables."""
import os
import sys
import numpy
import datetime
import matplotlib.pyplot as plt
import ilisa.observations.stationcontrol as stationcontrol

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
    caltabfilename = 'CalTable_'+stnid[2:]+'_mode'+rcumode+'.dat'
    if obsdatestr is not None:
        # Need to determine most appropriate caltab to use
        caltabarchive = os.path.join(caltabdirstn, 'old_data/')
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
        print ctfullpath
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
    """
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
        fin = open(caltabfile)
        headline = fin.readline().rstrip()
        if headline != 'HeaderStart':
          raise Exception("{} is not a CalTable file.".format(caltabfile))
    except :
        raise Exception("Cannot use {} as CalTable file.".format(caltabfile))
    observation = {}
    calibration = {}
    comment = []
    while True:
        headline = fin.readline().rstrip()
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
    caltab = numpy.fromfile(fin, dtype=('c16', (stationcontrol.nrofrcus,)),
                                 count=stationcontrol.TotNrOfsb)
    fin.close()
    header = {'Observation': observation,
              'Calibration': calibration,
              'Comment':     comment}

    return caltab, header


def plotcaltab(caltab, header):
    """Plot a calibration table."""
    plt.subplot(211)
    plt.imshow(numpy.abs(caltab.T))
    plt.xlabel('subband [#]')
    plt.ylabel('RCU [#]')
    plt.title('Gain (Abs)\nRCU mode: '+header['Observation']['Mode'])
    plt.colorbar()
    plt.subplot(212)
    plt.imshow(numpy.rad2deg(numpy.angle(caltab.T)))
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
    G = numpy.einsum('ij,il->ijl', caltab, numpy.conj(caltab))
    acccal = G*accunc
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
    G = numpy.einsum('ij,il->ijl', caltab, numpy.conj(caltab))
    xstcal = G[sb,:,:]*xstunc
    return xstcal


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
    print header
    plotcaltab(caltab, header)

