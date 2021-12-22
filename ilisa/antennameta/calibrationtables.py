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
import numpy
import datetime
import argparse
import warnings
import matplotlib.pyplot as plt

import ilisa.operations.modeparms

__version__ = '0.2'
CALTABDIRROOT = os.path.join(os.path.dirname(__file__), 'share/CalTables/')


def _default_caltab_filename(stnid, rcumode):
    """Get default filename of a LOFAR station calibration table.
    """
    caltabfilename = 'CalTable_'+stnid[2:]+'_mode'+rcumode+'.dat'
    return caltabfilename


def _findcaltabpath(rcumode, stnid, obsdatestr=None):
    """Find appropriate caltab file based on rcumode stnid and observation
    date.
    """
    def adddatestr(_cthist, _ctfilename, _ctheader):
        use_date = 'Calibration'
        if use_date == 'Observation':
            datestr = _ctheader['Observation']['Date']
            datestr = datestr[:8]+'T'+datestr[8:]+'00'
        else:
            datestr = _ctheader['Calibration']['Date']
            datestr = datestr+'T120000'
        _cthist[datestr] = _ctfilename
    caltabdirroot = CALTABDIRROOT
    caltabdirstn = os.path.join(caltabdirroot, stnid)
    # In practice rcumode 4 uses rcumode 3's caltab
    # while rcumode 6 uses rcumode 5's caltab.
    # So map accordingly
    if rcumode == '4':
        rcumode = '3'
        warnings.warn("Using caltab for rcumode 3 instead.")
    if rcumode == '6':
        rcumode = '5'
        warnings.warn("Using caltab for rcumode 5 instead.")
    caltabfilename = _default_caltab_filename(stnid, rcumode)
    caltabarchive = os.path.join(caltabdirstn, 'old_data/')
    caltabarchive_exists = os.path.isdir(caltabarchive)
    if obsdatestr is not None and caltabarchive_exists:
        # Need to determine most appropriate caltab to use
        caltabstndatestr = os.listdir(caltabarchive)
        cthist = {}
        for ctdir in caltabstndatestr:
            ctfullpath = os.path.join(caltabarchive, ctdir, caltabfilename)
            try:
                (cttable, ctheader) = readcaltab(ctfullpath)
            except RuntimeError:
                continue
            adddatestr(cthist, ctfullpath, ctheader)
        # Get latest:
        ctfullpath = os.path.join(caltabdirstn, 'data', caltabfilename)
        (caltab_latest, ctheader_latest) = readcaltab(ctfullpath)
        adddatestr(cthist, ctfullpath, ctheader_latest)
        obsdate = datetime.datetime.strptime(obsdatestr, "%Y-%m-%d")
        cthistdates = list(cthist.keys())
        caltabdates = [datetime.datetime.strptime(d, "%Y%m%dT%H%M%S")
                       for d in cthistdates]
        difobscal = [abs(obsdate - d) for d in caltabdates]
        caltabpath = cthist[cthistdates[difobscal.index(min(difobscal))]]
    else:
        if obsdatestr is not None and not caltabarchive_exists:
            warnings.warn('No archived caltables; getting default caltable.')
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
        Date of the observation. The format is 'YYYY-mm-dd'.
        If None (default), will use latest.

    Returns
    -------
    caltab : (512, 192) array
        The calibration (gains) table: 0-axis=subband, 1-axis=rcunr
    header : dict
        The header of the calibration table file.
    """
    rcumode = str(rcumode)
    caltabpath = _findcaltabpath(rcumode, stnid, obsdatestr)
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
            raise RuntimeError("{} is not a CalTable file.".format(caltabfile))
    except (OSError, RuntimeError):
        raise RuntimeError("Cannot use {} as CalTable file.".format(caltabfile))
    observation = {}
    calibration = {}
    comment = []
    while True:
        headline = fin.readline().decode('UTF-8').rstrip()
        if headline == 'HeaderStop':
            break
        (caltabheadmark, caltableheadline) = headline.split('.', 1)
        var, val = caltableheadline.split('=', 1)
        var = var.rstrip(' ')
        val = val.lstrip(' ')
        if var == 'Comment':
            comment.append(val)
        else:
            cat, key = var.split('.', 1)
            if cat == 'Observation':
                observation[key] = val
            elif cat == 'Calibration':
                calibration[key] = val
    # Assume that max nr of subbands (ilisa.operations.modeparms.TotNrOfsb)
    # is always 512, while nr of RCUs (ilisa.operations.modeparms.nrofrcus)
    # may differ:
    caltab = numpy.fromfile(fin, dtype='c16').reshape(
        (ilisa.monitorcontrol.modeparms.TotNrOfsb, -1))
    # nrrcus = caltab.shape[1]
    fin.close()
    header = {'Observation': observation,
              'Calibration': calibration,
              'Comment':     comment}
    return caltab, header


def writecaltab(filename, caltab, observation, calibration, comments):
    """Write a calibration table to file. Inverse of readcaltab()."""
    fout = open(filename, 'w')
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


def initcaltab(stnid, rcumode):
    """Initialize a calibration table file.
    Sets gains for all rcus and all subbands to 1.
    """
    # Default filename
    ctfname = _default_caltab_filename(stnid, rcumode)
    # Init caltab
    nrrcus = ilisa.monitorcontrol.modeparms.nrofrcus
    nrsbs = ilisa.monitorcontrol.modeparms.TotNrOfsb
    caltab = numpy.ones((nrrcus, nrsbs), dtype='c16')
    # Init header
    band = ilisa.monitorcontrol.modeparms.rcumode2band(rcumode)
    antset = ilisa.monitorcontrol.modeparms.band2antset(band)
    header_observation = {'Station': stnid, 'Mode': rcumode,
                          'AntennaSet': antset, 'Band': band, 'Source': '',
                          'Date': ''}
    user = os.environ.get('USER')
    now = datetime.datetime.utcnow()
    header_calibration = {'Version': '0', 'Name': user,
                          'Date': now.strftime('YYYYmmdd'),
                          'PPSDelay': '[]'}
    header_comment = ['Initial CalTable created with iLiSA']
    writecaltab(ctfname, caltab, header_observation, header_calibration,
                header_comment)


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
    delay_n_phasercu = numpy.dot(betainv, argsrcusb)
    reconstruct = True
    if reconstruct:
        caltabr = numpy.outer(numpy.ones((nrsbs,)), ampsrcu) * numpy.exp(1j*(
                delay_n_phasercu[1, :]+numpy.outer(fs, delay_n_phasercu[0, :])))
        assert numpy.allclose(caltab, caltabr)
    return ampsrcu, delay_n_phasercu


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


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers(dest='subparser_name',
                                       help='sub-command help')
    # 'show' command
    parser_show = subparsers.add_parser('show',
                                        help='Show contents of caltab file.')
    parser_show.add_argument('caltab_path',
                             help="LOFAR calibration table file")
    # 'find' command
    parser_find = subparsers.add_parser('find',
                                        help="""Find caltab file for
                                         rcumode, stnid, obsdate""")
    parser_find.add_argument('rcumode', help="Band of observation.")
    parser_find.add_argument('stnid', help="Station ID of observation.")
    parser_find.add_argument('date', nargs='?', default=None,
                             help="Date of observation. Format: YYYY-mm-dd")
    # 'create' command
    parser_create = subparsers.add_parser('create',
                                          help='Create a caltab file.')
    parser_create.add_argument('rcumode',
                               help="Rcumode of calibration table file")
    parser_create.add_argument('stnid', help="Station ID of observation.")

    args = parser.parse_args()
    if args.subparser_name == 'show':
        caltab_cur, header_cur = readcaltab(args.caltab_path)
        print(header_cur)
        print(caltab_cur)
        plotcaltab(caltab_cur, header_cur)
    elif args.subparser_name == 'find':
        print(_findcaltabpath(args.rcumode, args.stnid, args.date))
    elif args.subparser_name == 'create':
        initcaltab(args.stnid, args.rcumode)
