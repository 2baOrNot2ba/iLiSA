import datetime
import os
import warnings

import numpy

import ilisa
from ilisa.operations import modeparms

def datafolder_type(datafolderpath):
    """Determine what type of LOFAR local mode recording the datafolder is.

    Parameters
    ----------
    datafolderpath : str
        Path to LOFAR local mode datafolder.

    Returns
    -------
    lofardatatype : str
        LOFAR local mode datatype. Can be 'acc', 'bfs', 'bst', 'sst', 'xst'.
        Returns None if datafolderpath is not one of the recognized types.
    """
    datafolder = os.path.basename(os.path.normpath(datafolderpath))
    suf = datafolder.split('_')[-1]
    if suf == 'acc':
        datatype = 'acc'
    elif suf == 'bfs':
        datatype = 'bfs'
    elif suf.endswith('bst'):
        # Allows for extended bst eg '12bst' ie with less than 4 lanes
        datatype = 'bst'
    elif suf == 'bstc':
        datatype = 'bstc'
    elif suf == 'sst':
        datatype = 'sst'
    elif suf == 'xst':
        datatype = 'xst'
    else:
        datatype = None
    return datatype


def obsinfo2filefolder(obsinfo):
    """\
    Convert obsinfo dict to filefolder name

    Parameters
    ----------
    obsinfo: dict
        Some parameters used for observation. Has the following keys (optional
        keys marked with '*'):
            'duration_scan':
            'filenametime':
            'integration':
            'spw':
            'ldat_type' | 'datatype':
            'pointing':
            'subbands':
            'station_id':
            *'frequencies':
            *'antennaset':
            *'cal':
            *'model':

    Returns
    -------
    filefoldername : str
        Meta-data formatted name of file-folder.
        Name format is
            <station_id(6c)>[<ANTSET>]_<filenametime>_spw<rcumodes>_sb<subbands>\
            _int<integration>_dur<duration_scan>[_dir<pointing>][_cal|_mod]\
            _<ldat_type>
            where <ANTSET> := <LBA|HBA%<conf1>[%<conf2>]>
    """
    if obsinfo.get('datatype', None):
        ldat_type = obsinfo['datatype']
    else:
        ldat_type = obsinfo.get('ldat_type')
    filefoldername = obsinfo['station_id']
    antennaset = obsinfo.get('antennaset')
    if antennaset:
        antennaset = antennaset.replace('_', '-')
        filefoldername += antennaset
    filefoldername += '_' + obsinfo['filenametime']

    spwstr = ''
    if obsinfo['spw']:
        spwstr = \
            ''.join([str(spw) for spw in obsinfo['spw']])
    filefoldername += '_spw' + spwstr
    if ldat_type != 'sst' and ldat_type != 'acc':
        if obsinfo['subbands'] != [] and obsinfo['subbands'] != '':
            filefoldername += "_sb"
            filefoldername += modeparms.seqlists2slicestr(obsinfo['subbands'])
    if 'integration' in obsinfo and obsinfo['integration']:
        filefoldername += '_int' + str(obsinfo['integration'])
    if 'duration_scan' in obsinfo:
        filefoldername += '_dur' + str(int(obsinfo['duration_scan']))
    if ldat_type != 'sst':
        if obsinfo['pointing'] is not None and str(obsinfo['pointing']) != '':
            filefoldername += '_dir' + str(obsinfo['pointing'])
        else:
            filefoldername += '_dir,,'
    cal = obsinfo.get('cal', None)
    if cal:
        if ldat_type == 'acc' or ldat_type == 'xst':
            filefoldername += '_cal'
        else:
            warnings.warn('Only ACC and XST can be calibrated.')
    else:
        model = obsinfo.get('model', None)
        if model:
            filefoldername += '_mod'
    # filefoldername += "_" + obsinfo['source']
    # ldat_type extension
    filefoldername += "_" + ldat_type
    return filefoldername


def filefolder2obsinfo(filefolderpath):
    """\
    Parse filefolder name and return an obsinfo

    For description of filefolder naming and obsinfo see
    obsinfo2filefolder().

    Parameters
    ----------
    filefolderpath : str
        Path to filefolder.

    Returns
    -------
    obsinfo : dict
        Dict of metadata corresponding to file-folder.
    """
    obsinfo = {}
    filefolderpath = os.path.normpath(filefolderpath)
    filefoldername = os.path.basename(filefolderpath)
    ldat_type = datafolder_type(filefoldername)
    # Format:
    # stnid_Ymd_HMS_spwstr_intstr_durstr_dirstr_[cal*]_acc
    # stnidantset_Ymd_HMS_spwstr_intstr_durstr_dirstr_[cal*]_acc
    # stnid?_Ymd_HMS_rcustr_sbstr_intstr_durstr_dirstr_bst
    # stnid?_Ymd_HMS_rcustr_sbstr_intstr_durstr_dirstr_[str]bst
    # stnid?_Ymd_HMS_rcustr_intstr_durstr_sst
    # stnid?_Ymd_HMS_rcustr_sbstr_intstr_durstr_dirstr_xst
    # stnid_Ymd_HMS_rcustr_sbstr_durstr_dirstr_bfs
    filefoldersplit = filefoldername.split('_')
    # Take care of possible "cal*" or "mod" tag just before ldattype:
    if (filefoldersplit[-2].startswith('cal')
            or filefoldersplit[-2].startswith('mod')
            or filefoldersplit[-2].startswith('xtr')):
        # Special modality for this ldat either
        # calibration has been applied
        # or this is a visibility model for the ldat.
        modality = filefoldersplit.pop(-2)
        obsinfo[modality] = True
    ff_suffix = filefoldersplit.pop()  # Could have extra info about ldat_type
    # Add `int<str>` field if missing
    if ldat_type == 'bfs':
        # Does not have int<int> field so insert int 0:
        filefoldersplit.insert(-2, 'int0')
    # Get `dir<str>` and `sb<str>` fields if missing
    if ldat_type != 'sst' and ldat_type != 'acc':
        # Have a sb<str> field:
        dirstr = filefoldersplit.pop()
        sbstr = filefoldersplit.pop(-3)
    else:
        # Do not have a sb<str> field:
        if ldat_type != 'acc':
            dirstr = ",,"
        else:
            dirstr = filefoldersplit.pop()
        sbstr = 'sb0:511'
    # Get `<stnid><antset>`
    if len(filefoldersplit[0]) >= 5:
        stnidantset = filefoldersplit.pop(0)
        stnid = stnidantset[:5]
        antennaset = stnidantset[5:].replace('-', '_')
        obsinfo['antennaset'] = antennaset
    else:
        stnid = None
    # Main split of filefolder into obs parameters:
    (Ymd, HMS, spwstr, intstr, durstr) = filefoldersplit[:]

    obsinfo['station_id'] = stnid
    obsinfo['filenametime'] = Ymd + '_' + HMS
    obsinfo['datetime'] = datetime.datetime.strptime(Ymd + 'T' + HMS,
                                                     '%Y%m%dT%H%M%S')
    obsinfo['spw'] = spwstr[3:]
    obsinfo['subbands'] = sbstr[2:]  # slicestr2seqlists(sbstr[2:])
    obsinfo['integration'] = float(intstr[3:])
    obsinfo['duration_scan'] = int(durstr[3:])
    obsinfo['pointing'] = dirstr[3:]
    obsinfo['ldat_type'] = ldat_type
    # Split spw str of chars into list of chars
    obsinfo['spw'] = list(obsinfo['spw'])

    # Figure out 'antennaset' field
    if not obsinfo['antennaset']:
        # Assume a sensible value for antennaset based on spw:
        _antennaset = modeparms.rcumode2antset_eu(obsinfo['spw'][0])
        obsinfo['antennaset'] = _antennaset[:3]
        if len(obsinfo['spw']) > 1:
            # If multi-spw, then set antennaset to combined antennaset
            obsinfo['antennaset'] = 'LBA+HBA'

    # Figure out 'subbands' and 'frequencies' fields
    obsinfo['subbands'] = modeparms.slicestr2seqlists(obsinfo['subbands'])
    if type(obsinfo['subbands']) is not list:
        obsinfo['subbands'] = [obsinfo['subbands']]
    obsinfo['frequencies'] = numpy.empty(0)
    beamlets = []
    totnrsbs = 0
    for spw_nr, spw in enumerate(obsinfo['spw']):
        sblist = modeparms.seqarg2list(obsinfo['subbands'][spw_nr])
        nrsbs = len(sblist)
        sblo = sblist[0]
        sbhi = sblist[-1]
        nz = modeparms.rcumode2nyquistzone(spw)
        freqlo = modeparms.sb2freq(sblo, nz)
        freqhi = modeparms.sb2freq(sbhi, nz)
        obsinfo['frequencies'] = numpy.append(obsinfo['frequencies'],
                                              numpy.linspace(freqlo, freqhi,
                                                             nrsbs))
        bmltarg = modeparms.seqlists2slicestr(
            ','.join([str(_b) for _b in range(nrsbs)]))
        beamlets.append(bmltarg)
        totnrsbs += nrsbs

    # Figure out bits and lanes
    bits = 16
    if  ldat_type == 'bst' or ldat_type == 'bfs':
        if ldat_type == 'bst':
            maxnrbls = modeparms.NRBEAMLETSBYBITS[bits]
            # Get nr of lanes and bits from filefolder suffix if it exists
            if len(ff_suffix) > 3:
                _suffix_head = ff_suffix[:-3]
                _bm, _laneshex  = _suffix_head[0], _suffix_head[1:]
                bits = modeparms.bitmode2nrbits(_bm)
                maxnrbls = modeparms.NRBEAMLETSBYBITS[bits]
                _lanesbin = bin(int(_laneshex, 16))[2:]
                lane_start = _lanesbin.index('1')
                lane_stop = list(reversed(_lanesbin)).index('1')+len(_lanesbin)
                lanes_slc = slice(lane_start, lane_stop+1)
                nrlanes = _lanesbin.count('1')
                maxnrbls = (maxnrbls // modeparms.MAX_NRLANES) * nrlanes
        if ldat_type == 'bfs':
            # When the beamlets allocated is less than the maximum (given by bit
            # depth) the RSPs fill the remaining ones regardless. Hence we have to
            # account for them:
            if totnrsbs <= modeparms.BASE_NR_BEAMLETS:
                maxnrbls = modeparms.BASE_NR_BEAMLETS
                bits = 16
            elif totnrsbs <= modeparms.BASE_NR_BEAMLETS * 2:
                maxnrbls = modeparms.BASE_NR_BEAMLETS * 2
                bits = 8
            else:
                maxnrbls = modeparms.BASE_NR_BEAMLETS * 4
                bits = 4
        missing_nr_sbs = maxnrbls - totnrsbs
        print('mis',missing_nr_sbs, totnrsbs)
        if missing_nr_sbs > 0:
            nrsbs = missing_nr_sbs
            sblo = sbhi + 1
            sbhi = sblo + nrsbs - 1
            freqlo = modeparms.sb2freq(sblo, nz)
            freqhi = modeparms.sb2freq(sbhi, nz)
            obsinfo['frequencies'] = numpy.append(obsinfo['frequencies'],
                                        numpy.linspace(freqlo, freqhi, nrsbs))
        elif missing_nr_sbs < 0:
            # Not all lanes are being used so adjust subbands & freqs
            _nrsbs = 0
            allsbs = []
            for _beam in obsinfo['subbands']:
                allsbs.extend(modeparms.seqarg2list(_beam))
            bls_pl = maxnrbls // nrlanes
            if lanes_slc.start is None:
                lanes_slc = slice(lanes_slc.stop, lanes_slc.stop+1)
            print('slc',lanes_slc)
            lane_sbs = allsbs[lanes_slc.start*bls_pl:lanes_slc.stop*bls_pl]
            print('lnsbs',len(lane_sbs))
            obsinfo['frequencies'] = obsinfo['frequencies'][range(len(lane_sbs))]

        obsinfo['max_nr_bls'] = maxnrbls
        obsinfo['nrbits'] = bits

    # Assemble _cmds
    #    rcusetup_cmds
    rcusetup_cmds = modeparms.rcusetup_args2cmds(bits, 0)
    obsinfo['rcusetup_cmds'] = rcusetup_cmds
    #    beamctl_cmds
    beamctl_cmds = []
    for spw_nr, spw in enumerate(obsinfo['spw']):
        band = modeparms.rcumode2band(spw)
        try:
            anadigdir = ilisa.operations.directions.normalizebeamctldir(
                obsinfo['pointing'])
        except ValueError:
            anadigdir = obsinfo['pointing']
        beamctl_cmd = modeparms.beamctl_args2cmds(beamlets[spw_nr],
                                                  obsinfo['subbands'][spw_nr],
                                                  band, anadigdir)
        beamctl_cmds.append(beamctl_cmd)
    obsinfo['beamctl_cmds'] = beamctl_cmds
    #     rspctl_cmds
    rspctl_cmds = modeparms.rspctl_stats_args2cmds(obsinfo['ldat_type'],
                                                   obsinfo['integration'],
                                                   obsinfo['duration_scan'],
                                                   obsinfo['subbands'])
    obsinfo['rspctl_cmds'] = rspctl_cmds
    return obsinfo


def filename2obsparm(filepath):
    """
    Convert filename to observation parameters

    Sets the dimensions of the visibility cube, which most generally is:
    dimtimes * dimrcu0 * dimrcu1.

    The file name should have the format:
        Ymd_HMS_xst_spw3_sb123_96x96.dat for xst file
        Ymd_HMS_acc_nrsbsxnrrcus0xnrrcus1.dat for acc file

    Parameters
    ----------
    filepath: str

    Returns
    -------
    obsparm: dict
        Consists of at least the keys:
            ldat_type
            datetime
            sb
        Then depending on ldat_type.
        ACC has in addition:
            cvcdim0
            cvcdim1
        and XST has in addition:
            spw
    """
    filename = os.path.basename(filepath)
    (Ymd, HMS, cvcextrest) = filename.split('_', 2)
    datatype, restdat = cvcextrest[:3], cvcextrest[3:]
    (rest, _datstr) = restdat.split('.')
    rest = rest.lstrip('_')
    if datatype == 'acc':
        (nrsamps, nrrcus0, nrrcus1) = map(int, rest.split('x'))
    elif datatype == 'xst':
        spwnr, sbnr = rest.split('_')
        spw = int(spwnr[3:])
        sb = str(int(sbnr[2:]))
    filenamedatetime = datetime.datetime.strptime(Ymd + 'T' + HMS,
                                                  '%Y%m%dT%H%M%S')
    # NOTE: For ACC, filename is last obstime, while for XST, it is first.
    if datatype == 'acc':
        filebegindatetime = filenamedatetime - datetime.timedelta(
            seconds=nrsamps)
    else:
        filebegindatetime = filenamedatetime
    obsparm = {'ldat_type': datatype,
               'datetime': filebegindatetime}
    if datatype == 'acc':
        obsparm['sb'] = modeparms.list2seqarg(list(range(nrsamps)))
        obsparm['cvcdim0'] = nrrcus0
        obsparm['cvcdim1'] = nrrcus1
    elif datatype == 'xst':
        obsparm['spw'] = spw
        obsparm['sb']= sb
    return obsparm


