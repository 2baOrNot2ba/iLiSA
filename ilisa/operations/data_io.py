"""Provides I/O of LOFAR stand-alone data.

The main data products in the stand-alone mode are files of containing
the following type of data:
   * ACC: (512, 192, 192) array of complex. (Autocovariance cube)
   * BST: (8/bits*488, T) array. (Beamlet statistics)
   * SST: (512, T) array. (Subband statistics)
   * XST: (192, 192, T) array of complex. (Cross-correlation statistics)

This module assumes that the stand-alone data files have been placed in an
appropriately named folder. The folder name contains observational settings
such as start-time, rcumode, duration, integration and pointing direction.
These folders typically contain more than one datafile representing a dataset of
the corresponding datatype.
"""
import os
import shutil
import time
import datetime

import numpy
import yaml
import argparse
import warnings

# Set ignore FutureWarning (for h5py)
warnings.filterwarnings(action='ignore', category=FutureWarning) 
import h5py
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates

try:
    import dreambeam
    canuse_dreambeam = True
except ImportError:
    canuse_dreambeam = False
if canuse_dreambeam:
    from dreambeam.polarimetry import convertxy2stokes

import ilisa
import ilisa.operations
import ilisa.operations.directions
import ilisa.antennameta.antennafieldlib as antennafieldlib
from ilisa.calim.visibilities import cov_flat2polidx
from ilisa.operations import USER_CACHE_DIR
import ilisa.operations.modeparms as modeparms


_RCU_SB_SEP = "+"


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
    elif suf == 'bst':
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


def seqlists2slicestr(seqlists):
    """
    Convert a sequence list to slice format

    Instead of comma separated list format (e.g. 202,204,206), try to construct
    subband slice syntax (e.g. 202:2:206), if possible. One use-case is in the
    construction of file names containing subband selection, in order to avoid
    file names that are potentially longer than 255 chars.

    Parameters
    ----------
    seqlists : list or str
        List of strings with comma separated numbers that are monotonically
        increasing by a constant increment.

    Returns
    -------
    slicestr : str
        Slice expression string.

    Examples
    --------
    Simple string with comma separated numbers:
    >>> from ilisa.operations.data_io import seqlists2slicestr
    >>> seqlists2slicestr('2,3,4,5,6')
    '2:6'

    Lists of strings with comma separated numbers:
    >>> seqlists2slicestr(['1,2,3','11,12,13'])
    '1:3+11:13'

    If number sequences increment by more than 1:
    >>> seqlists2slicestr(['1,3,5','12,15,18'])
    '1:2:5+12:3:18'

    """
    def seqlist2slice(seqlist):
        seqlistcanon = []
        for seqel in seqlist.split(','):
            seqel = [int(el) for el in seqel.split(':')]
            seq = range(seqel[0], seqel[-1]+1)
            seqlistcanon.extend(seq)
        seqsteps = set(numpy.diff(seqlistcanon))
        if len(seqsteps) > 1:
            raise ValueError('Subband spec {} too complicated.'.format(seqlist))
        elif len(seqsteps) == 0:
            slicestr = "{}".format(seqlistcanon[0])
        else:
            seqstep = seqsteps.pop()
            seqstepstr = str(seqstep) + ':' if seqstep > 1 else ''
            slicestr = "{}:{}{}".format(seqlistcanon[0], seqstepstr,
                                        seqlistcanon[-1])
        return slicestr

    if type(seqlists) is list:
        slicestrlist = []
        for seqlist in seqlists:
            seqstr = seqlist2slice(seqlist)
            slicestrlist.append(seqstr)
        slicestr = _RCU_SB_SEP.join(slicestrlist)
    else:
        slicestr = seqlist2slice(seqlists)
    return slicestr


def slicestr2seqlists(slicestr):
    """
    Convert slicestr to seqlists

    Inverse of seqlists2slicestr().

    Parameters
    ----------
    slicestr: str
        A slice str e.g. '1:2:5+12:3:18'.

    Returns
    -------
    seqlists: list
        A list of sequence lists e.g. ['1,3,5','12,15,18'].
    """
    seqlists = []
    slicelist = slicestr.split(_RCU_SB_SEP)
    for sl in slicelist:
        sliceargs = [int(s) for s in sl.split(':')]
        slice_inc = 1
        if len(sliceargs) == 1:
            seqstr = str(sliceargs[0])
        else:
            if len(sliceargs)>2:
                slice_inc = sliceargs[1]
                sliceargs[-1] += slice_inc
            else:
                sliceargs[-1] += 1
            seqstr = ','.join(
                [str(s) for s in range(sliceargs[0], sliceargs[-1], slice_inc)])
        seqlists.append(seqstr)
    return seqlists


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
            filefoldername += seqlists2slicestr(obsinfo['subbands'])
    if 'integration' in obsinfo and obsinfo['integration']:
        filefoldername += '_int' + str(obsinfo['integration'])
    if 'duration_scan' in obsinfo:
        filefoldername += '_dur' + str(int(obsinfo['duration_scan']))
    if ldat_type != 'sst':
        if str(obsinfo['pointing']) != '':
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
    # Format:
    # stnid_Ymd_HMS_spwstr_intstr_durstr_dirstr_[cal*]_acc
    # stnidantset_Ymd_HMS_spwstr_intstr_durstr_dirstr_[cal*]_acc
    # stnid?_Ymd_HMS_rcustr_sbstr_intstr_durstr_dirstr_bst
    # stnid?_Ymd_HMS_rcustr_intstr_durstr_sst
    # stnid?_Ymd_HMS_rcustr_sbstr_intstr_durstr_dirstr_xst
    # stnid_Ymd_HMS_rcustr_sbstr_durstr_dirstr_bfs
    filefoldersplit = filefoldername.split('_')
    # Take care of possible "cal*" or "mod" tag just before ldattype:
    if filefoldersplit[-2].startswith('cal')\
            or filefoldersplit[-2].startswith('mod'):
        # Special modality for this ldat either
        # calibration has been applied
        # or this is a visibility model for the ldat.
        modality = filefoldersplit.pop(-2)
        obsinfo[modality] = True
    ldat_type = filefoldersplit.pop()
    if ldat_type == 'bfs':
        # Does not have int<int> field so insert int 0:
        filefoldersplit.insert(-2, 'int0')
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
    if len(filefoldersplit[0]) >= 5:
        stnidantset = filefoldersplit.pop(0)
        stnid = stnidantset[:5]
        antennaset = stnidantset[5:].replace('-', '_')
        obsinfo['antennaset'] = antennaset
    else:
        stnid = None
    (Ymd, HMS, spwstr, intstr, durstr) = filefoldersplit[:]

    obsinfo['station_id'] = stnid
    obsinfo['filenametime'] = Ymd + '_' + HMS
    obsinfo['datetime'] = datetime.datetime.strptime(Ymd + 'T' + HMS,
                                                     '%Y%m%dT%H%M%S')
    obsinfo['spw'] = spwstr[3:]
    if obsinfo['antennaset'] == '':
        # Assume a sensible value for antennaset based on spw:
        obsinfo['antennaset'] = modeparms.rcumode2antset_eu(obsinfo['spw'])[:3]
    obsinfo['subbands'] = sbstr[2:]  # slicestr2seqlists(sbstr[2:])
    obsinfo['integration'] = float(intstr[3:])
    obsinfo['duration_scan'] = int(durstr[3:])
    obsinfo['pointing'] = dirstr[3:]
    obsinfo['ldat_type'] = ldat_type

    if len(obsinfo['spw']) > 1:
        obsinfo['spw'] = list(obsinfo['spw'])
    obsinfo['subbands'] = slicestr2seqlists(obsinfo['subbands'])

    if type(obsinfo['spw']) is not list:
        obsinfo['spw'] = [obsinfo['spw']]
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
        bmltarg = seqlists2slicestr(
            ','.join([str(_b) for _b in range(nrsbs)]))
        beamlets.append(bmltarg)
        totnrsbs += nrsbs
    bits = 16

    if ldat_type == 'bst' or ldat_type == 'bfs' or ldat_type == 'bstc':
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

        if missing_nr_sbs > 0:
            nrsbs = missing_nr_sbs
            sblo = sbhi + 1
            sbhi = sblo + nrsbs - 1
            freqlo = modeparms.sb2freq(sblo, nz)
            freqhi = modeparms.sb2freq(sbhi, nz)
            obsinfo['frequencies'] = numpy.append(obsinfo['frequencies'],
                                                  numpy.linspace(freqlo, freqhi,
                                                                 nrsbs))
        obsinfo['max_nr_bls'] = maxnrbls

    # Assemble _cmds
    #    rcusetup_cmds
    rcusetup_cmds = modeparms.rcusetup_args2cmds(bits, 0)
    obsinfo['rcusetup_cmds'] = rcusetup_cmds
    #    beamctl_cmds
    beamctl_cmds = []
    for spw_nr, spw in enumerate(obsinfo['spw']):
        band = modeparms.rcumode2band(spw)
        anadigdir = ','.join(obsinfo['pointing'])
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
        datatype
        filebegindatetime
        cvcdim1
        cvcdim2
    """
    filename = os.path.basename(filepath)
    (Ymd, HMS, cvcextrest) = filename.split('_', 2)
    datatype, restdat = cvcextrest[:3], cvcextrest[3:]
    (rest, _datstr) = restdat.split('.')
    rest = rest.lstrip('_')
    if datatype == 'acc':
        (_nrsamps, nrrcus0, nrrcus1) = map(int, rest.split('x'))
    elif datatype == 'xst':
        spwnr, sbnr = rest.split('_')
        spw = int(spwnr[3:])
        sb = str(int(sbnr[2:]))
    filenamedatetime = datetime.datetime.strptime(Ymd + 'T' + HMS,
                                                  '%Y%m%dT%H%M%S')
    # NOTE: For ACC, filename is last obstime, while for XST, it is first.
    if datatype == 'acc':
        filebegindatetime = filenamedatetime - datetime.timedelta(
            seconds=_nrsamps)
    else:
        filebegindatetime = filenamedatetime
    obsparm = {'ldat_type': datatype,
               'datetime': filebegindatetime,
               'spw': spw,
               'sb': sb}
    return obsparm


def dataff_raw_model_cal(dataff):
    """
    Return nominal names of raw, model and calibrated data file-folders.

    Parameters
    dataff: str
        Path to a data file-folder.

    Returns
    -------
    dataff_raw: str
        Path to raw data file-folder.
    dataff_mod: str
        Path to model data file-folder.
    dataff_cal: str
        Path to calibrated data file-folder.
    """
    dataff = os.path.normpath(dataff)
    dataff_dir, dataff = os.path.split(dataff)
    lofar_datatype = datafolder_type(dataff)
    if lofar_datatype != 'acc' and lofar_datatype != 'xst':
        raise RuntimeError("Datafolder '{}'\n not ACC or XST type data."
                           .format(dataff))
    obsinfo_raw = filefolder2obsinfo(dataff)

    # Determine 'raw', 'mod', and 'cal' dataffs,
    obsinfo_raw.pop('cal', None)
    obsinfo_raw.pop('model', None)
    dataff_raw = obsinfo2filefolder(obsinfo_raw)
    dataff_raw = os.path.join(dataff_dir, dataff_raw)

    obsinfo_mod = dict(obsinfo_raw)
    obsinfo_mod['model'] = True
    dataff_mod = obsinfo2filefolder(obsinfo_mod)
    dataff_mod = os.path.join(dataff_dir, dataff_mod)

    obsinfo_cal = dict(obsinfo_raw)
    obsinfo_cal['cal'] = True
    dataff_cal = obsinfo2filefolder(obsinfo_cal)
    dataff_cal = os.path.join(dataff_dir, dataff_cal)

    return dataff_raw, dataff_mod, dataff_cal


class ScanRecInfo(object):
    """This class maintains info on a scan recording (scanrec), which is one
    of the results of an iLiSA scan. One scanrec is a group of one or more
    files of a unique LOFAR station data product (ldat), i.e. acc, bfs, bst,
    sst or xst. The info in a ScanRecInfo object consists of sufficient
    parameters to redo the data, namely the stn_id, the iLiSA scanrec
    parameters, and a list of LDatInfo objects called ldatinfos that maps to
    each ldat within the scanrec.
    """
    scanrecinfo_header = "SCANREC_INFO.yml"

    def __init__(self, stnid='', caltabinfos=[], mockdata='', sourcename=''):
        self.headerversion = 5
        self.ldatinfos = {}
        self._pointing = ''
        self.sourcename = sourcename
        self.caltabinfos = caltabinfos
        self.scanrecpath = None
        self.scanrecparms = {}
        self.obs_ids = []
        self.stnid = stnid
        # Data modalities (if none are set, then normal observed data):
        self.calibrationfile = ''
        self.gs_model = ''
        self.mockdata = mockdata
        self.comments = ''

    def add_obs(self, ldatinfo):
        """Add an LDatInfo object to this ScanRecInfo."""
        obs_id = ldatinfo.filenametime
        self.ldatinfos[obs_id] = ldatinfo

    def get_obs_ids(self):
        """
        Get the list of obs_ids

        An obs_id is a key to the ldatinfos list
        """
        _obs_ids = self.ldatinfos.keys()
        return sorted(_obs_ids)

    def list_ldatfiles(self):
        # Select only data files in folder (avoid CalTable*.dat files)
        ls = os.listdir(self.scanrecpath)
        filenames = [filename for filename in ls if filename.endswith('.dat')
                     and not filename.startswith('CalTable')]
        filenames.sort()  # This enforces chronological order
        return filenames

    def set_stnid(self, stnid):
        self.stnid = stnid

    def get_stnid(self):
        try:
            stnid = self.stnid
        except:
            try:
                stnid = self.ldatinfos[0].stnid
            except:
                try:
                    stnid = self.scanrecparms['station']
                except:
                    raise RuntimeError('Station id not found.')
        return stnid

    def set_scanrecparms(self, ldat_type, freqspec, duration,
                         direction="None,None,None", integration=1.0,
                         antset=''):
        """
        Set ScanRecInfo's dict attribute `scanrecparms` parameters

        Parameters
        ----------
        ldat_type: str
            Ldat type.
        freqspec: str
            Argument to FreqSetup()
        duration: float
            Duration in seconds.
        direction: str
            Comma separated direction triplet str.
        integration: int
            Integration in seconds.
        antset: str
            Antenna-set name.
        """
        self.scanrecparms['datatype'] = ldat_type
        self.scanrecparms['freqspec'] = freqspec
        self.scanrecparms['duration'] = duration
        self.scanrecparms['direction'] = direction
        self.scanrecparms['integration'] = integration
        if antset:
            self.scanrecparms['antennaset'] = antset

    def set_caltabinfos(self, caltabinfos):
        self.caltabinfos = caltabinfos

    def write_scanrec(self, datapath):
        with open(os.path.join(datapath, self.scanrecinfo_header), "w") as f:
            f.write("# Scan recording header\n")
            f.write("# Created by {} version {}\n".format("iLiSA",
                                                          ilisa.__version__))
            f.write("headerversion: {}\n".format(self.headerversion))
            f.write("station: {}\n".format(self.stnid))
            f.write("scanrecparms: {!r}\n".format(self.scanrecparms)
                    .replace('None', 'null'))
            if self.sourcename:
                f.write("sourcename: {}\n".format(self.sourcename))
            f.write("ldat_ids: {!r}\n".format(list(self.ldatinfos.keys())))
            # Data modalities:
            if self.caltabinfos != []:
                f.write("caltabinfos: {}\n".format(self.caltabinfos))
            if self.gs_model:
                f.write("gs_model: {}\n".format(self.gs_model))
            if self.mockdata:
                f.write("mockdata: true\n")

    def read_scanrec(self, datapath):
        """
        Read in scan-rec parameters from nominal file in datapath

        Parameters
        ----------
        datapath: str
            Path to data file-folder.
        """
        try:
            _h_path = os.path.join(datapath, self.scanrecinfo_header)
            with open(_h_path, 'r') as hf:
                scanrecfiledict = yaml.safe_load(hf)
        except Exception:
            warnings.warn(
                "Couldn't find Scanrecinfo file. Will try filefolder name...")
            self.read_scanrec_from_ff(datapath)
            return
        self.headerversion = scanrecfiledict['headerversion']
        self.stnid = scanrecfiledict['station']
        self.scanrecparms = scanrecfiledict['scanrecparms']
        self.sourcename = scanrecfiledict.get('sourcename')
        self.obs_ids = scanrecfiledict.get('ldat_ids', None)
        # Data modalities:
        self.calibrationfile = scanrecfiledict.get('calibrationfile', '')
        self.gs_model = scanrecfiledict.get('gs_model', '')
        self.mockdata = scanrecfiledict.get('mockdata', False)
        self.comments = scanrecfiledict.get('comments', '')
        return self

    def read_scanrec_from_ff(self, datapath):
        """
        Read in scan-rec parameters from file-folder name

        Parameters
        ----------
        datapath: str
             Path to data file-folder.
        """
        try:
            obsinfo = filefolder2obsinfo(datapath)
        except ValueError:
            print("Could not parse filefolder {}".format(datapath))
            # Will hope to read LDat header
            self.scanrecparms = None
        else:
            spw = obsinfo['spw'][0]
            nqz = modeparms.rcumode2nyquistzone(spw)
            sbs = modeparms.seqarg2list(obsinfo['subbands'][0])
            freqspec_hi = modeparms.sb2freq(sbs[-1], nqz)
            self.set_scanrecparms(obsinfo['ldat_type'], str(freqspec_hi),
                                  obsinfo['duration_scan'], obsinfo['pointing'],
                                  obsinfo['integration'], obsinfo['antennaset'])
            self.scanrecparms['rcumode'] = spw
            self.set_stnid(obsinfo['station_id'])
            self.calibrationfile = None
            print("Read in filefolder meta.")

    def set_postcalibration(self, caltabpath, scanrecpath):
        """Add the caltab file that was applied to this scanrec (typically CVC
        data) after station recording.

        Parameters
        ----------
        caltabpath: str
            Path to caltab file
        scanrecpath: str
            Path to scanrec folder
        """
        shutil.copy(caltabpath, scanrecpath)
        scanrecinfo_header_path = os.path.join(scanrecpath,
                                               self.scanrecinfo_header)
        with open(scanrecinfo_header_path, 'a') as h:
            h.write("calibrationfile: " + os.path.basename(caltabpath) + '\n')

    def set_model(self, gs_model):
        """
        Set the global-sky model used to generate this data

        Parameters
        ----------
        gs_model: str
            Global-sky model name.
        """
        scanrecinfo_header_path = os.path.join(self.scanrecpath,
                                               self.scanrecinfo_header)
        with open(scanrecinfo_header_path, 'a') as h:
            h.write("gs_model: " + gs_model + '\n')

    def get_datatype(self):
        return self.scanrecparms['datatype']

    def get_rcumode(self, filenr=0):
        try:
            rcumode = modeparms.FreqSetup(self.scanrecparms['freqspec'])._rcumodes[0]
        except:
            try:
                rcumode = self.ldatinfos[filenr].beamctl_cmd['rcumode']
            except:
                rcumode = self.scanrecparms['rcumode']
        return str(rcumode)

    def get_band(self):
        return modeparms.rcumode2band(self.get_rcumode())

    def get_bandarr(self):
        antset = modeparms.rcumode2antset_eu(self.get_rcumode())
        return antset.split('_')[0]

    def get_antset(self):
        antset  = self.scanrecparms.get('antennaset')
        if not antset:
            antset = self.ldatinfos[self.get_obs_ids()[0]].antset
        return antset

    def get_xcsubband(self, filenr=0):
        return int(self.ldatinfos[filenr].rspctl_cmd['xcsubband'])

    def get_integration(self):
        return self.scanrecparms['integration']

    def get_pointingstr(self, filenr=0):
        return self.scanrecparms['direction']

    def is_septon(self, filenr=0):
        obs_ids = self.get_obs_ids()
        try:
            self.ldatinfos[obs_ids[filenr]]
        except:
            if self.get_datatype().endswith('SEPTON'):
                return True
            else:
                return False
        else:
            if self.ldatinfos[obs_ids[filenr]].septonconf:
                return True
            else:
                return False

    def get_septon_elmap(self, filenr=0):
        obs_ids = self.get_obs_ids()
        elmap = modeparms.str2elementMap2(
                    self.ldatinfos[obs_ids[filenr]].septonconf)
        return elmap

    def get_ldat_filenames(self):
        datatype = self.get_datatype()
        if datatype == 'acc':
            nameformat = "{}_{}_512x192x192.dat"
        else:
            nameformat = "{}_{}.dat"
        return [nameformat.format(obs_id, datatype) for obs_id in self.obs_ids]
    
    def get_allsky(self):
        """Determine if allsky FoV"""
        band = self.get_band()
        septon = self.is_septon()
        if band == '10_90' or band == '30_90' or septon:
            allsky = True
        else:
            allsky = False
        return allsky


class LDatInfo(object):
    """Maintains a minimum of info for the physical interpretation of an LCU
    data product file (ldat). It consists of the LCU commands or settings that
    resulted in the ldat. An observation can update an object with LCU settings
    and the create a header file associated with the resulting ldat file of the
    observation. Later when ingesting the ldat file, an object can read-in the
    associated header file to get physically meaningful data.

    Specifically an object should have:

    filename:

    """

    def __init__(self, ldat_type, rcusetup_cmds, beamctl_cmds, rspctl_cmds,
                 septonconf=None, **kwargs):
        """Create observation info from parameters."""
        self.headerversion = '4'

        # ldat_type attr
        self.ldat_type = ldat_type
        # filenametime attr (Is set later, since it's know only after obs)
        self.filenametime = None

        # rcusetup_cmds attr
        if rcusetup_cmds == []:
            rcusetup_cmds = ['rspctl']
        self.rcusetup_cmds = rcusetup_cmds
        rcusetup_args = modeparms.parse_rspctl_args(self.rcusetup_cmds)
        self.attenuation = None
        if 'attentuation' in rcusetup_args:
            self.attenuation = rcusetup_args['attentuation']
        self.bits = int(rcusetup_args.get('bitmode', 16))  # 16 is default
        self.mode = rcusetup_args.get('mode', None)

        # beamctl_cmds related attr
        self.rcumode = []
        self.sb = []
        self.bl = []
        self.direction = None
        self.beamctl_cmds = beamctl_cmds
        if self.beamctl_cmds != []:
            digdir = None
            if type(self.beamctl_cmds) is not list:
                self.beamctl_cmds = [self.beamctl_cmds]
            for _beamctl_cmd in self.beamctl_cmds:
                (antset, _rcus, rcumode, beamlets, subbands, _anadir,
                    digdir) = modeparms.parse_beamctl_args(_beamctl_cmd)
                self.antset = antset
                self.rcumode.append(int(rcumode))
                self.sb.append(subbands)
                self.bl.append(beamlets)
            self.direction = digdir
        else:
            # No beamctl was issued so figure out LCU setup for other info
            if self.mode:
                self.rcumode = self.mode
                self.antset = modeparms.rcumode2antset_eu(self.mode)

        # rspctl_cmds attr
        if not rspctl_cmds:
            rspctl_cmds = ['rspctl']
        self.rspctl_cmds = rspctl_cmds
        rspctl_args = modeparms.parse_rspctl_args(self.rspctl_cmds)
        
        # septonconf attr
        self.septonconf = septonconf
        if self.septonconf is not None:
            self.rcumode = [5]
            self.antset = 'HBA'
        
        # attrs: integration, duration_scan
        if self.ldat_type != 'bfs':
            if self.ldat_type != 'acc':
                self.integration = float(rspctl_args['integration'])
                self.duration_subscan = float(rspctl_args['duration'])
            else:
                self.integration = 1.0
                self.duration_subscan = 512
        # attrs: sb
        if self.ldat_type == 'sst' or self.ldat_type == 'acc':
            self.sb = ""
        elif self.ldat_type.startswith('xst'):
            self.sb = str(rspctl_args['xcsubband'])
            # self.rcumode = self.rcumode[0]
        elif self.ldat_type == 'bst':
            self.sb = self.sb


    def write_ldat_header(self, datapath):
        """Create a header file for LOFAR standalone observation."""
        contents = {}
        contents['ldat_type'] = self.ldat_type
        contents['filenametime'] = self.filenametime
        contents['rcusetup_cmds'] = self.rcusetup_cmds
        contents['beamctl_cmds'] = self.beamctl_cmds
        contents['rspctl_cmds'] = self.rspctl_cmds
        if self.septonconf:
            contents['septonconf'] = self.septonconf

        if not modeparms.is_ldattype(self.ldat_type):
            raise ValueError("Unknown LOFAR statistic type {}."
                             .format(self.ldat_type))
        xtra = ''
        if self.ldat_type == 'acc':
            xtra = '_512x192x192'
        ldat_header_filename = (self.filenametime + '_' + self.ldat_type
                                + xtra + '.h')
        with open(os.path.join(datapath, ldat_header_filename), 'w') as f:
            f.write('# LCU obs settings, header file\n')
            f.write('# Header version'+' '+self.headerversion+'\n')
            yaml.dump(contents, f, default_flow_style=False, width=1000)

    def get_recfreq(self, sampnr=0):
        """Return data recording frequency in Hz

        Parameters
        ----------
        sampnr : int
            Sample number from start of file

        Returns
        -------
        freq : float
            Frequency for this covariance sample
        """
        if self.ldat_type == 'acc':
            sb = sampnr
        else:
            sb = self.sb
        if self.ldat_type != "xst-SEPTON" and not self.septonconf:
            rcumode = self.rcumode[0]
        else:
            rcumode = 5
        nz = modeparms.rcumode2nyquistzone(rcumode)
        return modeparms.sb2freq(sb, nz)

    def get_starttime(self):
        """Return the datetime when this obs started."""
        filetime = datetime.datetime.strptime(self.filenametime,
                                              "%Y%m%d_%H%M%S")
        if self.ldat_type != 'acc':
            starttime = filetime
        else:
            starttime = filetime - datetime.timedelta(seconds=512)
        return starttime

    def get_spw(self):
        """\
        Get Spectral Window

        Could either be rcumode (as set with beamctl) or mode (set with rspctl).
        """
        if self.rcumode:
            return self.rcumode
        elif self.mode:
            return self.mode
        else:
            return None

    @classmethod
    def from_filename(cls, datfilename, anadigdir=',,', integration=1.0,
                      duration_scan=1.0):
        """
        Create an LDatInfo from an ldat filename

        Parameters
        ----------
        datfilename: str
            Formatted name of ldat file.
        anadigdir: str
            Analog-digital direction str.
        integration: float
            Integration time.
        duration_scan: float
            Duration of scan.

        Returns
        -------
        ldatinfo: LDatInfo
            Corresponding LDatInfo object.
        """
        obsparm = filename2obsparm(datfilename)
        ldat_type = obsparm['ldat_type']
        sb = obsparm['sb']
        rcusetup_cmds = ''
        anadigdir = ilisa.operations.directions.normalizebeamctldir(anadigdir)
        beamlets = modeparms.alloc_beamlets(sb)[0]
        band = modeparms.rcumode2band(obsparm['spw'])
        beamctl_cmds = modeparms.beamctl_args2cmds(beamlets=beamlets,
                                                   subbands=sb,
                                                   band=band,
                                                   anadigdir=anadigdir)
        rspctl_cmds = modeparms.rspctl_stats_args2cmds(ldat_type,
                                                       integration,
                                                       duration_scan,
                                                       sb)
        ldatinfo = cls(ldat_type, rcusetup_cmds, beamctl_cmds, rspctl_cmds)
        ldatinfo.filenametime = obsparm['datetime'].strftime('%Y%m%d_%H%M%S')
        return ldatinfo

    @staticmethod
    def headerfromdatfile(datfile):
        """
        Determine header name from dat file name

        Pareamters
        ----------
        datfile: str
            Name of ldat file.

        Returns
        -------
        hfilename: str
            Name of header file.
        """
        # Try to get obsfile header
        (bfilename, _dat) = datfile.split('.')
        hfilename = bfilename + '.h'
        # Check if xst might have some extra stuff in name
        #   So first get ldattype
        ymd, hms, ldattype_full = bfilename.split('_', 2)
        if '_' in ldattype_full:
            ldattype, _rest = ldattype_full.split('_', 1)
        else:
            ldattype = ldattype_full
        if ldattype == 'xst':
            hfilename = ymd + '_' + hms + '_' + ldattype + '.h'
        return hfilename

    @classmethod
    def read_ldat_header(cls, headerpath):
        """Parse a ldat header file and return it as an LDatInfo."""
        # TODO extract CalTable info.
        if os.path.isdir(headerpath):
            files = os.listdir(headerpath)
            headerfiles = [f for f in files if f.endswith('.h')]
            headerfile = os.path.join(headerpath, headerfiles.pop())
        else:
            headerfile = headerpath
        stnid = None
        starttime = None
        headerversion = 0
        with open(headerfile, 'r') as hf:
            for hline in hf:
                if "Header version" in hline:
                    headerversion = hline.split()[-1]
        beamctl_line = ""
        contents = {}
        datatype = None
        with open(headerfile, 'r') as hf:
            if headerversion == '1':
                rspctl_lines = []
                for line in hf:
                    if "Observer" in line:
                        _label, _observer = line.split('=')
                    if "Project" in line:
                        _label, _project = line.split('=')
                    if "DataType" in line:
                        _label, datatype = line.split('=')
                    if "StationID" in line:
                        _label, stnid = line.split('=')
                        stnid = stnid.strip()
                    if "StartTime" in line:
                        _label, starttime = line.split('=')
                        starttime = starttime.strip()
                    if "beamctl" in line:
                        # HACK
                        beamctl_line = line
                    if "rspctl" in line:
                        rspctl_lines.append(line)
            elif headerversion == '2':
                contents = yaml.safe_load(hf)
                _observer = contents['Observer']
                _project = contents['Project']
                datatype = contents['DataType']
                stnid = contents['StationID']
                starttime = contents['StartTime']
                beamctl_line = contents['BeamctlCmds']
                rspctl_lines = contents['RspctlCmds'].split('\n')
            elif headerversion == '3':
                contents = yaml.safe_load(hf)
                datatype = contents['datatype']
                filenametime = contents['filetime']
                rcusetup_cmds = []
                beamctl_cmds = [contents['beamctl_cmds']]
                rspctl_cmds = [contents['rspctl_cmds'].strip(' ;')]
                if 'caltabinfos' in contents:
                    caltabinfos = contents['caltabinfos']
                else:
                    caltabinfos = []
                if 'septonconf' in contents:
                    septonconf = contents['septonconf']
                else:
                    septonconf = None
            else:
                # headerversion == '4':
                contents = yaml.safe_load(hf)
                datatype = contents['ldat_type']
                filenametime = contents['filenametime']
                # stnid = contents['station_id']
                rcusetup_cmds = contents['rcusetup_cmds']
                beamctl_cmds = contents['beamctl_cmds']
                rspctl_cmds = contents['rspctl_cmds']
                if 'caltabinfos' in contents:
                    caltabinfos = contents['caltabinfos']
                else:
                    caltabinfos = []
                if 'septonconf' in contents:
                    septonconf = contents['septonconf']
                else:
                    septonconf = None
        ldatinfo = cls(datatype, rcusetup_cmds, beamctl_cmds, rspctl_cmds,
                       caltabinfos=caltabinfos, septonconf=septonconf)
        ldatinfo.filenametime = filenametime
        return ldatinfo


def readbstfolder(bst_filefolder):
    """\
    Read a BST file-folder and return data

    Parameters
    ----------
    bst_filefolder : str
        Path to BST file-folder.

    Returns
    -------
    bst_data_x : array
        Data array from X antenna of power per beamlet over time.
    bst_data_y : array
        Data array from Y antenna of power per beamlet over time.
    bst_data_xy : array
        Complex cross-correlation of antenna's X & Y per beamlet over time.
    ts : array
        Sample epoch times.
    freqs : array
        Frequencies.
    obsinfo :
        Observation metadata.
    """
    obsinfo = filefolder2obsinfo(bst_filefolder)
    maxnrsbs = obsinfo['max_nr_bls']
    intg = obsinfo['integration']
    freqs = obsinfo['frequencies']
    bst_dirls = os.listdir(bst_filefolder)
    bst_files = sorted([f for f in bst_dirls if f.endswith('.dat')])

    # Now read the BST pol data
    bst_dtype = numpy.dtype(('f8', (maxnrsbs,)))
    bstc_dtype = numpy.dtype(('c16', (maxnrsbs,)))
    bst_data_xx, bst_data_yy, bst_data_xy = [], [], []
    ts = []

    for bst_polfile in bst_files:
        pol = bst_polfile.rstrip('.dat').split('_')[-1].lstrip('0')
        with open(os.path.join(bst_filefolder, bst_polfile), 'rb') as fin:
            if pol != 'XY':
                _bst_dtype = bst_dtype
            else:
                _bst_dtype = bstc_dtype
            filedata = numpy.fromfile(fin, dtype=_bst_dtype)
        if pol == 'X' or pol == 'XX':
            bst_data_xx.append(filedata)
            Ymd, HMS, _ = bst_polfile.split('_', 2)
            file_start_dt = datetime.datetime.strptime(Ymd + '_' + HMS,
                                                       '%Y%m%d_%H%M%S')
            file_dur = filedata.shape[0]*intg
            ts_rel = numpy.arange(0., file_dur, intg)
            file_ts =\
                [file_start_dt + datetime.timedelta(microseconds=t*1000000)
                 for t in ts_rel]
            ts.append(file_ts)
        elif pol == 'Y' or pol == 'YY':
            bst_data_yy.append(filedata)
            # Assumes 'Y' identical to 'X'
        elif pol == 'XY' or pol == 'YX':
            bst_data_xy.append(filedata)

    return bst_data_xx, bst_data_yy, bst_data_xy, ts, freqs, obsinfo


def readsstfolder(sstfolder):
    """Read-in SST datafile.

    Parameters
    ----------
    sstfolder : str
        The name of the folder which contains an SST datafile for each RCU.

    Returns
    -------
    sstdata_rcu : (192, N, S, 512)
        The SST data, where N is number of files per RCU
        and S is the number of time samples with a file.
    ts : (N, S)
        The datetimes over files and samples.
    freqs : (F)
        Array of the subband freqs.
    obsinfo : dict
        Observation metadata
    """
    obsinfo = filefolder2obsinfo(sstfolder)
    intg = obsinfo['integration']
    freqs = obsinfo['frequencies']
    files = os.listdir(sstfolder)
    sstfiles = [f for f in files if f.endswith('.dat')]
    sstfiles.sort()
    ts = []
    # For each RCU initialize a list for each (rcu) file SST data
    sstdata_rcu = [[] for _ in range(192)]
    for sstfile in sstfiles:
        # Read sst file
        sst_filepath = os.path.join(sstfolder, sstfile)
        try:
            (Ymd, HMS, _sststr, rcudatstr) = sstfile.split('_')
            file_start_dattim = datetime.datetime.strptime(Ymd + 'T' + HMS,
                                                           '%Y%m%dT%H%M%S')
            (rcu, _datext) = rcudatstr[3:].split('.')
            rcu = int(rcu)
        except:
            raise ValueError("File name {} not in sst format.".format(sstfile))
        # Now read the SST data
        sst_dtype = numpy.dtype(('f8', (512,)))

        # Normally just read all data into memory...
        sstfiledata = numpy.fromfile(sst_filepath, dtype=sst_dtype)
        # ...but in case data does not fit into memory (24h @ 1s) use memmap:
        # #sstfiledata = numpy.memmap(sst_filepath, dtype=sst_dtype, mode='r')
        # Note that using memmap can cause a linux error for too many open files
        # in which case use bash> ulimit -n <N>  where <N> is no. files

        # Append file data array for this rcu to full data array
        sstdata_rcu[rcu].append(sstfiledata)
        del sstfiledata
        # Assume time of samples is same for all RCU files;
        # deal only with 1st one
        if rcu == 0:
            Ymd, HMS, _ = sstfile.split('_', 2)
            file_nrsmps = sstdata_rcu[0][0].shape[0]
            file_dur = intg * file_nrsmps
            ts_rel = numpy.arange(0., file_dur, intg)
            file_ts = [file_start_dattim + datetime.timedelta(seconds=t)
                       for t in ts_rel]
            ts.append(file_ts)
    return sstdata_rcu, ts, freqs, obsinfo


class CVCfiles(object):
    """Provides functionality for covariance cube (CVC) files

    CVC is essentially visiblity cubes. CVC files from a LOFAR station includes
    ACC and XST files. The dataset files are accessed on demand through a
    __getitem__ call. Each item corresponds to one CVC file and is the actual
    covariance matrix cube with shape cvcdim0 x cvcdim1 x cvcdim2.

    Attributes
    ----------
    samptimeset: list of datetimes
        The datetime of the visibility matrix sample.
    freqset: list of floats
        The frequencies of the subbands used.
    """
    NRRCUS_EU = 192  # Default number of RCUs on EU stations

    def __init__(self, datapath):
        self._getitem_prev_file_nr_data_ = (None, None)
        self.samptimeset = []
        self.freqset = []

        self.cvcdim1 = self.NRRCUS_EU
        self.cvcdim2 = self.NRRCUS_EU

        datapath = os.path.abspath(datapath)
        if not os.path.isdir(datapath):
            if os.path.isfile(datapath):
                # FIXME:
                self._readcvcfile(datapath)
            else:
                raise ValueError('Path {} does not exist'.format(datapath))
        self.filefolder = datapath

        # Now scan CVC filefolder to determine datafile contents
        #
        # The filefolder name may have the format as specified in the
        # parse_cvcfolder() method. The contents of the data file is stored in
        # the object attribute:
        #       data : [(N,192,192), ... , (N,192,192)]
        # where N is nominally the number of time samples and the len of data is
        # the number of files in the folder.

        # Initialize
        scanrecinfo = ScanRecInfo()
        samptimeset = []
        freqset = []
        scanrecinfo.scanrecpath = self.filefolder
        scanrecinfo.read_scanrec(self.filefolder)
        stnid = scanrecinfo.get_stnid()
        nrrcus = modeparms.nrrcus_stnid(stnid)
        self.cvcdim1 = nrrcus
        self.cvcdim2 = nrrcus
        filenames = scanrecinfo.list_ldatfiles()
        # Scan filefolder for ldat file and process
        for cvcfile in filenames:
            cvcdim_t = (os.path.getsize(os.path.join(self.filefolder, cvcfile))
                        // self.__get_cvc_dtype().itemsize)
            hfilename = LDatInfo.headerfromdatfile(cvcfile)
            hfilepath = os.path.join(self.filefolder, hfilename)
            try:
                ldatinfo = LDatInfo.read_ldat_header(hfilepath)
            except:
                warnings.warn(
                    "Couldn't find a header file for {}".format(cvcfile))
                ldatinfo = LDatInfo.from_filename(cvcfile,
                                                  duration_scan=
                                                  scanrecinfo.scanrecparms[
                                                      'duration'])
            scanrecinfo.add_obs(ldatinfo)
            _datatype, t_begin = self._parse_cvcfile(
                os.path.join(self.filefolder, cvcfile))

            # Compute time of each autocovariance matrix sample per subband
            integration = scanrecinfo.get_integration()
            obscvm_datetimes = [None] * cvcdim_t
            for t_idx in range(cvcdim_t):
                t_delta = datetime.timedelta(
                    seconds=t_idx * integration
                )
                obscvm_datetimes[t_idx] = t_begin + t_delta
            samptimeset.append(obscvm_datetimes)

            # Compute frequency of corresponding time sample
            rcumode = scanrecinfo.get_rcumode()
            nz = modeparms.rcumode2nyquistzone(rcumode)
            if scanrecinfo.get_datatype() == 'acc':
                freqs = modeparms.rcumode2sbfreqs(rcumode)
            else:
                sb = ldatinfo.sb
                freq = modeparms.sb2freq(sb, nz)
                freqs = [freq] * cvcdim_t
            freqset.append(freqs)
        (self.scanrecinfo, self.filenames, self.samptimeset, self.freqset
         ) = scanrecinfo, filenames, samptimeset, freqset
        # Get/Compute ant positions
        antset = self.scanrecinfo.get_antset()
        self.stn_pos, self.stn_rot, self.stn_antpos, self.stn_intilepos \
            = antennafieldlib.get_antset_params(stnid, antset)

        # Account for SEPTON antenna positions
        septon = self.scanrecinfo.is_septon()
        if septon:
            elmap = self.scanrecinfo.get_septon_elmap()
            for tile, elem in enumerate(elmap):
                self.stn_antpos[tile] += self.stn_intilepos[elem]

    def __getitem__(self, filenr):
        """Get data from a CVC file (in this set of files) by filenr

        (Caches file data)
        """
        prev_filenr, prev_datafromfile = self._getitem_prev_file_nr_data_
        if filenr == prev_filenr:
            datafromfile = prev_datafromfile
        else:
            cvcfile = self.filenames[filenr]
            cvcpath = os.path.join(self.filefolder, cvcfile)
            datafromfile, _t_begin = self._readcvcfile(cvcpath)
            self._getitem_prev_file_nr_data_ = (filenr, datafromfile)
        return datafromfile

    def __setitem__(self, filenr, data_arr):
        """Save CVC data item filenr data into an ldat file."""
        cvcfile = self.filenames[filenr]
        cvcpath = os.path.join(self.filefolder, cvcfile)
        data_arr.tofile(cvcpath)

    def __get_cvc_dtype(self, cvcdim1=None, cvcdim2=None):
        if cvcdim1 is None:
            cvcdim1 = self.cvcdim1
        if cvcdim2 is None:
            cvcdim2 = self.cvcdim2
        cvc_dtype = numpy.dtype(('c16', (cvcdim1, cvcdim2)))
        return cvc_dtype

    def _parse_cvcfile(self, cvcfilepath):
        """Parse the cvc file name.

        Sets the dimensions of the visibility cube, which most generally is:
        dimtimes * dimrcu0 * dimrcu1.

        The file name should have the format:
            Ymd_HMS_xst_spw3_sb123_96x96.dat for xst file
            Ymd_HMS_acc_nrsbsxnrrcus0xnrrcus1.dat for acc file

        :param cvcfilepath: str
        :return: datatype, filebegindatetime, cvcdim1, cvcdim2
        """
        cvcfilename = os.path.basename(cvcfilepath)
        (Ymd, HMS, cvcextrest) = cvcfilename.split('_', 2)
        datatype, restdat = cvcextrest[0:3], cvcextrest[3:]
        (rest, _datstr) = restdat.split('.')
        _nr512 = 512
        if datatype == 'acc':
            rest = rest.lstrip('_')
            (_nr512, nrrcus0, nrrcus1) = map(int, rest.split('x'))
        filenamedatetime = datetime.datetime.strptime(Ymd + 'T' + HMS,
                                                      '%Y%m%dT%H%M%S')
        # NOTE: For ACC, filename is last obstime, while for XST, it is first.
        if datatype == 'acc':
            filebegindatetime = filenamedatetime - datetime.timedelta(
                seconds=_nr512)
        else:
            filebegindatetime = filenamedatetime
        return datatype, filebegindatetime

    def _readcvcfile(self, cvcfilepath):
        """Reads in a single acc or xst data file by filepath and creates
        corresponding sample times.

        The contents of the data file is appended to the object attribute list
        `data`.

        Parameters
        ----------
        cvcfilepath : str
            Full path to cvc file

        Returns
        -------
        datafromfile, t_begin
            Data contents of file. Time stamp.
        """
        _datatype, t_begin = self._parse_cvcfile(cvcfilepath)
        # Get cvc data from file.
        cvc_dtype = self.__get_cvc_dtype()
        print("Reading cvcfile: {}".format(cvcfilepath))
        with open(cvcfilepath, 'rb') as fin:
            datafromfile = numpy.fromfile(fin, dtype=cvc_dtype)
        return datafromfile, t_begin

    def getnrfiles(self):
        """Return number of data files in this filefolder."""
        return len(self.filenames)

    def getfreqs(self):
        """Return tuple of the frequency axis of the data"""
        return tuple(sorted(set([freq for f in self.freqset for freq in f])))

    def getfsamps(self):
        """Return number of data samples per frequency"""
        fsamps = {}
        freqset_flat = [freq for fil in self.freqset for freq in fil]
        freqs = sorted(set(freqset_flat))  # freqs = self.getfreqs()
        for freq in freqs:
            fsamps[freq] = freqset_flat.count(freq)
        return fsamps

    def get_positions_ITRF(self):
        """Return antenna array positions in ITRF coordinates

        Returns
        -------
        pos_itrf : array
            (nrants, 3) array of absolute 3D cartesian element positions in
            meters with respect to ITRF. It is computed by adding the station
            ITRF position to the relative ITRF positions of the elements.
        """
        pos_itrf = self.stn_pos.squeeze() + self.stn_antpos
        return pos_itrf

    def as_array(self):
        """Return data as a numpy array

        Returns
        -------
        data : array
            Complex numpy array with indices:
            ['filenr', 'sampnr', 'antnr', 'antnr']
        """
        data = []
        for fnr in range(self.getnrfiles()):
            data.append(self[fnr])
        data = numpy.asarray(data)
        return data


def readacc2bst(anacc2bstfilepath, datformat='hdf'):
    """\
    Read an acc2bst file.

    The fileformat can be either hdf or numpy.

    Parameters
    ----------
    anacc2bstfilepath: str
        Full path to an acc2bst file.
    datformat: str, optional
        Dataformat of file to be read. Can be be 'hdf' or 'npy' (default: 'hdf')

    Returns
    -------
    bst_pols: tuple of floatarray
        Tuple contains polarized BST-like data: (bst_xx, bst_xy, bst_yy).
    filestarttimes: list of datetime
        List of start times of the file partitioned data.
    freqs: list of floats
        List of frequencies.
    calrunstarttime: datetime
        Start time of this calibration run.
    calrunduration: float
        Duration of this calibration run.
    calsrc: str
        Source pointed to and which is intended as a calibration source.
    caltab_id: str
        Calibration table ID used for calibrating the input visibilities.
    stnid: str
        ID string of the station used.
    used_autocorr: bool
        Was autocorrelation used in beamforming or not?

    See Also
    --------
    saveacc2bst: converse of this function.
    """
    anacc2bstfilepath = os.path.abspath(anacc2bstfilepath)
    acc2bstfiledir = os.path.dirname(anacc2bstfilepath)
    anacc2bstfilename = os.path.basename(anacc2bstfilepath)
    (stnid, begin_utc_str, rcuarg, calsrc, durarg, caltabdate, acc2bst, _version
     ) = anacc2bstfilename.split('_')
    calrunstarttime = datetime.datetime.strptime(begin_utc_str, "%Y%m%dT%H%M%S")
    rcumode = rcuarg[3]
    calrunduration = durarg[3:]
    acc2bstvars = {}
    if datformat == 'hdf':
        hf = h5py.File(anacc2bstfilepath, 'r')
        acc2bstvars['XX'] = hf['XX']
        acc2bstvars['XY'] = hf['XY']
        acc2bstvars['YY'] = hf['YY']
        acc2bstvars['times'] = hf['timeaccstart']
        bst_pols = (hf['XX'], hf['XY'], hf['YY'])
        filestarttimes = hf['timeaccstart']
        freqs = hf['frequency']
        # Use calrunstarttime in filename
        # Use calrunduration in filename
        # Use calsrc in filename
        caltab_id = hf.attrs['calibrationTableDate']
        used_autocorr = hf.attrs['use_ac']
    else:
        acc2bstfilename = '_'.join((begin_utc_str, acc2bst, rcuarg, calsrc,
                                    durarg, caltabdate))
        acc2bstfilepath = acc2bstfiledir + '/' + acc2bstfilename
        bst_xx = numpy.load(acc2bstfilepath + '_XX' + '.npy')
        bst_xy = numpy.load(acc2bstfilepath + '_XY' + '.npy')
        bst_yy = numpy.load(acc2bstfilepath + '_YY' + '.npy')
        bst_pols = (bst_xx, bst_xy, bst_yy)
        filestarttimes = numpy.load(acc2bstfilepath + '_times' + '.npy')
        freqs = numpy.load(acc2bstfilepath + '_freqs' + '.npy')
        # Use calrunstarttime in filename
        # Use calrunduration in filename
        # Use calsrc in filename
        caltab_id = ""
        used_autocorr = None
    return (bst_pols, filestarttimes, freqs, calrunstarttime, calrunduration,
            calsrc, caltab_id, stnid, used_autocorr)


def saveacc2bst(bst_pols, filestarttimes, freqs, calrunstarttime,
                calrunduration, calsrc, caltab_id, stnid,
                used_autocorr, saveformat="hdf5"):
    """\
    Save acc2bst data to file.

    Dataformat can be hdf or numpy.

    Parameters
    ----------
    bst_pols: tuple of floatarray
        Tuple contains polarized BST-like data: (bst_xx, bst_xy, bst_yy).
    filestarttimes: list of datetime
        List of start times of the file partitioned data.
    freqs: list of floats
        List of frequencies.
    calrunstarttime: datetime
        Start time of this calibration run.
    calrunduration: float
        Duration of this calibration run.
    calsrc: str
        Source pointed to and which is intended as a calibration source.
    caltab_id: str
        Calibration table ID used for calibrating the input visibilities.
    stnid: str
        ID string of the station used.
    used_autocorr: bool
        Was autocorrelation used in beamforming or not?
    saveformat : str
        Format to use for saving.

    See Also
    --------
    readacc2bst: converse of this function.
    """
    (bstXX, bstXY, bstYY) = bst_pols
    version = '6'  # Version of this dataformat
    calrundurationstr = str(int(calrunduration.total_seconds()))
    # Calculate start of ACC run.
    # Form self describing filename.
    dtlabel = 'acc2bst'
    _sb, nqz = modeparms.freq2sb(freqs[0])
    rcumode = modeparms.nqz2rcumode(nqz)
    acc2bstbase = "{}_{}_spw{}_{}_dur{}_ct{}_v{}_{}".format(
        stnid, calrunstarttime.strftime("%Y%m%dT%H%M%S"), rcumode, calsrc,
        calrundurationstr, caltab_id, version, dtlabel)
    pntstr = ilisa.operations.directions.normalizebeamctldir(calsrc)
    # Write out the data.
    if saveformat == 'hdf5':
        hf = h5py.File(acc2bstbase + ".hdf5", "w")
        hf.attrs['DataDescription'] = 'LOFAR acc2bst data'
        hf.attrs['StationID'] = stnid
        hf.attrs['calibrationSource'] = calsrc
        hf.attrs['pointing'] = pntstr
        hf.attrs['ObservationStart'] = calrunstarttime.isoformat()
        hf.attrs['ObservationDuration'] = calrundurationstr
        hf.attrs['calibrationTableDate'] = caltab_id
        hf.attrs['version'] = version
        hf.attrs['use_ac'] = used_autocorr
        hf['frequency'] = freqs
        hf['frequency'].attrs['unit'] = "Hz"
        hf['timeaccstart'] = filestarttimes.view('<i8')
        hf['timeaccstart'].attrs['unit'] = "s"

        hf['XX'] = bstXX
        hf['XX'].attrs['unit'] = "arb. power"
        hf['XY'] = bstXY
        hf['XY'].attrs['unit'] = "arb. complex power"
        hf['YY'] = bstYY
        hf['YY'].attrs['unit'] = "arb. power"

        hf['XX'].dims.create_scale(hf['timeaccstart'])
        hf['XX'].dims.create_scale(hf['frequency'])
        hf['XY'].dims.create_scale(hf['timeaccstart'])
        hf['XY'].dims.create_scale(hf['frequency'])
        hf['YY'].dims.create_scale(hf['timeaccstart'])
        hf['YY'].dims.create_scale(hf['frequency'])
        hf['XX'].dims[0].attach_scale(hf['timeaccstart'])
        hf['XX'].dims[1].attach_scale(hf['frequency'])
        hf['XY'].dims[0].attach_scale(hf['timeaccstart'])
        hf['XY'].dims[1].attach_scale(hf['frequency'])
        hf['YY'].dims[0].attach_scale(hf['timeaccstart'])
        hf['YY'].dims[1].attach_scale(hf['frequency'])
        hf.close()
    else:
        numpy.save(acc2bstbase + '_times', filestarttimes)
        numpy.save(acc2bstbase + '_freqs', freqs)
        numpy.save(acc2bstbase + '_XX', bstXX)
        numpy.save(acc2bstbase + '_XY', bstXY)
        numpy.save(acc2bstbase + '_YY', bstYY)
    return acc2bstbase + "." + saveformat


def viewbst(bstff, pol_stokes=True, printout=False, update_wait=False):
    """\
    View BST data

    Parameters
    ----------
    bstff : str
        Path to BST file-folder.
    pol_stokes : bool
        Use Stokes representation for polarization?
    printout : bool
        Just print-out data rather than plot it?
    """
    (bst_datas_x, bst_datas_y, bst_datas_xy, ts_list, freqs, obsinfo
            ) = readbstfolder(bstff)
    stnid = obsinfo['station_id']
    starttime = obsinfo['datetime']
    intg = obsinfo['integration']
    dur_tot = float(obsinfo['duration_scan'])
    pointing = obsinfo['pointing']
    max_nr_bls = obsinfo['max_nr_bls']

    # Get nrfiles and nrsamps per file:
    nrfiles = len(ts_list)
    nrsamps = len(ts_list[0])
    print('# nrfiles:', nrfiles)
    print('# nrsamps:', nrsamps)

    # Squash list of data arrays (no padding between files)
    file_dur = (ts_list[0][-1]-ts_list[0][0]).total_seconds()
    ts = numpy.concatenate(ts_list)
    dur_cur = (ts[-1]-ts[0]).total_seconds()
    t_left = dur_tot - dur_cur
    if t_left > 0:
        print("Data duration {}s less than requested. (Rec not finished?)"
              .format(t_left))
    bst_data_x = numpy.concatenate(bst_datas_x).reshape(-1, max_nr_bls).T
    bst_data_y = numpy.concatenate(bst_datas_y).reshape(-1, max_nr_bls).T
    if bst_datas_xy:
        bst_data_xy =  numpy.concatenate(bst_datas_xy).reshape(-1, max_nr_bls).T
        data2view_pq_name, data2view_pq = 'Abs(XY)', numpy.abs(bst_data_xy)
        data2view_qp_name, data2view_qp = 'Arg(XY)', numpy.rad2deg(numpy.angle(bst_data_xy))
        data2view_pq_unit = 'Flux [arb. units]'
        data2view_qp_unit = 'Phase [deg]'
        norm_pq, norm_qp = colors.LogNorm(), None
        cmap_pq, cmap_qp = None, 'hsv'
    data2view_pp_name, data2view_pp = 'X-pol', bst_data_x
    data2view_qq_name, data2view_qq = 'Y-pol', bst_data_x
    data2view_qq_unit, data2view_qq_unit = 'Flux [arb. units]', 'Flux [arb. units]'
    norm_pp, norm_qq = colors.LogNorm(), colors.LogNorm()
    cmap_pp, cmap_qq = None, None
    if pol_stokes:
        # Stokes I
        data2view_pp_name = 'Stokes I'
        data2view_pp = bst_data_x + bst_data_y
        cmap_pp = None
        norm_pp = colors.LogNorm()

        # Stokes Q
        data2view_qq_name = '(antenna) Stokes Q'
        data2view_qq = bst_data_x - bst_data_y
        data2view_qq_unit = 'Signed flux [arb. units]'
        cmap_qq = 'RdBu_r'
        norm_qq = colors.SymLogNorm(linthresh=1e5, base=10)

        if bst_datas_xy:
            # Stokes U
            data2view_pq_name = 'Stokes U'
            data2view_pq = numpy.real(bst_data_xy + numpy.conj(bst_data_xy))
            data2view_pq_unit = 'Signed flux [arb. units]'
            cmap_pq = 'RdBu_r'
            norm_pq = colors.SymLogNorm(linthresh=1e5, base=10)
            # Stokes V
            data2view_qp_name = 'Stokes V'
            data2view_qp = numpy.imag(bst_data_xy - numpy.conj(bst_data_xy))
            data2view_qp_unit = 'Signed flux [arb. units]'
            cmap_qp = 'RdBu_r'
            norm_qp = colors.SymLogNorm(linthresh=1e5, base=10)
        stokes_norm = False
        if stokes_norm:
            data2view_qq = data2view_qq / data2view_pp
            data2view_q_name = '(antenna) Stokes q'
            data2view_q_unit = 'Signed relative flux []'
            norm_qq = None  # colors.SymLogNorm(linthresh=1e-3, base=10)
    if not printout:
        _xtr_corrs = 1
        if bst_datas_xy:
            _xtr_corrs = 2
        fig, axs = plt.subplots(2, _xtr_corrs, sharex=True, sharey=True)
        if _xtr_corrs == 1:
            (ax_pp, ax_qq) = axs
        else:
            ((ax_pp, ax_pq), (ax_qq, ax_qp)) = axs
        fr_grd, ts_grd = numpy.meshgrid(freqs/1e6, ts)

        # Plot quantity pp:
        bstplt_pp = ax_pp.pcolormesh(ts_grd, fr_grd, data2view_pp.T,  # numpy.log10(data2view_pp.T),
                                     cmap=cmap_pp, norm=norm_pp,
                                     shading='nearest')
        cbar_pp = fig.colorbar(bstplt_pp, ax=ax_pp)
        cbar_pp.set_label('Flux [arb. units]')
        ax_pp.set_ylabel('Frequency [MHz]')
        ax_pp.set_title('{}'.format(data2view_pp_name))

        # Plot quantity qq:
        bstplt_qq = ax_qq.pcolormesh(ts_grd, fr_grd, data2view_qq.T,
                                     cmap=cmap_qq, norm=norm_qq,
                                     shading='nearest')
        cbar_qq = fig.colorbar(bstplt_qq, ax=ax_qq)
        cbar_qq.set_label(data2view_qq_unit)
        ax_qq.set_title('{}'.format(data2view_qq_name))
        fig.autofmt_xdate()

        ax_qq.xaxis.set_major_formatter(mdates.DateFormatter('%H:%M:%S'))
        ax_qq.set_xlabel('Datetime [UT]  Starts: {}'.format(starttime))
        ax_qq.set_ylabel('Frequency [MHz]')

        if bst_datas_xy:
            # Plot quantity pq:
            bstplt_pq = ax_pq.pcolormesh(ts_grd, fr_grd, data2view_pq.T,
                                         cmap=cmap_pq, norm=norm_pq,
                                         shading='nearest')
            cbar_pq = fig.colorbar(bstplt_pq, ax=ax_pq)
            cbar_pq.set_label(data2view_pq_unit)
            ax_pq.set_title('{}'.format(data2view_pq_name))

            # Plot quantity qp:
            bstplt_qp = ax_qp.pcolormesh(ts_grd, fr_grd, data2view_qp.T,
                                         cmap=cmap_qp, norm=norm_qp,
                                         shading='nearest')
            cbar_qp = fig.colorbar(bstplt_qp, ax=ax_qp)
            cbar_qp.set_label(data2view_qp_unit)
            ax_qp.set_title('{}'.format(data2view_qp_name))

        supertitle = ('{} BST intg: {}s dur: {}s'.format(stnid, intg, dur_cur)
                      + ' pointing: {}'.format(pointing))
        plt.suptitle(supertitle)
        if update_wait:
            plt.draw()
            print("Waiting for data update (press any key to unblock)")
            plt.waitforbuttonpress(file_dur)
            plt.close()
        else:
            plt.show()
    else:
        # CSV style:
        #   Header
        t_prev = ts[0]
        print("## H:M:S since {} UT".format(t_prev.isoformat()), "Freq[MHz]",
              data2view_pp_name, data2view_qq_name, sep=' ')
        #   Data
        for ti, t in enumerate(ts):
            for freqi, freq in enumerate(freqs):
                dataval_p, dataval_q = (data2view_pp[freqi, ti],
                                        data2view_qq[freqi, ti])
                del_t = t - t_prev
                print(del_t, freq/1e6, dataval_p, dataval_q, sep=' ')


def viewsst(sstff, freqreq, sample_nr=None, rcu_sel=None, printout=False):
    """\
    View SST data

    If printout is False it produces a plot otherwise it prints data to stdout.

    SST data vary over frequency, time sample, and RCU. Various types of plots
    of it can be produced:
      * *persbs*: per frequency, plot over time waterfall of RCUs
      * *mean*: mean dynamic spectra averaged over RCU
      * *dynspec*: dynamic spectra averaged over selected RCUs
      * *overlay*: overlay spectrum of given RCU slice
      * *ssmosaic*: snapshot mosaic at time samp of spectra for all RCUs

    Which plot is selected depends on argument settings (T,F,X = Set, Unset,
    Arbitrary), see Table:
    | freqreq | sample_nr | rcu_sel | Plot     |
    |---------|-----------|---------|----------|
    | T       | X         | X       | persb    |
    | F       | F         | F       | mean     |
    | F       | F         | T       | dynspec  |
    | F       | T         | T       | overlay  |
    | F       | T         | F       | ssmosaic |

    Parameters
    ----------
    sstff : str
        Path to SST file-folder.
    freqreq : float
        Requested frequency in Hz.
    sample_nr : int
        Sample number.
    rcu_sel : str
        RCU number or range of numbers.
    printout : bool
        Print out data instead of plotting it.
    """
    sstdata_rcu, ts_list, freqs, obsinfo = readsstfolder(sstff)
    starttime = obsinfo['datetime']

    # Squash file_nr and in file intg index to just samples
    sstdata = numpy.array(sstdata_rcu).reshape((modeparms.nrofrcus, -1,
                                                modeparms.TotNrOfsb))
    nrrcus, nrsamps, nrsbs = sstdata.shape
    nrfiles = len(sstdata_rcu[0])
    sampsperfile = sstdata_rcu[0][0].shape[0]
    print('# nrfiles:', nrfiles)
    print('# samps/file', sampsperfile)
    print('# nrsamps:', nrsamps)
    print('# nrrcus:', nrrcus)
    print('# nrsbs:', nrsbs)
    ts = numpy.ravel(ts_list)
    sbreq = None
    if freqreq:
        sbreq = int(numpy.argmin(numpy.abs(freqs-freqreq)))
    if rcu_sel is not None:
        if ':' in rcu_sel:
            rcu_sel = rcu_sel.split(':')
            rcu_sel = slice(*[int(_a) for _a in rcu_sel])
        else:
            rcu_sel = slice(int(rcu_sel), int(rcu_sel) + 1)
    if printout:
        if rcu_sel is None:
            rcus_str = ' '.join(
                ['P_rcu'+str(_rcu) for _rcu in range(modeparms.nrofrcus)])
        else:
            rcus_str = 'P_rcu' + str(rcu_sel)
        print('## UT Freq[Hz] ' + rcus_str)
        for tidx, t in enumerate(ts):
            if sample_nr is not None:
                if tidx != sample_nr:
                    continue
            for frqidx, freq in enumerate(freqs):
                if sbreq is not None:
                    if sbreq != frqidx:
                        continue
                print(modeparms.astimestr(t), freq, end=' ')
                _out = sstdata[:, tidx, frqidx]
                if rcu_sel is not None:
                    _out = [_out[rcu_sel]]
                print(*_out)
        return
    # Plot
    if freqreq:
        show = 'persb'
    else:
        if sample_nr is None and rcu_sel is None:
            show = 'mean'
        if sample_nr is None and rcu_sel is not None:
            show = 'dynspec'
        if sample_nr is not None and rcu_sel is not None:
            show = 'overlay'
        if sample_nr is not None and rcu_sel is None:
            show = 'ssmosaic'
    if show == 'mean':
        # Show mean over RCUs
        meandynspec = numpy.mean(sstdata, axis=0)
        res = meandynspec
        if res.shape[0] > 1:
            plt.pcolormesh(freqs/1e6, ts, res, norm=colors.LogNorm(),
                           shading='nearest')
            plt.colorbar()
            plt.title('Mean (over RCUs) dynamicspectrum\n'
                      + 'Starttime: {} Station: {}'
                      .format(starttime, obsinfo['station_id']))
            plt.xlabel('Frequency [MHz]')
            plt.ylabel('Time [h]')
        else:
            # Only one integration so show it as 2D spectrum
            plt.plot(freqs/1e6, res[0, :])
            plt.yscale('log')
            plt.xlabel('Frequency [MHz]')
            plt.ylabel('Power [arb. unit]')
    elif show == 'dynspec':
        # Show dynamic spectrum of RCU
        dynspec = numpy.mean(sstdata[rcu_sel], axis=0)
        res = dynspec
        if res.shape[0] > 1:
            plt.pcolormesh(freqs/1e6, ts, res, norm=colors.LogNorm(),
                           shading='nearest')
            plt.colorbar()
            plt.title(('Dynamicspectrum of RCU {}\n'
                      + 'Starttime: {} Station: {}')
                      .format(rcu_sel, starttime, obsinfo['station_id']))
            plt.xlabel('Frequency [MHz]')
            plt.ylabel('Time [h]')
        else:
            # Only one integration so show it as 2D spectrum
            plt.plot(freqs/1e6, res[0, :])
            plt.yscale('log')
            plt.xlabel('Frequency [MHz]')
            plt.ylabel('Power [arb. unit]')
    elif show == 'persb':
        ampVStime = True
        res = sstdata[:, :, sbreq]
        resX = res[0::2, :]
        resY = res[1::2, :]
        plt.subplot(211)
        if ampVStime:
            plt.plot(ts, numpy.transpose(resX))
        else:
            plt.pcolormesh(ts, numpy.arange(96), resX, norm=colors.LogNorm())
        plt.title('X pol')
        plt.subplot(212)
        if ampVStime:
            plt.plot(ts, numpy.transpose(resY))
            plt.gcf().autofmt_xdate()
        else:
            plt.pcolormesh(ts, numpy.arange(96), resY, norm=colors.LogNorm())
            plt.ylabel('rcu [nr]')
        plt.xlabel('Time [UT]')
        plt.title('Y pol')
        plt.suptitle('Freq {} MHz'.format(freqs[sbreq]/1e6))
    elif show == 'ssmosaic':
        axdim1 = 8
        axdim0 = 192 // axdim1 // 2
        _fig, axs = plt.subplots(axdim0, axdim1)
        for rcu_nr in range(0, 192, 2):
            sbplt_nr = rcu_nr // 2
            # Plot X
            idx_axs = numpy.unravel_index(sbplt_nr, (axdim0, axdim1))
            axs[idx_axs].semilogy(freqs, sstdata[rcu_nr + 0, sample_nr, :],
                                  label='X')
            # & Y in same subplot
            axs[idx_axs].semilogy(freqs, sstdata[rcu_nr + 1, sample_nr, :],
                                  label='Y')
            axs[idx_axs].set_title('{},{} (Ant{})'.format(rcu_nr, rcu_nr + 1,
                                   modeparms.rcu2antpol(rcu_nr)[0]), y=0.65,
                                   fontsize=8)
            axs[idx_axs].legend(fontsize=5)  # prop={'size': 6})
        plt.suptitle('RCU spectra @ {} UT, station {}'.format(
            ts[sample_nr], obsinfo['station_id']))
    elif show == 'overlay':

        res = sstdata[rcu_sel, sample_nr, :].squeeze()
        plt.semilogy(freqs, numpy.transpose(res))
        rcus = range(rcu_sel.start, rcu_sel.stop)
        plt.legend(rcus)
        plt.title('RCU spectrum @ {} UT, station {}'.format(
            ts[sample_nr], obsinfo['station_id']))
    plt.show()


def plotcmplxmat(cm, cmplxrep='ReIm', xylabels='', title='Complex matrix'):
    """\
    Plot a complex matrix

    Parameters
    ----------
    cm : 2D array
        Complex matrix to plot
    cmplxrep : str
        Representation to use for complex matrix. Default is 'ReIm' which
        uses real & imaginary parts, 'AbsArg' or  'polar' plots absolute and
        arguments, and 'HermSD' is useful to check the hermitianity of the
        matrix by plotting the real and imaginary parts of the sum and difference
        of the matrix with its hermitian transpose.
    xylabels : str
        Labels to use for axes.
    """
    cm = numpy.asmatrix(cm)
    mats = (numpy.real(cm), numpy.imag(cm))
    ptit = ('Re', 'Im')
    cmaps = ('coolwarm', 'coolwarm')
    nrcomps = 2
    if cmplxrep == 'AbsArg' or cmplxrep == 'polar':
        # Use AbsArg
        mats = (numpy.abs(cm), numpy.angle(cm))
        ptit = ('Abs', 'Arg')
        cmaps = ('viridis', 'hsv')
    elif cmplxrep == 'HermSD':
        nrcomps = 4
        mats = (numpy.real(cm+cm.H), numpy.imag(cm+cm.H),
                numpy.real(cm-cm.H), numpy.imag(cm-cm.H))
        ptit = ('Re(mat+mat.H)','Im(mat+mat.H)', 'Re(mat-mat.H)', 'Im(mat-mat.H)')
        cmaps = ('viridis', 'viridis', 'viridis', 'viridis')
    fig, axs = plt.subplots(nrcomps//2, 2, sharey=True)
    ar_is = range(nrcomps) if nrcomps//2 == 1 else numpy.ndindex(nrcomps//2, 2)
    for cmp_i, ar_i in enumerate(ar_is):
        spi = axs[ar_i].matshow(mats[cmp_i], cmap=cmaps[cmp_i])
        axs[ar_i].set_xlabel(xylabels)
        axs[ar_i].set_ylabel(xylabels)
        axs[ar_i].set_title(ptit[cmp_i])
        plt.colorbar(spi, ax=axs[ar_i])
    fig.suptitle(title)


def viewxst(xstsampdata, samptime, freq, titleingress, printout=False,
            poltype=None, cmplxrep='ReIm'):
    """
    View XST data

    Parameters
    ----------
    xstff : str
        Path to a CVC data filefolder.

    printout : bool
        Should the selected data be printed out, otherwise it will be plotted.
    poltype : str
        Specifies how polarized data should be represented. Choices are:
        'sto' for Stokes I,Q,U,V; 'lin' for XX,XY,YY; 'mix' for one matrix
        with both X,Y mixed in alternating indices.
    cmplxrep :  str
        Specifies how complex data should be represented. Choices are:
        'ReIm' for cartesian; 'AbsArg' for polar representation.
    """
    cvpol = cov_flat2polidx(xstsampdata)

    title = """{}, freq={} MHz @ {} UT""".format(
        titleingress, round(freq / 1e6, 2), samptime)

    if not printout:
        colorscale = None  # Colorscale for xst data plot (default None)
        if colorscale == 'log':
            normcolor = colors.LogNorm()
        else:
            normcolor = None
    else:
        import sys

    if poltype == 'sto':
        # Avoiding dependence on dreambeam's convertxy2stokes():
        cv_serial = {'I': +1* (cvpol[0, 0] + cvpol[1, 1]),
                     'Q': +1* (cvpol[0, 0] - cvpol[1, 1]),
                     'U': +1* (cvpol[0, 1] + cvpol[1, 0]),
                     'V': -1j*(cvpol[0, 1] - cvpol[1, 0])}
    elif poltype == 'lin':
        cv_serial = {"XX'": cvpol[0][0],
                     "XY'": cvpol[0][1],
                     "YX'": cvpol[1][0],
                     "YY'": cvpol[1][1]}
    else:
        cv_serial = {'Cov': xstsampdata}
    if printout: print("# " + title)
    for polk in cv_serial:
        if printout:
            print(polk)
            numpy.savetxt(sys.stdout, cv_serial[polk], fmt='%.6e')
        else:
            plotcmplxmat(cv_serial[polk], cmplxrep=cmplxrep,
                         title=title+"\n"+polk, xylabels='RCU [#]')
    if not printout:
        print("Kill plot window for next plot...")
        plt.show()
    return cv_serial


def latest_scanrec_path():
    """\
    Yield latest scanrec path

    Yields
    ------
    latest_scanrecpath : str
        Path on DRU of latest ScanRec data.
    """
    latestdatafile = sorted(filter(lambda _f: _f.startswith('latest'),
                                   os.listdir(USER_CACHE_DIR)),
                            key=os.path.getmtimeget)[-1]
    with open(latestdatafile, 'r') as f:
        lines_in = f.readlines()
    _contents = [_l.rstrip() for _l in lines_in]
    if _contents[0] != 'ONGOING':
        contents = _contents
        yield from contents
    else:
        scanrec_nr = None
        while True:
            with open(latestdatafile, 'r') as f:
                lines_in = f.readlines()
            _contents = [_l.rstrip() for _l in lines_in]
            if _contents[0] == 'ONGOING':
                contents = _contents[1:]
            else:
                contents = _contents
            if scanrec_nr is None and len(contents) == 0:
                scanrecpath = None
            elif scanrec_nr is None and len(contents) > 0:
                # Start with 1st scanrec
                scanrec_nr = 0
                scanrecpath = contents[scanrec_nr]
            elif scanrec_nr+1 < len(contents):
                # Goto next scanrec
                scanrec_nr += 1
                scanrecpath = contents[scanrec_nr]
            elif scanrec_nr+1 == len(contents):
                # No new datafiles yet (keep last scanrecpath)
                # scanrecpath = None
                pass
            if scanrecpath is None and _contents[0] != 'ONGOING':
                break
            yield scanrecpath


def view_bsxst(dataff, filenr, sampnr, freq, printout=False, poltype=None,
               cmplxrep=None):
    """\
    View BST, SST, XST statistics data files

    Nominal use-case is plotting data after completed recording.
    It this case, a file-folder is given, which could be one of
    the 'STatistics' data: BST, SST or XST. Each of these
    would then be displayed in a matter appropriate to it.

    BST is plotted as beamlets versus time with all (chunked) files collated,
    albeit with time gaps.

    SST is plotted as RCU versus time with all (chunked) files collated,
    albeit with time gaps.

    XST is plotted as a sequence of plots of RCU versus RCU per sampnr
    per filenr.

    This function can also be used to plot real-time data.

    Parameters
    ----------
    dataff : str
        Data file-folder containing the bst, sst or xst.
    filenr : str
        File number selection.
    sampnr : str
        Sample number selection.
    freq : float
        Frequency selection.
    printout : bool
        Print out data rather than plot.
    poltype : str
        Format for polarization: 'sto', 'lin' or 'mix'.
    cmplxrep : str
        Complex representation of complex values: 'ReIm' or 'AbsArg',
        for real-imaginary (cartesian) or absolute-argument (polar)
        representation.
    """
    if not dataff:
        dataff = latest_scanrec_path()
        print('dataff', next(dataff))
    dataff = os.path.normpath(dataff)
    scnrecinfo = ScanRecInfo().read_scanrec(dataff)
    if scnrecinfo and scnrecinfo.comments:
        print('# comments:\n', scnrecinfo.comments)
    lofar_datatype = datafolder_type(dataff)
    if lofar_datatype == 'sst':
        rcunr = filenr
        filenr = None
    filenrs = [0]
    if filenr is not None:
        if ':' in filenr:
            filenrs = range(*map(int, filenr.split(':')))
        elif filenr.isdigit():
            filenrs = [int(filenr)]

    sampnrs = [0]
    if sampnr is not None:
        if ':' in sampnr:
            sampnrs = range(*map(int, sampnr.split(':')))
        elif sampnr.isdigit():
            sampnr = int(sampnr)
            sampnrs = [sampnr]

    for filenr in filenrs:
        for sampidx in sampnrs:
            if dataff:
                dataffs = [dataff]
            else:
                dataffs = latest_scanrec_path()
            for dataff in dataffs:
                if not dataff:
                    print('Waiting 1s for data...')
                    time.sleep(1.0)
                    continue
                dataff = os.path.normpath(dataff)
                lofar_datatype = datafolder_type(dataff)
                if printout:
                    # Provide a sort of header
                    print('# datafile:', dataff)
                    if lofar_datatype == 'sst':
                        print('# rcunr:', rcunr)
                    else:
                        print('# filenr:', filenr)
                    print('# sampnr:', sampidx)
                    print('# poltype:', poltype)
                if not freq:
                    if lofar_datatype == 'acc':
                        freq = 0.0
                        if sampidx:
                            freq = modeparms.sb2freq(sampidx, 0)  # Assume LBA band
                else:
                    freq = float(freq)
                if (lofar_datatype == 'bst' or lofar_datatype == 'bst-357'
                        or lofar_datatype == 'bstc'):
                    _pol_stokes = False
                    if poltype == 'sto':
                        _pol_stokes = True
                    viewbst(dataff, pol_stokes=_pol_stokes,
                            printout=printout)
                elif lofar_datatype == 'sst':
                    viewsst(dataff, freq, sampnr, rcunr, printout)
                elif lofar_datatype == 'xst' or lofar_datatype == 'xst-SEPTON' \
                        or lofar_datatype=='acc':
                    # Get selected data
                    xstobj = CVCfiles(dataff)
                    nrfiles = xstobj.getnrfiles()
                    nrsamps = len(xstobj.samptimeset[0])
                    print('# nrfiles', nrfiles)
                    print('# nrsamps', nrsamps)
                    xstfiledata = xstobj[filenr]
                    xstsampdata = xstfiledata[sampidx]

                    # Get metadata
                    obs_ids = xstobj.scanrecinfo.get_obs_ids()
                    times_in_filetimes = xstobj.samptimeset[filenr]
                    ldatinfo = xstobj.scanrecinfo.ldatinfos[obs_ids[filenr]]
                    samptime = times_in_filetimes[sampidx]
                    freq = ldatinfo.get_recfreq(sampidx)
                    stnid = xstobj.scanrecinfo.stnid
                    sourcename = None
                    if scnrecinfo:
                        sourcename = scnrecinfo.sourcename
                    titleingress = "LOFAR {}, src {}".format(stnid, sourcename)

                    xstdata = viewxst(xstsampdata, samptime, freq, titleingress,
                                      printout, poltype, cmplxrep)
                else:
                    raise RuntimeError("Not a bst, sst, or xst filefolder")


def export_ldat(dataff):
    """Export LOFAR data files

    Parameters
    ----------
    dataff : str
        name of filefolder with lofar data
    save : bool
        If True, save data to numpy formated file (filename derived from dataff)

    Returns
    -------
    data_arrs : tuple of arrays
        Tuple of complex numpy array pf data from dataff. Tuple index represents
        file number ('filenr'). Array indices depends on
        lofar data-type:
          'bst': ['sampnr', 'sbnr', 'polpolnr']
          'sst': ['sampnr', 'sbnr', 'rcunr']
          'xst' | 'acc': ['sampnr', 'polnr', 'polnr', 'antnr', antnr']
    metadata : dict
        Dataset dict. Consists of
        'ID_scanrec' : str
            ID string for observation, equal to the ldat filefolder name.
        'datatype': str
            The ldat data-type: ['bst', 'sst', 'xst', 'acc'].
        'positions' : (nrants, 3) float array
            The element layout absolute positions w.r.t. ITRF.
        'stn_rot' : (3, 3) float array
            The rotation matrix between the ITRF and the station local frame.
            It is used like this:
                'positions_stn' = 'positions' @ 'stn_rot'.
        'start_datetime' : numpy.datetime64
            The datetime of the first data sample.
        'delta_secs' : float array
            The delta-time in seconds of the sample w.r.t the start_datetime.
            Has shape ['filenr', 'sampnr'].
        'frequencies' : float array
            The center frequency in Hertz of the observation subband.
            Can either be a scalar (for observation at fixed frequency)
            or have shape ['filenr', 'sampnr'] like 'delta_secs'.
        'station' : str
            ID of station.
        'pointing' : str
            String format of 3-tuple of pointing direction: '<azi>,<elv>,<ref>'.
    """
    lofardatatype = datafolder_type(dataff)
    if not lofardatatype:
        raise TypeError('Cannot export unknow LOFAR data-type: {}'
                        .format(lofardatatype))
    id_scanrec = os.path.basename(os.path.normpath(dataff))
    stn_rot = None
    if lofardatatype == 'bst':
        (bst_dat_xx, bst_dat_yy, bst_dat_xy, ts_list, freqs, obsinfo
         ) = readbstfolder(dataff)
        bst_dat = [bst_dat_xx, bst_dat_yy]
        if bst_dat_xy:
            bst_dat.append(bst_dat_xy)
        data_arr = numpy.asarray(bst_dat)
    elif lofardatatype == 'sst':
        sstdata_rcu, ts_list, freqs, obsinfo = readsstfolder(dataff)
        data_arr = numpy.asarray(sstdata_rcu)
    elif lofardatatype == 'xst' or lofardatatype == 'acc':
        cvcobj = CVCfiles(dataff)
        cvc_array = cvcobj.as_array()
        data_arr = cov_flat2polidx(cvc_array)
        positions = cvcobj.get_positions_ITRF()
        ts_list = cvcobj.samptimeset
        freqs = numpy.asarray(cvcobj.freqset)
        station_id = cvcobj.scanrecinfo.get_stnid()
        pointing = cvcobj.scanrecinfo.get_pointingstr()
        # sourcename  = cvcobj.scanrecinfo.sourcename
        stn_rot = cvcobj.stn_rot
    if lofardatatype == 'bst' or lofardatatype == 'sst':
        data_arr = numpy.moveaxis(data_arr, 0, -1)
        station_id = obsinfo['station_id']
        positions, _names, _xyzlbls \
            = antennafieldlib.get_tier_layouts(station_id,
                                               obsinfo['antennaset'][:3],
                                               coordsys='ITRF')
        pointing = obsinfo['pointing']
    ts = numpy.vectorize(numpy.datetime64)(ts_list)
    start_datetime = ts[0][0]
    delta_secs = (ts - start_datetime) / numpy.timedelta64(1, 's')
    if numpy.all(freqs == numpy.take(freqs, 0)):
        # If all frequencies are the same then use only that frequency
        # (mainly for single freq XST)
        freqs = numpy.take(freqs, 0)
    metadata = {'ID_scanrec': id_scanrec,
               'station': station_id,
               'datatype': lofardatatype,
               'positions': positions,
               'start_datetime': start_datetime,
               'delta_secs': delta_secs,
               'frequencies': freqs}
    # if sourcename:
    #     dataset['source'] = sourcename
    if pointing:
        metadata['pointing'] = pointing
    if stn_rot is not None:
        metadata['stn_rot'] = stn_rot
    return tuple(data_arr), metadata


def cli_export():
    parser = argparse.ArgumentParser(
        description="Export lofar data to other formats",
        formatter_class=argparse.RawTextHelpFormatter
    )
    parser.add_argument('dataff', nargs='?', default=None,
                        help="acc, bst, sst or xst filefolder")
    parser.add_argument('-d', '--dataformat', default='npz',
                        help="""\
Output data-formats:
    * 'npz'
        The npz data consists of observational meta-data
        and data paylooad. In particular:
            * 'ID_scanrec': the ID of observatino.
            * 'station': Name of station.
            * 'datatype': lofar datatype, can be 'ACC', 'BST', 'SST' or 'XST'.
            * 'positions': Absolute ITRF positions of layout in meters.
            * 'start_datetime': Start of observation as UTC (numpy.datetime64)
            * 'delta_secs': Time delta from start of sample in seconds.
            * 'frequencies': frequencies in Hertz.
            * 'arr_'+<arr_nr>: Observed data split into indexed chunks.                    
"""
)
    args = parser.parse_args()
    data_arrs, metadata = export_ldat(args.dataff)
    if args.dataformat == 'npz':
        numpy.savez_compressed(metadata['ID_scanrec'], *data_arrs, **metadata)


def cli_view():
    parser = argparse.ArgumentParser(description="Plot or print lofar data")
    # filenr : int
    #     Selects data file based on its ordinal number during recording.
    # sampnr : int
    #     Ordinal number of the record sample within the filenr-th file.
    parser.add_argument('-n', '--filenr', type=str, default=None,
                        help="Can be file # or RCU #, or also a range #:#")
    parser.add_argument('-s', '--sampnr', type=str, default=None,
                        help="Sample # or a range #:#")
    parser.add_argument('-f', '--freq', type=float, default=None,
                        help='Frequency in Hz')
    parser.add_argument('-p', '--pol', type=str, default=None,
                        help='Set pol. format: "lin", "sto", def: "None"')
    parser.add_argument('-c', '--cmplx', type=str, default='ReIm',
                        help='Complex format: "ReIm" or "AbsArg", def: "ReIm"')
    parser.add_argument('-o', '--printout', action="store_true",
                        help='Print out data to stdout, else plot')
    parser.add_argument('dataff', nargs='?', default=None,
                        help="acc, bst, sst or xst filefolder")
    args = parser.parse_args()

    view_bsxst(args.dataff, args.filenr, args.sampnr, args.freq, args.printout,
               args.pol, args.cmplx)


if __name__ == "__main__":
    cli_view()
