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
import numpy
import re
import datetime
import h5py
import yaml
import warnings
import ilisa
import ilisa.monitorcontrol.directions
import ilisa.antennameta.antennafieldlib as antennafieldlib
try:
    import dreambeam
    canuse_dreambeam = True
except ImportError:
    canuse_dreambeam = False
if canuse_dreambeam:
    from dreambeam.polarimetry import convertxy2stokes
    from dreambeam.polarimetry import cov_lin2cir


_RCU_SB_SEP = "+"

regex_ACCfolder = (
    r"^(?P<stnid>\w{5})_(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})"
    r"_(?P<hour>\d{2})(?P<minute>\d{2})(?P<second>\d{2})"
    r"_rcu(?P<rcumode>\d+)_dur(?P<duration_tot>\d+)(_(?P<calsrc>\w+))?_acc$")
regex_ACCfilename = (
    r"^(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})"
    r"_(?P<hour>\d{2})(?P<minute>\d{2})(?P<second>\d{2})"
    r"_acc_(?P<totnrsb>\d+)x(?P<nrrcu0>\d+)x(?P<nrrcu1>\d+)"
    r".dat$")
regex_xstfilename = (
    r"^(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})"
    r"_(?P<hour>\d{2})(?P<minute>\d{2})(?P<second>\d{2})"
    r"_rcu(?P<rcumode>\d+)_int(?P<integration>\d+)_dur(?P<duration_scan>\d+)"
    r"_dir(?P<RAint>\d+).(?P<RAdecimal>\d+),(?P<DECint>\d+).(?P<DECdecimal>\d+),(?P<ref>\s+)"
    r"_xst.dat$")


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
    >>> from ilisa.monitorcontrol.data_io import seqlists2slicestr
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


def obsfileinfo2filefolder(obsfileinfo):
    """\
    Convert obsfileinfo dict to filefolder name

    obsfileinfo:
        duration_tot:
        filenametime:
        integration:
        rcumode:
        ldat_type or datatype:
        mode:
        pointing:
        sb:
        station_id:

    Name format is
        <station_id>_<filenametime>_spw<rcumodes>_sb<sb>_int<integration>\
        _dur<duration_tot>[_dir<pointing>]_<ldat_type>
    """
    if obsfileinfo.get('datatype', None):
        ldat_type = obsfileinfo['datatype']
    else:
        ldat_type = obsfileinfo.get('ldat_type')
    filefoldername = '{}_{}'.format(obsfileinfo['station_id'],
                                    obsfileinfo['filenametime'])
    if ldat_type == "bst-357":
        filefoldername += "_rcu357"
    else:
        if obsfileinfo['rcumode'] != []:
            rcumodestr = \
                ''.join([str(rcumode) for rcumode in obsfileinfo['rcumode']])
        else:
            rcumodestr = str(obsfileinfo['mode'])
        filefoldername += "_spw" + rcumodestr
    if obsfileinfo['sb'] != [] and obsfileinfo['sb'] != '':
        filefoldername += "_sb"
        filefoldername += seqlists2slicestr(obsfileinfo['sb'])
    if 'integration' in obsfileinfo:
        filefoldername += "_int" + str(int(obsfileinfo['integration']))
    if 'duration_tot' in obsfileinfo:
        filefoldername += "_dur" + str(int(obsfileinfo['duration_tot']))
    if ldat_type != 'sst':
        if str(obsfileinfo['pointing']) != "":
            filefoldername += "_dir" + str(obsfileinfo['pointing'])
        else:
            filefoldername += "_dir,,"
    # filefoldername += "_" + obsfileinfo['source']
    # ldat_type extension
    filefoldername += "_" + ldat_type
    return filefoldername


def filefolder2obsfileinfo(filefolderpath):
    """Parse filefolder name and return an obsfileinfo"""
    filefolderpath = os.path.normpath(filefolderpath)
    filefoldername = os.path.basename(filefolderpath)
    # Format:
    # stnid_Ymd_HMS_spwstr_intstr_durstr_dirstr_[cal*]_acc
    # stnid?_Ymd_HMS_rcustr_sbstr_intstr_durstr_dirstr_bst
    # stnid?_Ymd_HMS_rcustr_intstr_durstr_sst
    # stnid?_Ymd_HMS_rcustr_sbstr_intstr_durstr_dirstr_xst
    filefoldersplit = filefoldername.split('_')
    # Take care of possible "cal*" tag just before ldattype:
    if filefoldersplit[-2].startswith('cal'):
        # Calibration applied to this ldat
        filefoldersplit.pop(-2)
    ldat_type = filefoldersplit.pop()
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
    if len(filefoldersplit) == 6:
        stnid = filefoldersplit.pop(0)
    else:
        stnid = None
    (Ymd, HMS, rcustr, intstr, durstr) = filefoldersplit[:]

    obsfileinfo = {}
    obsfileinfo['station_id'] = stnid
    obsfileinfo['filenametime'] = Ymd + '_' + HMS
    obsfileinfo['datetime'] = datetime.datetime.strptime(Ymd + 'T' + HMS,
                                                         '%Y%m%dT%H%M%S')
    obsfileinfo['rcumode'] = rcustr[3:]
    obsfileinfo['subbands'] = sbstr[2:]
    obsfileinfo['integration'] = int(intstr[3:])
    obsfileinfo['duration_scan'] = int(durstr[3:])
    obsfileinfo['pointing'] = dirstr[3:].split(',')
    obsfileinfo['ldat_type'] = ldat_type

    if len(obsfileinfo['rcumode']) > 1:
        obsfileinfo['rcumode'] = list(obsfileinfo['rcumode'])
    if _RCU_SB_SEP in obsfileinfo['subbands']:
        obsfileinfo['subbands'] = obsfileinfo['subbands'].split(
            _RCU_SB_SEP)

    if type(obsfileinfo['rcumode']) is not list:
        obsfileinfo['rcumode'] = [obsfileinfo['rcumode']]
    if type(obsfileinfo['subbands']) is not list:
        obsfileinfo['subbands'] = [obsfileinfo['subbands']]
    obsfileinfo['frequencies'] = numpy.empty(0)
    beamlets = []
    totnrsbs = 0
    for spw, rcumode in enumerate(obsfileinfo['rcumode']):
        sblist = modeparms.seqarg2list(obsfileinfo['subbands'][spw])
        nrsbs = len(sblist)
        sblo = sblist[0]
        sbhi = sblist[-1]
        nz = modeparms.rcumode2nyquistzone(rcumode)
        freqlo = modeparms.sb2freq(sblo, nz)
        freqhi = modeparms.sb2freq(sbhi, nz)
        obsfileinfo['frequencies'] = numpy.append(obsfileinfo['frequencies'],
                                                  numpy.linspace(freqlo,
                                                                 freqhi,
                                                                 nrsbs))
        bmltarg = seqlists2slicestr(
            ','.join([str(_b) for _b in range(nrsbs)]))
        beamlets.append(bmltarg)
        totnrsbs += nrsbs
    bits = 16

    if ldat_type == 'bst':
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
            obsfileinfo['frequencies'] = numpy.append(obsfileinfo['frequencies'],
                                                      numpy.linspace(freqlo,
                                                                     freqhi,
                                                                     nrsbs))
        obsfileinfo['max_nr_bls'] = maxnrbls

    # Assemble _cmds
    #    rcusetup_cmds
    rcusetup_cmds = modeparms.rcusetup_args2cmds(bits, 0)
    obsfileinfo['rcusetup_cmds'] = rcusetup_cmds
    #    beamctl_cmds
    beamctl_cmds = []
    for spw, rcumode in enumerate(obsfileinfo['rcumode']):
        band = modeparms.rcumode2band(rcumode)
        anadigdir = ','.join(obsfileinfo['pointing'])
        beamctl_cmd = modeparms.beamctl_args2cmds(beamlets[spw],
                                                  obsfileinfo['subbands'][spw],
                                                  band, anadigdir)
        beamctl_cmds.append(beamctl_cmd)
    obsfileinfo['beamctl_cmds'] = beamctl_cmds
    #     rspctl_cmds
    rspctl_cmds = modeparms.rspctl_stats_args2cmds(obsfileinfo['ldat_type'],
                                                   obsfileinfo['integration'],
                                                   obsfileinfo['duration_scan'],
                                                   obsfileinfo['subbands'])
    obsfileinfo['rspctl_cmds'] = rspctl_cmds

    return obsfileinfo


class ScanRecInfo(object):
    """This class maintains info on a scan recording (scanrec), which is one
    of the results of an iLiSA scan. One scanrec is a group of one or more
    files of a unique LOFAR station data product (ldat), i.e. acc, bfs, bst,
    sst or xst. The info in a ScanRecInfo object consists of sufficient
    parameters to redo the data, namely the stn_id, the iLiSA scanrec
    parameters, and a list of LDatInfo objects called obsinfos that maps to
    each ldat within the scanrec.
    """
    scanrecinfo_header = "SCANREC_INFO.yml"

    def __init__(self):
        self.headerversion = 2
        self.obsinfos = {}
        self.sourcename = ''

    def add_obs(self, obsinfo):
        """Add an LDatInfo object to this ScanRecInfo."""
        obs_id = obsinfo.filenametime
        self.obsinfos[obs_id] = obsinfo

    def get_obs_ids(self):
        """Get list of obs_ids.
        A obs_id is a key to the obsinfos list."""
        return sorted(self.obsinfos.keys())

    def set_stnid(self, stnid):
        self.stnid = stnid

    def get_stnid(self):
        try:
            stnid = self.stnid
        except:
            try:
                stnid = self.obsinfos[0].stnid
            except:
                try:
                    stnid = self.scanrecparms['stnid']
                except:
                    raise RuntimeError('Station id not found.')
        return stnid

    def set_sourcename(self, sourcename):
        self.sourcename = sourcename

    def set_scanrecparms(self, datatype, freqband, duration_tot,
                         pointing="None,None,None", integration=None):
        """Record parameters used as arguments to record_scan program."""
        self.scanrecparms = {}
        self.scanrecparms['datatype'] = datatype
        self.scanrecparms['freqband'] = freqband
        self.scanrecparms['duration_tot'] = duration_tot
        self.scanrecparms['pointing'] = pointing
        self.scanrecparms['integration'] = integration

    def write_scanrec(self, datapath):
        with open(os.path.join(datapath, self.scanrecinfo_header), "w") as f:
            f.write("# LOFAR local station project\n")
            f.write("# Created by {} version {}\n".format("iLiSA",
                                                          ilisa.__version__))
            f.write("headerversion: {}\n".format(self.headerversion))
            f.write("stnid: {}\n".format(self.stnid))
            f.write("scanrec: {!r}\n".format(self.scanrecparms))
            f.write("source: {}\n".format(self.sourcename))
            f.write("obs_ids: {!r}\n".format(list(self.obsinfos.keys())))

    def read_scanrec(self, datapath):
        try:
            _h_path = os.path.join(datapath, self.scanrecinfo_header)
            with open(_h_path, 'r') as hf:
                scanrecfiledict = yaml.safe_load(hf)
        except Exception:
            raise RuntimeError()
        self.headerversion = scanrecfiledict['headerversion']
        self.stnid = scanrecfiledict['stnid']
        self.scanrecparms = scanrecfiledict['scanrec']
        self.sourcename = scanrecfiledict.get('source')
        self.obs_ids = scanrecfiledict['obs_ids']
        try:
            scanrecfiledict['calibrationfile']
        except KeyError:
            self.calibrationfile = None
        else:
            self.calibrationfile = scanrecfiledict['calibrationfile']

    def write(self, scanrecpath=None):
        """Write scanrecinfo file and all obsinfo headers."""
        if not scanrecpath:
            scanrecpath = self.get_scanrecpath()
        self.write_scanrec(scanrecpath)
        for obs_id in self.obsinfos:
            self.obsinfos[obs_id].write_ldat_header(scanrecpath)

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
            h.write("calibrationfile: " + os.path.basename(caltabpath))

    def set_scanpath(self, scanpath):
        """Set path where this scanrec is stored. scanpath is the path to
        parent folder.
        """
        self.scanpath = scanpath

    def get_scanrecpath(self):
        """\
        Return path to this scanrec.

        Create name and destination path for folders (on the DPU) in
        which to save the various LOFAR data products.
        """

        start_key = min(self.obsinfos)
        ofi = self.scanrecparms
        ofi.update({'station_id': self.get_stnid()})
        ofi.update({'filenametime': start_key})
        ofi.update({'rcumode': self.obsinfos[start_key].rcumode})
        sb = modeparms.FreqSetup(self.scanrecparms['freqband']).subbands_spw
        ofi.update({'sb': sb})
        folder_name_beamctl_type = True
        if not folder_name_beamctl_type:
            scanrecname = ofi['filenametime']
            scanrecname += "_" + ofi['datatype']
        else:
            scanrecname = obsfileinfo2filefolder(ofi)
        scanrecpath = os.path.join(self.scanpath, scanrecname)
        return scanrecpath

    def get_datatype(self):
        return self.scanrecparms['datatype']

    def get_rcumode(self, filenr=0):
        try:
            rcumode = modeparms.FreqSetup(self.scanrecparms['freqband']
                                          ).rcumodes[0]
        except:
            try:
                rcumode = self.obsinfos[filenr].beamctl_cmd['rcumode']
            except:
                rcumode = self.scanrecparms['rcumode']
        return str(rcumode)

    def get_band(self):
        return modeparms.rcumode2band(self.get_rcumode())

    def get_bandarr(self):
        antset = modeparms.rcumode2antset_eu(self.get_rcumode())
        return antset.split('_')[0]

    def get_xcsubband(self, filenr=0):
        return int(self.obsinfos[filenr].rspctl_cmd['xcsubband'])

    def get_integration(self):
        return self.scanrecparms['integration']

    def get_pointingstr(self, filenr=0):
        return self.scanrecparms['pointing']

    def is_septon(self, filenr=0):
        obs_ids = self.get_obs_ids()
        try:
            self.obsinfos[obs_ids[filenr]]
        except:
            if self.get_datatype().endswith('SEPTON'):
                return True
            else:
                return False
        else:
            if self.obsinfos[obs_ids[filenr]].septonconf:
                return True
            else:
                return False

    def get_septon_elmap(self, filenr=0):
        obs_ids = self.get_obs_ids()
        elmap = modeparms.str2elementMap2(
                    self.obsinfos[obs_ids[filenr]].septonconf)
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

    def __init__(self, ldat_type, station_id, rcusetup_cmds, beamctl_cmds,
                 rspctl_cmds, caltabinfos=[], septonconf=None,
                 **kwargs):
        """Create observation info from parameters."""
        self.headerversion = '4'

        # ldat_type attr
        self.ldat_type = ldat_type

        # filenametime attr (Is set later, since it's know only after obs)
        self.filenametime = None

        # station_id attr
        self.station_id = station_id

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
        self.pointing = ""
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
            self.pointing = digdir

        # rspctl_cmds attr
        if rspctl_cmds == []:
            rspctl_cmds = ['rspctl']
        self.rspctl_cmds = rspctl_cmds
        rspctl_args = modeparms.parse_rspctl_args(self.rspctl_cmds)
        
        # septonconf attr
        self.septonconf = septonconf
        if self.septonconf is not None:
            self.rcumode = [5]
        
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
        
        # caltabinfos attr
        if self.ldat_type != 'bst':
            # Only need caltab info if it's BST
            caltabinfos = []
        self.caltabinfos = caltabinfos

    def isLOFARdatatype(self, obsdatatype):
        """Test if a string 'obsdatatype' is one of iLiSA's recognized LOFAR
        data types"""
        if (obsdatatype == 'acc' or
                obsdatatype == 'bst' or
                obsdatatype == 'bst-357' or
                obsdatatype == 'sst' or
                obsdatatype == 'xst' or
                obsdatatype == 'xst-SEPTON' or
                obsdatatype == 'bfs'):
            return True
        else:
            return False

    def write_ldat_header(self, datapath):
        """Create a header file for LOFAR standalone observation."""
        contents = {}
        contents['ldat_type'] = self.ldat_type
        contents['filenametime'] = self.filenametime
        contents['station_id'] = self.station_id
        contents['rcusetup_cmds'] = self.rcusetup_cmds
        contents['beamctl_cmds'] = self.beamctl_cmds
        contents['rspctl_cmds'] = self.rspctl_cmds
        if self.caltabinfos != []:
            contents['caltabinfos'] = self.caltabinfos
        if self.septonconf:
            contents['septonconf'] = self.septonconf

        if not self.isLOFARdatatype(self.ldat_type):
            raise ValueError("Unknown LOFAR statistic type {}."\
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

    def get_recfreq(self):
        """Return data recording frequency in Hz."""
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
            else:
                # headerversion == '4':
                contents = yaml.safe_load(hf)
                datatype = contents['ldat_type']
                filenametime = contents['filenametime']
                stnid = contents['station_id']
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
        obsinfo = cls(datatype, stnid, rcusetup_cmds, beamctl_cmds, rspctl_cmds,
                      caltabinfos=caltabinfos, septonconf=septonconf)
        obsinfo.filenametime = filenametime
        return obsinfo


# BEGIN BST related code
def readbstfolder(BSTfilefolder):
    """Read a BST filefolder"""
    obsfileinfo = filefolder2obsfileinfo(BSTfilefolder)
    maxnrsbs = obsfileinfo['max_nr_bls']
    BSTdirls = os.listdir(BSTfilefolder)
    BSTfiles = [f for f in BSTdirls if f.endswith('.dat')]

    # Now read the BST pol data
    BST_dtype = numpy.dtype(('f8', (maxnrsbs,)))
    BSTdata = {}
    for BSTpolfile in BSTfiles:
        pol = BSTpolfile[-5]
        with open(os.path.join(BSTfilefolder, BSTpolfile), "rb") as fin:
            BSTdata[pol] = numpy.fromfile(fin, dtype=BST_dtype)
    return BSTdata, obsfileinfo


def parse_sstfolder(SSTfolderpath):
    SSTfoldername = os.path.basename(os.path.normpath(SSTfolderpath))
    obsfolderinfo = {}
    try:
        (Ymd, HMS, rcustr, intstr, durstr, _sststr) = SSTfoldername.split('_')
        obsfolderinfo['datetime'] = datetime.datetime.strptime(Ymd + 'T' + HMS,
                                                               '%Y%m%dT%H%M%S')
        obsfolderinfo['rcumode'] = rcustr[3:]
        obsfolderinfo['integration'] = int(intstr[3:])
        obsfolderinfo['duration'] = int(durstr[3:])
    except:
        raise ValueError("Folder name not in sst_ext format.")
    return obsfolderinfo


def parse_sstfilename(SSTfilepath):
    SSTfilename = os.path.basename(SSTfilepath)
    obsfileinfo = {}
    try:
        (Ymd, HMS, _sststr, rcudatstr) = SSTfilename.split('_')
        obsfileinfo['datetime'] = datetime.datetime.strptime(Ymd + 'T' + HMS,
                                                             '%Y%m%dT%H%M%S')
        (rcu, _datext) = rcudatstr[3:].split('.')
        obsfileinfo['rcu'] = int(rcu)
    except:
        raise ValueError("File name not in sst format.")
    return obsfileinfo


def readsst(SSTfile):
    """Read-in SST datafile.

    Parameters
    ----------
    SSTfile : str
        Name of SST datafile.

    Returns
    -------
    SSTdata : (512, N)
        The SST data, where N is the number of time samples.
    """
    obsfileinfo = parse_sstfilename(SSTfile)
    # Now read the SST data
    SST_dtype = numpy.dtype(('f8', (512,)))
    with open(SSTfile, "rb") as fin:
        SSTdata = numpy.fromfile(fin, dtype=SST_dtype)
    return SSTdata, obsfileinfo


def readsstfolder(SSTfolder):
    """Read-in SST datafile.

    Parameters
    ----------
    SSTfolder : str
        The name of the folder which contains an SST datafile for each RCU.

    Returns
    -------
    SSTdatarcu : (192, 512, N)
        The SST data, where N is the number of time samples.
    """
    obsfolderinfo = filefolder2obsfileinfo(SSTfolder)
    #obsfolderinfo = parse_sstfolder(SSTfolder)
    files = os.listdir(SSTfolder)
    SSTfiles = [f for f in files if f.endswith('.dat')]
    SSTfiles.sort()
    SSTdatarcu = [None] * len(SSTfiles)
    for sstfile in SSTfiles:
        SSTrcudata, obsfileinfo = readsst(os.path.join(SSTfolder, sstfile))
        SSTdatarcu[obsfileinfo['rcu']] = SSTrcudata
    return SSTdatarcu, obsfolderinfo


class CVCfiles(object):
    """Provides functionality for covariance cube (CVC) files. (CVC is
    essentially visiblity cubes.) CVC files from a LOFAR station includes ACC
    and XST files.

    Attributes
    ----------
    dataset: list of array_like
        Each item in list corresponds to one CVC file and is the actual
        covariance matrix cube with shape cvcdim0 x cvcdim1 x cvcdim2.
    samptimeset: list of datetimes
        The datetime of the visibility matrix sample.
    """
    NRRCUS_EU = 192  # Default number of RCUs on EU stations

    def __init__(self, datapath):
        self._getitem_prev_file_nr_data_ = (None, None)
        self.samptimeset = []
        self.freqset = []

        self.cvcdim1 = self.NRRCUS_EU
        self.cvcdim2 = self.NRRCUS_EU

        datapath = os.path.abspath(datapath)
        if os.path.isdir(datapath):
            obsfileinfo = filefolder2obsfileinfo(datapath)
            stnid = obsfileinfo['station_id']
            nrrcus = modeparms.nrrcus_stnid(stnid)
            self.cvcdim1 = nrrcus
            self.cvcdim2 = nrrcus
            self.filefolder = datapath
            (self.scanrecinfo, self.filenames, self.samptimeset, self.freqset
             ) = self._readcvcfolder()
        elif os.path.isfile(datapath):
            # FIXME:
            self._readcvcfile(datapath)
        else:
            raise ValueError('Path does not exist')
        # Get/Compute ant positions
        stnid = self.scanrecinfo.get_stnid()
        antset = next(iter(self.scanrecinfo.obsinfos.values())).antset
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

    def _parse_cvcfolder(self, cvcfolderpath):
        """Parse the cvc filefolder.

        The filefolder should have the format:
            stnid_Ymd_HMS_rcumode_subband_integration_duration_pointing_cvc

        Parameters
        ----------
        cvcfolderpath: str
            Path of folder with CVC files.

        Returns
        -------
        scanrecparms
        """
        cvcfoldername = os.path.basename(os.path.abspath(cvcfolderpath))
        obsfolderinfo = {}
        cvcextstr = cvcfoldername.split('_')[-1]
        if cvcextstr == 'xst' or cvcextstr == 'xst-SEPTON':
            cvcfoldername_split = cvcfoldername.split('_')
            try:
                (stnid, Ymd, HMS, rcustr, sbstr, intstr, durstr, dirstr, cvcextstr
                 ) = cvcfoldername_split
                obsfolderinfo['stnid'] = stnid
                obsfolderinfo['datetime'] = datetime.datetime.strptime(
                    Ymd + 'T' + HMS, '%Y%m%dT%H%M%S')
                obsfolderinfo['rcumode'] = rcustr[3:]
                obsfolderinfo['subband'] = sbstr[2:]
                obsfolderinfo['integration'] = float(intstr[3:])
                obsfolderinfo['duration_tot'] = float(durstr[3:])
                obsfolderinfo['pointing'] = dirstr[3:].split(',')
            except:
                raise ValueError("Foldername not in xst_ext format.")
        elif cvcextstr == 'acc':
            dirpat = re.compile(regex_ACCfolder)
            obsdirinfo_m = dirpat.match(cvcfoldername)
            if obsdirinfo_m is None:
                print("Cal error")
                raise ValueError(
                        "Calibration directory does not have correct syntax.")
            obsdirinfo = obsdirinfo_m.groupdict()
            obsfolderinfo['stnid'] = obsdirinfo['stnid']
            d0 = datetime.datetime(int(obsdirinfo['year']),
                                   int(obsdirinfo['month']),
                                   int(obsdirinfo['day']),
                                   int(obsdirinfo['hour']),
                                   int(obsdirinfo['minute']),
                                   int(obsdirinfo['second']))
            obsfolderinfo['sessiontimeid'] = d0
            obsfolderinfo['rcumode'] = obsdirinfo['rcumode']
            obsfolderinfo['subband'] = '0:511'
            obsfolderinfo['integration'] = 1.0
            obsfolderinfo['duration_tot'] = int(obsdirinfo['duration_tot'])
            obsfolderinfo['source'] = obsdirinfo['calsrc']
            obsfolderinfo['pointing'] = \
                ilisa.monitorcontrol.directions.std_pointings(
                    obsfolderinfo['source'])
        else:
            raise ValueError("Folder not expected xst or acc format.")
        obsfolderinfo['datatype'] = cvcextstr
        return obsfolderinfo

    def _readcvcfolder(self):
        """Read in CVC data from the filefolder.

        The filefolder name may have the format as specified in the
        parse_cvcfolder() method. The contents of the data file is stored in
        the object attribute:
           data : [(N,192,192), ... , (N,192,192)]
        where N is nominally the number of time samples and the len of data is
        the number of files in the folder.
        """
        # Initialize
        scanrecinfo = ScanRecInfo()
        samptimeset = []
        freqset = []
        try:
            scanrecinfo.read_scanrec(self.filefolder)
        except Exception:
            warnings.warn("Could not read session header."
                          +" Will try filefolder name...")
            try:
                obsfolderinfo = self._parse_cvcfolder(self.filefolder)
            except ValueError as er:
                print(er)
                scanrecinfo.scanrecparms = None
            else:
                spw = obsfolderinfo['rcumode']
                nqz = modeparms.rcumode2nyquistzone(spw)
                sbs = modeparms.seqarg2list(obsfolderinfo['subband'])
                freqspec_hi = modeparms.sb2freq(sbs[-1], nqz)
                scanrecinfo.set_scanrecparms(obsfolderinfo['datatype'],
                                             str(freqspec_hi),
                                             obsfolderinfo['duration_tot'],
                                             obsfolderinfo['pointing'],
                                             obsfolderinfo['integration'])
                scanrecinfo.scanrecparms['rcumode'] = spw
                scanrecinfo.set_stnid(obsfolderinfo['stnid'])
                scanrecinfo.calibrationfile = None
                print("Read in filefolder meta.")
        # Select only data files in folder (avoid CalTable*.dat files)
        ls = os.listdir(self.filefolder)
        filenames = [filename for filename in ls if filename.endswith('.dat')
                     and not filename.startswith('CalTable')]
        filenames.sort()  # This enforces chronological order
        for cvcfile in filenames:
            cvcdim_t = (os.path.getsize(os.path.join(self.filefolder, cvcfile))
                        // self.__get_cvc_dtype().itemsize)
            # Try to get obsfile header
            try:
                (bfilename, _dat) = cvcfile.split('.')
                hfilename = bfilename + '.h'
                # Check if xst might have some extra stuff in name
                #   So first get ldattype
                ymd, hms, ldattype_full = bfilename.split('_', 2)
                if '_' in ldattype_full:
                    ldattype, _rest = ldattype_full.split('_', 1)
                else:
                    ldattype = ldattype_full
                if ldattype == 'xst':
                    hfilename = ymd+'_'+hms+'_'+ldattype+'.h'
                hfilepath = os.path.join(self.filefolder, hfilename)
                obsinfo = LDatInfo.read_ldat_header(hfilepath)
                scanrecinfo.add_obs(obsinfo)
            except:
                warnings.warn(
                    "Couldn't find a header file for {}".format(cvcfile))
            _datatype, t_begin = self._parse_cvcfile(os.path.join(self.filefolder, cvcfile))

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
                sb = obsinfo.sb
                freq = modeparms.sb2freq(sb, nz)
                freqs = [freq] * cvcdim_t
            freqset.append(freqs)
        return scanrecinfo, filenames, samptimeset, freqset

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


def cov_flat2polidx(cvc, parity_ord=True):
    """
    Convert flat array covariance matrix (visibilities) to polarization
    component indexed array covariance matrix.

    Parameters
    ----------
    cvc : array_like
        Covariance cube to be converted. Shape should be (..., 2*N, 2*N),
        where N is the number of dual-polarized elements.
    parity_ord : Boolean
        If True (default) then the polarization component is determined from
        the index's parity: even index maps to component 0, odd to 1.
        If False the baseline indices are split into a first and second half
        and mapped to pol component 0,1 respectively.

    Returns
    -------
    cvcpol : array_like
        Polarization index array covariance cube. Shape will be 
        (2, 2, ..., N, N). Polarization component order from flat
        visibility will 

    Notes
    -----
    This function is agnostic to whether the polarization basis is
    linear or circular. Also the ordering is conserved from the
    flat ordering, so the mapping of index 0,1 to say L,R or X,Y
    components is determined from the original ordering.

    Examples
    --------
    Minimal example:
    >>> cov_flat2polidx(numpy.arange(4*4).reshape((4,4))
    (2.0, 0.0, 2.0, 2.0)
    """
    if parity_ord:
        pp = cvc[..., ::2, ::2]
        qq = cvc[..., 1::2, 1::2]
        pq = cvc[..., ::2, 1::2]
        qp = cvc[..., 1::2, ::2]
    else:
        # First-half, second-half order
        n =  cvc.shape[-1]/2
        pp = cvc[..., :n, :n]
        qq = cvc[..., n:, n:]
        pq = cvc[..., :n, n:]
        qp = cvc[..., :n, n:]
    cvpol = numpy.array([[pp, pq], [qp, qq]])
    return cvpol


def cvc2polrep(cvc, crlpolrep='lin'):
    """
    Convert a flat-indexed polarized covariance cube `cvc` into an array
    indexed according to correlated polarization representation
    `crlpolrep`.

    Parameters
    ----------
    cvc : (T,2*N,2*N) array (usually T=512 and 2*N=196)
        The Covariance Cube array produced by an International LOFAR
        station when it is in calibration mode. It is the covariance
        matrices of the 196 rcus (98 X-polarized & 98 Y-polarized 
        interleaved) over 512subbands.
    crlpolrep : str
        Correlated polarization representation can be 'lin', 'cir', 'sto' where:
            'lin' is linear pol rep
            'cir' is circular pol rep
            'sto' is Stokes pol rep.
    
    Returns
    -------
    cvpol : array
        Polarized visibilities. Shape depends on `crlpolrep`:
            * (2,2,T,N,N) for 'lin' or 'cir'. Basis is indexed with
              0,1 corresponding to X,Y or L,R resp.
            * (4,T,N,N) for 'sto'. Basis is indexed with
              0,1,2,3 corresponding to Stokes I,Q,U,V resp.
    
    Notes
    -----
    Requires dreamBeam.
    """
    cvcpolidx = cov_flat2polidx(cvc)
    if crlpolrep == 'lin':
        cvpol = cvcpolidx
    elif crlpolrep == 'cir':
        cvpol = cov_lin2cir(cvcpolidx)
    elif crlpolrep == 'sto':
        print("Hej")
        (S0, S1, S2, S3) =\
            convertxy2stokes(cvcpolidx[0][0], cvcpolidx[0][1], cvcpolidx[1][0],
                             cvcpolidx[1][1])
        cvpol = numpy.array([S0, S1, S2, S3])
    else:
        raise ValueError("No such correlated pol rep: '{}'".format(crlpolrep))
    return cvpol


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
    pntstr = ilisa.monitorcontrol.directions.normalizebeamctldir(calsrc)
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


import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.dates as mdates
import ilisa.monitorcontrol.modeparms as modeparms


def viewbst(bstff, pol_stokes=True, printout=False):
    """Plot BST data."""
    BSTdata, obsfileinfo = readbstfolder(bstff)
    stnid = obsfileinfo['station_id']
    starttime = obsfileinfo['datetime']
    intg = obsfileinfo['integration']
    dur = obsfileinfo['duration_scan']
    freqs = obsfileinfo['frequencies']
    pointing = obsfileinfo['pointing']

    ts = numpy.arange(0., dur, intg)
    ts = [starttime+datetime.timedelta(seconds=t) for t in ts]
    fig, (ax_p, ax_q) = plt.subplots(2, 1, sharex=True, sharey=True)
    data2plot_p_name, data2plot_p = 'X-pol', BSTdata['X']
    data2plot_q_name, data2plot_q = 'Y-pol', BSTdata['Y']
    data2plot_q_unit = 'Flux [arb. units]'
    norm_p, norm_q = None, None
    cmap_q = None
    if pol_stokes:
        data2plot_p_name = 'Stokes I'
        data2plot_p = BSTdata['X'] + BSTdata['Y']
        norm_p = colors.LogNorm()
        data2plot_q_name = '(antenna) Stokes Q'
        data2plot_q = BSTdata['X'] - BSTdata['Y']
        data2plot_q_unit = 'Signed flux [arb. units]'
        cmap_q = 'RdBu_r'
        norm_q = colors.SymLogNorm(linthresh=1e2)
        stokes_norm = True
        if stokes_norm:
            data2plot_q = data2plot_q / data2plot_p
            data2plot_q_name = '(antenna) Stokes q'
            data2plot_q_unit = 'Signed relative flux [%]'
            norm_q = colors.SymLogNorm(linthresh=1e-3)
    if not printout:
        # Plot quantity p:
        bstplt_p = ax_p.pcolormesh(ts, freqs/1e6, data2plot_p.T, norm=norm_p)
        cbar_p = fig.colorbar(bstplt_p, ax=ax_p)
        cbar_p.set_label('Flux [arb. units]')
        ax_p.set_ylabel('Frequency [MHz]')
        ax_p.set_title('{}'.format(data2plot_p_name))
        # Plot quantity q:
        bstplt_q = ax_q.pcolormesh(ts, freqs / 1e6, data2plot_q.T, cmap=cmap_q,
                                   norm=norm_q)
        cbar_q = fig.colorbar(bstplt_q, ax=ax_q)
        cbar_q.set_label(data2plot_q_unit)
        ax_q.set_title('{}'.format(data2plot_q_name))
        fig.autofmt_xdate()

        ax_q.xaxis.set_major_formatter( mdates.DateFormatter('%H:%M:%S'))
        ax_q.set_xlabel('Datetime [UT]  Starts: {}'.format(starttime))
        ax_q.set_ylabel('Frequency [MHz]')

        supertitle = ('{} BST intg: {}s dur: {}s'.format(stnid, intg, dur)
                      + ' pointing: {},{},{}'.format(*pointing))
        plt.suptitle(supertitle)
        plt.show()
    else:
        # CSV style:
        #   Header
        t_prev = ts[0]
        print("#H:M:S since {} UT".format(t_prev.isoformat()), "Freq[MHz]",
              data2plot_p_name, data2plot_q_name, sep=', ')
        #   Data
        for ti, t in enumerate(ts):
            for freqi, freq in enumerate(freqs):
                dataval_p, dataval_q= data2plot_p[ti, freqi], data2plot_q[ti, freqi]
                del_t = t - t_prev
                print(del_t, freq/1e6, dataval_p, dataval_q, sep=', ')


def plotsst(sstff, freqreq):
    """Plot SST data."""
    SSTdata, obsfolderinfo = readsstfolder(sstff)
    scanrecinfo = ScanRecInfo()
    scanrecinfo.read_scanrec(sstff)
    starttime = obsfolderinfo['datetime']
    SSTdata = numpy.array(SSTdata)
    freqs = obsfolderinfo['frequencies']
    sbreq = None
    if freqreq:
        sbreq = int(numpy.argmin(numpy.abs(freqs-freqreq)))
        show = 'persb'
    else:
        show = 'mean'
    intg = obsfolderinfo['integration']
    dur = obsfolderinfo['duration_scan']
    ts = [starttime + datetime.timedelta(seconds=td) for td in
          numpy.arange(0., dur, intg)]
    if show == 'mean':
        meandynspec = numpy.mean(SSTdata, axis=0)
        res = meandynspec
        if res.shape[0] > 1:
            plt.pcolormesh(freqs/1e6, ts, res, norm=colors.LogNorm())
            plt.colorbar()
            plt.title('Mean (over RCUs) dynamicspectrum\n'
                      + 'Starttime: {} Station: {}'
                      .format(starttime, scanrecinfo.stnid))
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
        res = SSTdata[:, :, sbreq]
        resX = res[0::2, :]
        resY = res[1::2, :]
        plt.subplot(211)
        if ampVStime:
            plt.plot(ts, numpy.transpose(resX))
            # plt.gcf().autofmt_xdate()
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
    plt.show()


def plotxst(xstff):
    """Plot XST data."""
    colorscale = None  # Colorscale for xst data plot (default None)
    if colorscale == 'log':
        normcolor = colors.LogNorm()
    else:
        normcolor = None
    xstobj = CVCfiles(xstff)
    obs_ids = xstobj.scanrecinfo.get_obs_ids()
    # Assume freq sweep over nr of files, so filenr is also sb.
    for sbstepidx in range(xstobj.getnrfiles()):
        obsinfo = xstobj.scanrecinfo.obsinfos[obs_ids[sbstepidx]]
        intg = obsinfo.integration
        dur = obsinfo.duration_subscan

        freq = obsinfo.get_recfreq()
        ts = numpy.arange(0., dur, intg)
        xstdata = xstobj[sbstepidx]
        for tidx in range(xstdata.shape[0]):
            print("Kill plot window for next plot...")
            plt.imshow(numpy.abs(xstdata[tidx,...]), norm=normcolor,
                       interpolation='none')
            plt.title("""Time (from start {}) {}s
                      @ freq={} MHz""".format(obsinfo.get_starttime(),
                      ts[tidx], freq/1e6))
            plt.xlabel('RCU [#]')
            plt.ylabel('RCU [#]')
            plt.colorbar()
            plt.show()


def plotacc(accff, freqreq=None):
    """Plot of ACC folder files."""
    dataobj = CVCfiles(accff)
    if freqreq is None:
        freqreq = 0.0
    sb, _nqzone = modeparms.freq2sb(freqreq)
    for fileidx in range(0, dataobj.getnrfiles()):
        filecvc = dataobj[fileidx]
        while sb < 512:
            fig, (ax1, ax2) = plt.subplots(2, 1, sharex=True, sharey=True)
            absdatplt = ax1.pcolormesh(numpy.abs(filecvc[sb]))
            ax1.set_title('Abs value')
            ax1.set_ylabel('RCU [#]')
            fig.colorbar(absdatplt, ax=ax1)
            angdatplt = ax2.pcolormesh(numpy.angle(filecvc[sb]),
                                       cmap=plt.get_cmap('hsv'))
            ax2.set_title('Phase value')
            ax2.set_xlabel('RCU [#]')
            ax2.set_ylabel('RCU [#]')
            fig.colorbar(angdatplt, ax=ax2)
            plt.suptitle('Station element covariance. Time: {}UT, SB: {}'\
                                .format(dataobj.samptimeset[fileidx][sb], sb))
            plt.show()
            sb += 1


def view_bsxst(args):
    lofar_datatype = datafolder_type(args.dataff)
    if lofar_datatype =='acc':
        plotacc(args.dataff, args.freq)
    if lofar_datatype=='bst' or lofar_datatype=='bst-357':
        viewbst(args.dataff, pol_stokes=not(args.linear),
                printout=args.printout)
    elif lofar_datatype=='sst':
        plotsst(args.dataff, args.freq)
    elif lofar_datatype=='xst' or lofar_datatype=='xst-SEPTON':
        plotxst(args.dataff)
    else:
        raise RuntimeError("Not a bst, sst, or xst filefolder")


import argparse


def main():
    parser = argparse.ArgumentParser()

    parser.add_argument('-n', '--filenr', type=int, default=0)
    parser.add_argument('-s', '--sampnr', type=int, default=0)
    parser.add_argument('-l', '--linear', action="store_true",
                        help='Use linear X,Y polarization rather than Stokes')
    parser.add_argument('-p', '--printout', action="store_true",
                        help='Print out data to stdout, else plot')
    parser.add_argument('dataff', help="acc, bst, sst or xst filefolder")
    parser.add_argument('freq', type=float, nargs='?', default=None)

    args = parser.parse_args()
    args.dataff = os.path.normpath(args.dataff)
    view_bsxst(args)


if __name__ == "__main__":
    main()
