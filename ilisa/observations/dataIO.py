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
import ilisa
import ilisa.observations.directions
import ilisa.observations.modeparms as modeparms

regex_ACCfolder=(
"^(?P<stnid>\w{5})_(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})"
"_(?P<hour>\d{2})(?P<minute>\d{2})(?P<second>\d{2})"
"_rcu(?P<rcumode>\d+)_dur(?P<duration_tot>\d+)(_(?P<calsrc>\w+))?_acc$")
regex_ACCfilename=(
"^(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})"
"_(?P<hour>\d{2})(?P<minute>\d{2})(?P<second>\d{2})"
"_acc_(?P<totnrsb>\d+)x(?P<nrrcu0>\d+)x(?P<nrrcu1>\d+)"
".dat$")
regex_xstfilename=(
"^(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})"
"_(?P<hour>\d{2})(?P<minute>\d{2})(?P<second>\d{2})"
"_rcu(?P<rcumode>\d+)_int(?P<integration>\d+)_dur(?P<duration_scan>\d+)"
"_dir(?P<RAint>\d+).(?P<RAdecimal>\d+),(?P<DECint>\d+).(?P<DECdecimal>\d+),(?P<ref>\s+)"
"_xst.dat$")


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


class ScanRecInfo(object):
    """This class maintains info on a scan recording (scanrec), which is one of the
    results of an iLiSA scan. One scanrec is a group of one or more files of a unique
    LOFAR station data product (ldat), i.e. acc, bfs, bst, sst or xst.
    The info in a ScanRecInfo object consists of sufficient parameters to redo the
    data, namely the stn_id, the iLiSA scanrec parameters, and a list of ObsInfo objects
    called obsinfos that maps to each ldat within the scanrec.
    """
    scanrecinfo_header = "SCANREC_INFO.yml"

    def __init__(self):
        self.headerversion = 2
        self.obsinfos = {}

    def add_obs(self, obsinfo):
        """Add an ObsInfo object to this ScanRecInfo."""
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

    def set_scanrecparms(self, datatype, freqband, duration_tot,
                         pointing="None,None,None", integration=None, allsky=False):
        """Record parameters used as arguments to record_scan program."""
        self.scanrecparms = {}
        self.scanrecparms['datatype'] = datatype
        self.scanrecparms['freqband'] =  freqband
        self.scanrecparms['duration_tot'] = duration_tot
        self.scanrecparms['pointing'] = pointing
        self.scanrecparms['integration'] = integration
        self.scanrecparms['allsky'] = allsky

    def write_scanrec(self, datapath):
        with open(os.path.join(datapath, self.scanrecinfo_header), "w") as f:
            f.write("# LOFAR local station project\n")
            f.write("# Created by {} version {}\n".format("iLiSA", ilisa.__version__))
            f.write("headerversion: {}\n".format(self.headerversion))
            f.write("stnid: {}\n".format(self.stnid))
            f.write("scanrec: {!r}\n".format(self.scanrecparms))
            f.write("obs_ids: {!r}\n".format(self.obsinfos.keys()))

    def read_scanrec(self, datapath):
        with open(os.path.join(datapath, self.scanrecinfo_header), 'r') as hf:
            try:
                scanrecfiledict = yaml.load(hf)
            except Exception as e:
                print("Couldn't load yaml formatted scan header file.")
        self.headerversion = scanrecfiledict['headerversion']
        self.stnid = scanrecfiledict['stnid']
        self.scanrecparms = scanrecfiledict['scanrec']
        self.obs_ids = scanrecfiledict['obs_ids']
        try:
            scanrecfiledict['calibrationfile']
        except KeyError:
            self.calibrationfile = None
        else:
            self.calibrationfile = scanrecfiledict['calibrationfile']

    def write(self, scanrecpath):
        """Write scanrecinfo file and all obsinfo headers."""
        self.write_scanrec(scanrecpath)
        for obs_id in self.obsinfos:
            self.obsinfos[obs_id].write_ldat_header(scanrecpath)

    def set_postcalibration(self, caltabpath, scanrecpath):
        """Add the caltab file that was applied to this scanrec (typically CVC data)
         after station recording.

        Parameters
        ----------
        caltabpath: str
            Path to caltab file
        scanrecpath: str
            Path to scanrec folder
        """
        shutil.copy(caltabpath, scanrecpath)
        scanrecinfo_header_path  = os.path.join(scanrecpath, self.scanrecinfo_header)
        with open(scanrecinfo_header_path, 'a') as h:
            h.write("calibrationfile: "+os.path.basename(caltabpath))

    def scanrecfolder(self):
        start_key = min(self.obsinfos)
        return self.obsinfos[start_key].obsfoldername()

    def get_datatype(self):
        return self.scanrecparms['datatype']

    def get_rcumode(self,filenr=0):
        try:
            rcumode = modeparms.FrequencyBand(self.scanrecparms['freqband']).rcumodes[0]
        except:
            try:
                rcumode = self.obsinfos[filenr].beamctl_cmd['rcumode']
            except:
                rcumode = self.scanrecparms['rcumode']
        return str(rcumode)

    def get_band(self):
        return modeparms.rcumode2band(self.get_rcumode())

    def get_bandarr(self):
        antset = modeparms.rcumode2antset(self.get_rcumode())
        return antset.split('_')[0]

    def get_xcsubband(self, filenr=0):
        return int(self.obsinfos[filenr].rspctl_cmd['xcsubband'])

    def get_integration(self):

        return self.scanrecparms['integration']

    def get_pointingstr(self, filenr=0):
        return self.scanrecparms['pointing']

    def is_septon(self, filenr=0):
        try:
            self.obsinfos[filenr]
        except:
            if self.get_datatype().endswith('SEPTON'):
                return True
            else:
                return False
        else:
            if self.obsinfos[filenr].septonconf:
                return True
            else:
                return False

    def get_septon_elmap(self, filenr=0):
        elmap = modeparms.str2elementMap2(self.obsinfos[filenr].septonconf)
        return elmap

    def get_ldat_filenames(self):
        datatype = self.get_datatype()
        return ["{}_{}.dat".format(obs_id, datatype) for obs_id in self.obs_ids]


class ObsInfo(object):
    """Maintains a minimum of info for the physical interpretation of an LCU data product
    file (ldat). It consists of the LCU commands or settings that resulted in the ldat.
    An observation can update an object with LCU settings and the create a header file
    associated with the resulting ldat file of the observation. Later when ingesting
    the ldat file, an object can read-in the associated header file to get physically
    meaningful data.

    Specifically an object should have:

    filename:

    """
    def __init__(self, lofardatatype, filenametime, beamctl_cmd, rspctl_cmd,
                 caltabinfos="", septonconf=""):
        """Create observation info from parameters."""
        self.ldat_type = lofardatatype

        self.filenametime = filenametime
        self.beamctl_cmd = beamctl_cmd
        if self.beamctl_cmd != "" and self.beamctl_cmd is not None:
            # FIXME better support for multiline beamctl cmds.
            if type(self.beamctl_cmd) is list:
                self.rcumode = []
                self.sb = []
                self.bl = []
                for beamctl_cmd in self.beamctl_cmd:
                    (antset, rcus, rcumode, beamlets, subbands, anadir, digdir
                     ) = modeparms.parse_beamctl_args(beamctl_cmd)
                    self.rcumode.append(int(rcumode))
                    self.sb.append(subbands)
                    self.bl.append(beamlets)
            else:
                (antset, rcus, rcumode, beamlets, subbands, anadir, digdir
                ) = modeparms.parse_beamctl_args(beamctl_cmd)
                self.rcumode = [rcumode]
                self.sb = [subbands]
            self.pointing = digdir
        else:
            self.pointing = ""
        if rspctl_cmd == '':
            rspctl_cmd = 'rspctl'
        self.rspctl_cmd = rspctl_cmd
        rspctl_args = modeparms.parse_rspctl_args(self.rspctl_cmd)
        if self.ldat_type != 'bfs' and self.ldat_type != 'acc':
            self.integration = float(rspctl_args['integration'])
            self.duration_scan = float(rspctl_args['duration'])
        if self.ldat_type == 'sst':
            self.sb = ""
        elif self.ldat_type.startswith('xst'):
            self.sb = str(rspctl_args['xcsubband'])
        elif self.ldat_type == 'bst':
            self.sb = self.sb
        self.caltabinfos = caltabinfos
        self.septonconf = septonconf
        if self.septonconf != "":
            self.rcumode = 5

    def obsfoldername(self, folder_name_beamctl_type = True):
        """Create name and destination path for folders (on the DPU) in
        which to save the various LOFAR data products.
        """
        obsextname = self.filenametime
        if folder_name_beamctl_type:
            if self.ldat_type == "bst-357":
                obsextname += "_rcu357"
            else:
                if type(self.rcumode) is list:
                    rcumodestr = ''.join([str(rcumode) for rcumode in self.rcumode])
                else:
                    rcumodestr = str(self.rcumode)
                obsextname += "_rcu"+rcumodestr
            if self.sb != [] and self.sb != '':
                obsextname += "_sb"
                obsextname += modeparms.seqlists2slicestr(self.sb)
            if hasattr(self, 'integration'):
                obsextname += "_int"+str(int(self.integration))
            if hasattr(self, 'duration_scan'):
                obsextname += "_dur"+str(int(self.duration_scan))
            if self.ldat_type != 'sst':
                if str(self.pointing) != "":
                    obsextname += "_dir"+str(self.pointing)
                else:
                    obsextname += "_dir,,"
        obsextname += "_"+self.ldat_type
        return obsextname


    def isLOFARdatatype(self, obsdatatype):
        """Test if a string 'obsdatatype' is one of iLiSA's recognized LOFAR data types"""
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
        ldat_type = self.ldat_type
        YMD_hms = self.filenametime
        if type(self.beamctl_cmd) is not list:
            beamctl_CMD = [self.beamctl_cmd]
        else:
            beamctl_CMD = self.beamctl_cmd
        beamctl_CMD = '\n'.join(beamctl_CMD)
        rspctl_CMD = self.rspctl_cmd
        caltabinfos = self.caltabinfos
        septonconfig = self.septonconf
        def indenttext(txt):
            indentstr = "  "
            return indentstr+txt.replace("\n","\n"+indentstr)
        headerversion = "3"
        if not self.isLOFARdatatype(ldat_type):
            raise ValueError("Unknown LOFAR statistic type {}.\
                              ".format(ldat_type))
        xtra = ''
        if ldat_type == 'acc':
            xtra = '_512x192x192'
        ldat_header_filename = YMD_hms+"_"+ldat_type+xtra+".h"
        f = open(os.path.join(datapath, ldat_header_filename), "w")
        f.write("# HeaderType: bsxSTdata (YAML)\n")
        f.write("# Header version {}\n".format(headerversion))
        f.write("datatype: {}\n".format(ldat_type))
        filetime = YMD_hms[0:4]+'-'+YMD_hms[4:6]+'-'\
                        + YMD_hms[6:8]+'T'+YMD_hms[9:11]+':'\
                        + YMD_hms[11:13]+':'+YMD_hms[13:15]
        f.write("filetime: "+filetime+"\n")
        if septonconfig is not "":
            f.write("SEPTONconfig: {}\n".format(septonconfig))
        f.write("beamctl_cmds: |-\n")
        f.write(indenttext(beamctl_CMD)+"\n")
        # f.write(rspsetup_CMD+"\n")
        # FIX separation of beamctl and rspsetup
        # (Currently rspsetup is in beamctl)
        if rspctl_CMD != "":
            f.write("rspctl_cmds: |-\n")
            f.write(indenttext(rspctl_CMD)+"\n")
        if ldat_type == 'bst' or ldat_type == 'bfs':
            f.write("caltabinfos:\n")
            #f.write(indenttext(caltableInfo))
            for caltabinfo in caltabinfos:
                f.write("  - ")
                f.write(str(caltabinfo))
                f.write("\n")
        f.close()

    def get_recfreq(self):
        """Return data recording frequency in Hz."""
        sb = self.sb
        if self.ldat_type != "xst-SEPTON" and  not self.septonconf:
            rcumode = self.rcumode
        else:
            rcumode = 5
        nz = modeparms.rcumode2nyquistzone(rcumode)
        return modeparms.sb2freq(sb, nz)

    def get_starttime(self):
        """Return the datetime when this obs started."""
        filetime = datetime.datetime.strptime(self.filenametime, "%Y%m%d_%H%M%S")
        if self.ldat_type != 'acc':
            starttime = filetime
        else:
            starttime =  filetime - datetime.timedelta(seconds=512)
        return starttime

    @classmethod
    def read_ldat_header(cls, headerpath):
        """Parse a ldat header file and return it as an ObsInfo."""
        # TODO extract CalTable info.
        if os.path.isdir(headerpath):
            files = os.listdir(headerpath)
            headerfiles = [f for f in files if f.endswith('.h')]
            headerfile = os.path.join(headerpath,headerfiles.pop())
        else:
            headerfile = headerpath
        stnid = None
        starttime = None
        with open(headerfile,'r') as hf:
            for hline in hf:
                if "Header version" in hline:
                    headerversion = hline.split()[-1]
        with open(headerfile, 'r') as hf:
            if headerversion == '1':
                rspctl_lines = []
                for line in hf:
                    if "Observer" in line:
                        label, observer = line.split('=')
                    if "Project" in line:
                        label, project = line.split('=')
                    if "DataType" in line:
                        label, datatype = line.split('=')
                    if "StationID" in line:
                        label, stnid = line.split('=')
                        stnid = stnid.strip()
                    if "StartTime" in line:
                        label, starttime = line.split('=')
                        starttime = starttime.strip()
                    if "beamctl" in line:
                        # HACK
                        beamctl_line = line
                    if "rspctl" in line:
                        rspctl_lines.append(line)
            elif headerversion == '2':
                contents = yaml.load(hf)
                observer = contents['Observer']
                project = contents['Project']
                datatype = contents['DataType']
                stnid = contents['StationID']
                starttime = contents['StartTime']
                beamctl_line = contents['BeamctlCmds']
                rspctl_lines = contents['RspctlCmds'].split('\n')
            else:
                contents = yaml.load(hf)
                datatype = contents['datatype']
                starttime = contents['filetime']
                beamctl_line = contents['beamctl_cmds']
                rspctl_lines_raw = contents['rspctl_cmds']
                rspctl_lines  = rspctl_lines_raw.split('\n')
        multishellcmds = beamctl_line.split('&')
        beamctl_cmd = multishellcmds[0]
        if beamctl_cmd is not "":
            (antennaset, rcus, rcumode, beamlets, subbands, anadir, digdir) \
             = modeparms.parse_beamctl_args(beamctl_cmd)
            beamctl_cmd = {'antennaset': antennaset,
                       'rcus': rcus,
                       'rcumode': rcumode,
                       'beamlets': beamlets,
                       'subbands': subbands,
                       'anadir': anadir,
                       'digdir': digdir}
            septonconf = None
        elif 'SEPTONconfig' in contents:
            beamctl_cmd = ""
            septonconf = contents['SEPTONconfig']
        rspctl_cmd = {}
        if rspctl_lines is not "":
            for rspctl_line in rspctl_lines:
                rspctl_args = modeparms.parse_rspctl_args(rspctl_line)
                rspctl_cmd.update(rspctl_args)

        filenametime = starttime.strftime('%Y%m%d_%H%M%S')
        obsinfo = cls(datatype, filenametime, beamctl_line, rspctl_lines_raw,
                      caltabinfos="", septonconf=septonconf)
        return obsinfo


# BEGIN BST related code
def parse_bstfolder(BSTfilepath):
    BSTfilepath = os.path.normpath(BSTfilepath)
    BSTfilename = os.path.basename(BSTfilepath)
    obsfileinfo = {}
    try:
        (Ymd, HMS, rcustr, sbstr, intstr, durstr, dirstr, bststr
        ) = BSTfilename.split('_')
        obsfileinfo['datetime'] = datetime.datetime.strptime(Ymd+'T'+HMS,'%Y%m%dT%H%M%S')
        obsfileinfo['rcumode'] =     rcustr[3:]
        obsfileinfo['subbands'] =    sbstr[2:]
        obsfileinfo['integration'] = int(intstr[3:])
        obsfileinfo['duration'] =    int(durstr[3:])
        obsfileinfo['pointing'] =    dirstr[3:].split(',')
    except:
        raise ValueError("Filename not in bst_ext format.")
    if len(obsfileinfo['rcumode']) > 1:
        obsfileinfo['rcumode'] = list(obsfileinfo['rcumode'])
    if modeparms.rcusbsep in obsfileinfo['subbands']:
        obsfileinfo['subbands'] = obsfileinfo['subbands'].split(modeparms.rcusbsep)
    return obsfileinfo


def readbstfolder(BSTfilefolder):
    obsfileinfo = parse_bstfolder(BSTfilefolder)
    if type(obsfileinfo['rcumode']) is not list:
        obsfileinfo['rcumode'] = [obsfileinfo['rcumode']]
    if type(obsfileinfo['subbands']) is not list:
        obsfileinfo['subbands'] = [obsfileinfo['subbands']]
    obsfileinfo['frequencies'] = numpy.empty(0)
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
                                                  numpy.linspace(freqlo, freqhi, nrsbs))
        totnrsbs += nrsbs

    # When the beamlets allocated is less than the maximum (given by bit depth) the
    # RSPs fill the remaining ones regardless. Hence we have to account for them:
    if totnrsbs <= modeparms.BASE_NR_BEAMLETS:
        maxnrsbs = modeparms.BASE_NR_BEAMLETS
    elif totnrsbs <= modeparms.BASE_NR_BEAMLETS * 2:
        maxnrsbs = modeparms.BASE_NR_BEAMLETS * 2
    else:
        maxnrsbs = modeparms.BASE_NR_BEAMLETS * 4
    missing_nr_sbs = maxnrsbs - totnrsbs
    if  missing_nr_sbs > 0:
        nrsbs = missing_nr_sbs
        sblo = sbhi + 1
        sbhi = sblo + nrsbs - 1
        freqlo = modeparms.sb2freq(sblo, nz)
        freqhi = modeparms.sb2freq(sbhi, nz)
        obsfileinfo['frequencies'] = numpy.append(obsfileinfo['frequencies'],
                                                  numpy.linspace(freqlo, freqhi, nrsbs))
        totnrsbs += nrsbs
    BSTdirls = os.listdir(BSTfilefolder)
    BSTfiles = [ f for f in BSTdirls if f.endswith('.dat')]

    # Now read the BST pol data
    BST_dtype = numpy.dtype(('f8', (maxnrsbs,)))
    BSTdata = {}
    for BSTpolfile in BSTfiles:
        pol = BSTpolfile[-5]
        with open(os.path.join(BSTfilefolder, BSTpolfile), "rb") as fin:
            BSTdata[pol] = numpy.fromfile(fin, dtype=BST_dtype)
    return BSTdata, obsfileinfo
# END BST related code

# BEGIN SST related code
def parse_sstfolder(SSTfolderpath):
    SSTfoldername = os.path.basename(os.path.normpath(SSTfolderpath))
    obsfolderinfo = {}
    try:
        (Ymd, HMS, rcustr, intstr, durstr, sststr) = SSTfoldername.split('_')
        obsfolderinfo['datetime'] = datetime.datetime.strptime(Ymd+'T'+HMS,'%Y%m%dT%H%M%S')
        obsfolderinfo['rcumode'] =     rcustr[3:]
        obsfolderinfo['integration'] = int(intstr[3:])
        obsfolderinfo['duration'] =    int(durstr[3:])
    except:
        raise ValueError("Folder name not in sst_ext format.")
    return obsfolderinfo

def parse_sstfilename(SSTfilepath):
    SSTfilename = os.path.basename(SSTfilepath)
    obsfileinfo = {}
    try:
        (Ymd, HMS, sststr, rcudatstr) = SSTfilename.split('_')
        obsfileinfo['datetime'] = datetime.datetime.strptime(Ymd+'T'+HMS,'%Y%m%dT%H%M%S')
        #scanrecparms['sst'] = sststr
        (rcu, datext) = rcudatstr[3:].split('.')
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
    obsfolderinfo = parse_sstfolder(SSTfolder)
    files = os.listdir(SSTfolder)
    SSTfiles = [f for f in files if f.endswith('.dat')]
    SSTfiles.sort()
    SSTdatarcu = [None]*len(SSTfiles)
    for sstfile in SSTfiles:
        SSTrcudata, obsfileinfo = readsst(os.path.join(SSTfolder,sstfile))
        SSTdatarcu[obsfileinfo['rcu']] = SSTrcudata
    return SSTdatarcu, obsfolderinfo
# END SST related code

class CVCfiles(object):
    """Provides functionality for covariance cube (CVC) files. (CVC is essentially
    visiblity cubes.) CVC files from a LOFAR station includes ACC and XST files.

    Attributes
    ----------
    dataset: list of array_like
        Each item in list corresponds to one CVC file and is the actual covariance matrix
        cube with shape cvcdim0 x cvcdim1 x cvcdim2.
    fileobstimes: list of str
        A list of datetimes as specified by the CVC filename.
    samptimes: list of datetimes
        The datetime of the visibility matrix sample.

    """
    def __init__(self, datapath):
        self.fileobstimes = []
        self.dataset = []
        self.samptimeset = []
        self.freqset = []
        self.scanrecinfo = ScanRecInfo()
        datapath = os.path.abspath(datapath)
        if os.path.isdir(datapath):
            self.filefolder = datapath
            self._readcvcfolder()
        elif os.path.isfile(datapath):
            self._readcvcfile(datapath)
        else:
            raise ValueError('Path does not exist')

    def get_cvc_dtype(self):
        cvc_dtype = numpy.dtype(('c16', (self.cvcdim1, self.cvcdim2)))
        return cvc_dtype

    def _parse_cvcfile(self, cvcfilepath):
        """Parse the cvc file name. Sets fileobstimes and sets the dimensions of the
        visibility cube, which most generally is: dimtimes * dimrcu0 * dimrcu1.

        The file name should have the format:
            Ymd_HMS_xst.dat for xst file
            Ymd_HMS_acc_nrsbsxnrrcus0xnrrcus1.dat for acc file

        :param cvcfilepath: str
        :return: datatype, filebegindatetime, cvcdim1, cvcdim2
        """
        cvcfilename = os.path.basename(cvcfilepath)
        (Ymd, HMS, cvcextrest) = cvcfilename.split('_',2)
        datatype, restdat = cvcextrest[0:3], cvcextrest[3:]
        (rest, datstr) = restdat.split('.')
        if datatype == 'acc':
            rest = rest.lstrip('_')
            (_nr512, nrrcus0, nrrcus1) = map(int, rest.split('x'))
        else:
            nrrcus0 = 192
            nrrcus1 = 192
        self.cvcdim1 = nrrcus0
        self.cvcdim2 = nrrcus1
        filenamedatetime = datetime.datetime.strptime(Ymd + 'T' + HMS, '%Y%m%dT%H%M%S')
        # NOTE: For ACC, filename is last obstime, while for XST, it is first.
        if datatype == 'acc':
            filebegindatetime = filenamedatetime - datetime.timedelta(seconds=_nr512)
        else:
            filebegindatetime = filenamedatetime
        # Save the observation datetime of file and the cvc dims.
        self.fileobstimes.append(filebegindatetime)
        return datatype, filebegindatetime, self.cvcdim1, self.cvcdim2

    def _parse_cvcfolder(self, cvcfolderpath):
        """Parse the cvc filefolder.

        The filefolder should have the format:
            stnid_Ymd_HMS_rcumode_subband_integration_duration_pointing_cvc

        :param cvcfilepath: str
        :return: scanrecparms
        """
        cvcfoldername = os.path.basename(os.path.abspath(cvcfolderpath))
        obsfolderinfo = {}
        cvcextstr = cvcfoldername.split('_')[-1]
        if cvcextstr == 'xst' or cvcextstr == 'xst-SEPTON':
            try:
                (Ymd, HMS, rcustr, sbstr, intstr, durstr, dirstr, cvcextstr
                ) = cvcfoldername.split('_')
                obsfolderinfo['datetime'] = datetime.datetime.strptime(
                                                Ymd+'T'+HMS, '%Y%m%dT%H%M%S')
                obsfolderinfo['rcumode'] =     rcustr[3:]
                obsfolderinfo['subband'] =     sbstr[2:]
                obsfolderinfo['integration'] = float(intstr[3:])
                obsfolderinfo['duration'] =    float(durstr[3:])
                obsfolderinfo['pointing'] =    dirstr[3:].split(',')
            except:
                raise ValueError("Foldername not in xst_ext format.")
        elif cvcextstr == 'acc':
            dirpat = re.compile(regex_ACCfolder)
            obsdirinfo_m = dirpat.match(cvcfoldername)
            if obsdirinfo_m is None:
                print("Cal error")
                raise ValueError("Calibration directory does not have correct syntax.")
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
            obsfolderinfo['pointing'] = ilisa.observations.directions.std_pointings(
                obsfolderinfo['source'])
        else:
            raise(ValueError, "Folder not expected xst or acc format.")
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
        try:
            self.scanrecinfo.read_scanrec(self.filefolder)
        except Exception as e:
            print(e.message)
            print(e.__doc__)
            print("Warning: Could not read session header. Will try filefolder name...")
            try:
                self.scanrecinfo.scanrecparms = self._parse_cvcfolder(self.filefolder)
            except ValueError:
                self.scanrecinfo.scanrecparms = None
            else:
                print("Read in filefolder meta.")
        cvcdirls = os.listdir(self.filefolder)
        # Select only data files in folder
        self.filenames = self.scanrecinfo.get_ldat_filenames()
        #self.filenames = [f for f in cvcdirls if f.endswith('.dat')]
        self.filenames.sort()  # This enforces chronological order
        for cvcfile in self.filenames:
            # Try to get obsfile header
            try:
                (bfilename, _dat) = cvcfile.split('.')
                hfilename = bfilename + '.h'
                hfilepath = os.path.join(self.filefolder, hfilename)
                obsinfo = ObsInfo.read_ldat_header(hfilepath)
                self.scanrecinfo.add_obs(obsinfo)
            except:
                print("Warning: Couldn't find a header file for {}".format(cvcfile))
            print("Reading cvcfile: {}".format(cvcfile))
            datafromfile, t_begin = self._readcvcfile(
                os.path.join(self.filefolder,cvcfile))
            cvcdim_t, cvcdim_rcu1, cvcdim_rcu2 = datafromfile.shape
            self.dataset.append(datafromfile)

            # Compute time of each autocovariance matrix sample per subband
            integration = self.scanrecinfo.get_integration() #
            obscvm_datetimes = [None] * cvcdim_t
            for t_idx in range(cvcdim_t):
                t_delta = datetime.timedelta(
                    seconds= t_idx * integration
                )
                obscvm_datetimes[t_idx] = t_begin + t_delta
            self.samptimeset.append(obscvm_datetimes)

            # Compute frequency of corresponding time sample
            rcumode = self.scanrecinfo.get_rcumode()
            nz = modeparms.rcumode2nyquistzone(rcumode)
            if self.scanrecinfo.get_datatype() == 'acc':
                freqs = modeparms.rcumode2sbfreqs(rcumode)
            else:
                sb = obsinfo.sb
                freq = modeparms.sb2freq(sb, nz)
                freqs = [freq]*cvcdim_t
            self.freqset.append(freqs)

    def _readcvcfile(self, cvcfilepath):
        """Reads in a single acc or xst data file by filepath and creates
        corresponding sample times.

        The contents of the data file is appended to the object attribute list
        `data`.

        Parameters
        ----------
        cvcfilepath : str
        """
        datatype, filenamedatetime, cvcdim_rcu1, cvcdim_rcu2 =\
            self._parse_cvcfile(cvcfilepath)
        t_begin = filenamedatetime
        # Get cvc data from file.
        cvc_dtype = self.get_cvc_dtype()
        with open(cvcfilepath, 'rb') as fin:
            datafromfile = numpy.fromfile(fin, dtype=cvc_dtype)
        return datafromfile, t_begin

    def save_ldat(self):
        """Save CVC dataset as ldat file."""
        for filenr, cvcfile in enumerate(self.filenames):
            cvcpath = os.path.join(self.filefolder, cvcfile)
            self.dataset[filenr].tofile(cvcpath)

    def getnrfiles(self):
        """Return number of data files in this filefolder."""
        return len(self.filenames)

    def getdata(self, filenr=None):
        """Return the data payload of the filefolder. For ACC each file is a sweep through
        512 frequency. For XST they represent another observation."""
        if filenr is None:
            return self.dataset
        else:
            return self.dataset[filenr]


def cvc2cvpol(cvc):
    """Convert a covariance cube into an array indexed by polarization channels.

    Parameters
    ----------
    cvc: (M,N,N) array (usually M=512 & N=196)
        The Covariance Cube array produced by an International LOFAR
        station when it is in calibration mode. It is the covariance matrices
        of the 196 rcus (98 X-polarized & 98 Y-polarized interleaved) over 512
        subbands.

    Returns
    -------
    cvpol: (2,2,M,N/2,N/2) array
        The same data but indexed into X & Y polarizations. X,Y is index 0,1 resp.
    """
    XX=cvc[:, ::2, ::2]
    YY=cvc[:,1::2,1::2]
    XY=cvc[:, ::2,1::2]
    YX=cvc[:,1::2, ::2]
    cvpol=numpy.array([[XX, XY],[YX,YY]])
    return cvpol


def readacc2bst(anacc2bstfilepath, datformat = 'hdf'):
    """Read an acc2bst file. The fileformat can be either hdf or numpy."""
    anacc2bstfilepath = os.path.abspath(anacc2bstfilepath)
    acc2bstfiledir = os.path.dirname(anacc2bstfilepath)
    anacc2bstfilename = os.path.basename(anacc2bstfilepath)
    (stnid, beginUTCstr, rcuarg, calsrc, durarg, caltabdate, acc2bst, version
     ) = anacc2bstfilename.split('_')
    beginUTC = datetime.datetime.strptime(beginUTCstr, "%Y%m%dT%H%M%S")
    rcumode = rcuarg[3]
    dur = durarg[3:]
    acc2bstvarstr = ['XX', 'YY', 'XY', 'times']
    acc2bstvars = {}
    if datformat == 'hdf':
        hf = h5py.File(anacc2bstfilepath,'r')
        acc2bstvars['XX'] = hf['XX']
        acc2bstvars['XY'] = hf['XY']
        acc2bstvars['YY'] = hf['YY']
        acc2bstvars['times'] = hf['timeaccstart']
        #freqs = hf['frequency']
    else:
        for varstr in acc2bstvarstr:
            acc2bstfilename = '_'.join((beginUTCstr, acc2bst, rcuarg, calsrc, durarg, caltabdate))
            acc2bstfilename += '_'+varstr
            acc2bstfilename += '.npy'
            acc2bstvars[varstr] = numpy.load(acc2bstfiledir+'/'+acc2bstfilename)
    return acc2bstvars, beginUTC, rcumode, calsrc, dur, caltabdate, stnid


def saveacc2bst(bst_pols, filestarttimes, calrunstarttime,
                calrunduration, rcumode, calsrc, calibmeta, stnid, used_autocorr,
                saveformat = "hdf5"):
    """Save acc2bst data to file. Dataformat can be hdf or numpy."""
    (bstXX, bstXY, bstYY) = bst_pols
    version = '5'  # Version of this dataformat
    calrundurationstr = str(int(calrunduration.total_seconds()))
    caltabID = calibmeta['Date']
    # Calculate start of ACC run.
    # Form self describing filename.
    dtlabel = 'acc2bst'
    acc2bstbase = "{}_{}_rcu{}_{}_dur{}_ct{}_v{}_{}".format(
        stnid, calrunstarttime.strftime("%Y%m%dT%H%M%S"), rcumode, calsrc,
        calrundurationstr, caltabID, version, dtlabel)
    #acc2bstsuffix = '.dat'
    pntstr = ilisa.observations.directions.std_pointings(calsrc)
    # Write out the data.
    if saveformat == 'hdf5':
        hf = h5py.File(acc2bstbase+".hdf5", "w")
        freqs = modeparms.rcumode2sbfreqs(rcumode)
        hf.attrs['DataDescription'] = 'LOFAR acc2bst data'
        hf.attrs['StationID'] = stnid
        hf.attrs['calibrationSource'] = calsrc
        hf.attrs['pointing'] = pntstr
        hf.attrs['ObservationStart'] = calrunstarttime.isoformat()
        hf.attrs['ObservationDuration'] = calrundurationstr
        hf.attrs['calibrationTableDate'] = caltabID
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
        numpy.save(acc2bstbase+'_times', filestarttimes)
        numpy.save(acc2bstbase+'_XX', bstXX)
        numpy.save(acc2bstbase+'_XY', bstXY)
        numpy.save(acc2bstbase+'_YY', bstYY)
    return acc2bstbase + "." + saveformat
