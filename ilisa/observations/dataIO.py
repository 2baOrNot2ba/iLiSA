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
import numpy
import re
import datetime
import h5py
import yaml

import ilisa.observations.modeparms
import ilisa.observations.stationinterface as stationcontrol
import ilisa.observations.modeparms as modeparms

regex_ACCfolder=(
"^(?P<stnid>\w{5})_(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})"
"_(?P<hour>\d{2})(?P<minute>\d{2})(?P<second>\d{2})"
"_rcu(?P<rcumode>\d+)_dur(?P<duration>\d+)(_(?P<calsrc>\w+))?_acc$")
regex_ACCfilename=(
"^(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})"
"_(?P<hour>\d{2})(?P<minute>\d{2})(?P<second>\d{2})"
"_acc_(?P<totnrsb>\d+)x(?P<nrrcu0>\d+)x(?P<nrrcu1>\d+)"
".dat$")
regex_xstfilename=(
"^(?P<year>\d{4})(?P<month>\d{2})(?P<day>\d{2})"
"_(?P<hour>\d{2})(?P<minute>\d{2})(?P<second>\d{2})"
"_rcu(?P<rcumode>\d+)_int(?P<integration>\d+)_dur(?P<duration>\d+)"
"_dir(?P<RAint>\d+).(?P<RAdecimal>\d+),(?P<DECint>\d+).(?P<DECdecimal>\d+),(?P<ref>\s+)"
"_xst.dat$")


def write_acc_header(datapath, starttime, rcuband, pointing):
    accfolderhdrfile = "{}_acc.h".format(starttime)
    with open(os.path.join(datapath, accfolderhdrfile), "w") as f:
        f.write("DataType: acc\n")
        f.write("StartTime: {}\n".format(starttime))
        f.write("RcuBand: {}\n".format(rcuband))
        f.write("Integration: 1.0\n")
        f.write("Duration: 512.0\n")
        f.write("RCUs: [192, 192]\n")
        f.write("Pointing: {}\n".format(pointing))


def write_project_header(datapath, stnid, project, observer):
    projhdrfile = "PROJECT_HEADER.yml"
    with open(os.path.join(datapath, projhdrfile), "w") as f:
        f.write("# LOFAR local station project\n")
        f.write("# Created by {} version {}\n".format("iLiSA", ilisa.__version__))
        f.write("Telescope: {}\n".format("LOFAR"))
        f.write("Station: {}\n".format(stnid))
        f.write("Project: {}\n".format(project))
        f.write("Observer: {}\n".format(observer))


class ObsInfo(object):
    """Contains most import technical information of on an observation."""
    def __init__(self):
        pass

    def setobsinfo(self, LOFARdatTYPE, datetime, rcumode, sb,
                         integration, duration, pointing):
        self.LOFARdatTYPE = LOFARdatTYPE
        self.datetime = datetime
        self.rcumode = str(rcumode)
        self.sb = int(sb)
        self.integration = float(integration)
        self.duration = float(duration)
        self.pointing = pointing

    def setobsinfo_fromname(self, obsdatapath):
        foldername = os.path.basename(os.path.abspath(obsdatapath))
        obsfolderinfo = {}
        dataextstr = foldername.split('_')[-1]
        if dataextstr == 'xst':
            try:
                (Ymd, HMS, rcustr, sbstr, intstr, durstr, dirstr, cvcextstr
                ) = foldername.split('_')
                obsfolderinfo['datetime'] = datetime.datetime.strptime(
                                                Ymd+'T'+HMS, '%Y%m%dT%H%M%S')
                obsfolderinfo['rcumode'] =     rcustr[3:]
                obsfolderinfo['subband'] =     int(sbstr[2:])
                obsfolderinfo['integration'] = float(intstr[3:])
                obsfolderinfo['duration'] =    float(durstr[3:])
                obsfolderinfo['pointing'] =    dirstr[3:].split(',')
                obsfolderinfo['LOFARdatType'] = dataextstr
            except:
                raise ValueError, "Foldername not in correct format."
        elif dataextstr == 'acc':
            dirpat = re.compile(regex_ACCfolder)
            obsdirinfo_m = dirpat.match(foldername)
            if obsdirinfo_m is None:
                raise ValueError, "Calibration directory does not have correct syntax."
            obsdirinfo = obsdirinfo_m.groupdict()
            d0 = datetime.datetime(int(obsdirinfo['year']),
                                   int(obsdirinfo['month']),
                                   int(obsdirinfo['day']),
                                   int(obsdirinfo['hour']),
                                   int(obsdirinfo['minute']),
                                   int(obsdirinfo['second']))
            obsfolderinfo['datetime'] = d0
            obsfolderinfo['rcumode'] = obsdirinfo['rcumode']
            obsfolderinfo['calsrc'] = obsdirinfo['calsrc']
            obsfolderinfo['duration'] = int(obsdirinfo['duration'])
            obsfolderinfo['stnid'] = obsdirinfo['stnid']
        return obsfolderinfo

    def setobsinfo_fromparams(self, lofardatatype, obsdatetime_stamp, beamctl_cmd,
                              rspctl_cmd, caltabinfo="", septonconf=""):
        """Set observation info from parameters"""
        self.LOFARdatTYPE = lofardatatype
        self.datetime = obsdatetime_stamp
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
                #self.rcumode = ''.join(self.rcumode)
                #self.sb = rcusbsep.join(self.sb)
                #self.bl = ','.join(self.bl)
                #self.beamctl_cmd = '\n'.join(self.beamctl_cmd)
            else:
                (antset, rcus, rcumode, beamlets, subbands, anadir, digdir
                ) = modeparms.parse_beamctl_args(beamctl_cmd)
                self.rcumode = [rcumode]
                self.sb = [subbands]
            self.pointing = digdir
        else:
            self.pointing = ""
        self.rspctl_cmd = rspctl_cmd
        rspctl_args = modeparms.parse_rspctl_args(self.rspctl_cmd)
        if self.LOFARdatTYPE != 'bfs':
            self.integration = float(rspctl_args['integration'])
            self.duration = float(rspctl_args['duration'])
        if self.LOFARdatTYPE == 'sst':
            self.sb = ""
        elif self.LOFARdatTYPE.startswith('xst'):
            self.sb = str(rspctl_args['xcsubband'])
        elif self.LOFARdatTYPE == 'bst':
            self.sb = self.sb
        self.caltabinfo = caltabinfo
        self.septonconf = septonconf
        if self.septonconf != "":
            self.rcumode = 5

    def getobsdatapath(self, LOFARdataArchive, folder_name_beamctl_type = True):
        """Create name and destination path for folders (on the DPU) in
        which to save the various LOFAR data products.
        """
        stDataArchive = os.path.join(LOFARdataArchive, self.LOFARdatTYPE)
        stObsEpoch = self.datetime
        st_extName = stObsEpoch
        if folder_name_beamctl_type:
            if self.LOFARdatTYPE == "bst-357":
                st_extName += "_rcu357"
            else:
                if type(self.rcumode) is list:
                    rcumodestr = ''.join([str(rcumode) for rcumode in self.rcumode])
                else:
                    rcumodestr = str(self.rcumode)
                st_extName += "_rcu"+rcumodestr
            if self.sb != [] and self.sb != '':
                st_extName += "_sb"
                st_extName += modeparms.seqlists2slicestr(self.sb)
            if hasattr(self, 'integration'):
                st_extName += "_int"+str(int(self.integration))
            if hasattr(self, 'duration'):
                st_extName += "_dur"+str(int(self.duration))
            if self.LOFARdatTYPE != 'sst':
                if str(self.pointing) != "":
                    st_extName += "_dir"+str(self.pointing)
                else:
                    st_extName += "_dir,,"
        st_extName += "_"+self.LOFARdatTYPE
        datapath = os.path.join(stDataArchive, st_extName)
        return stObsEpoch, datapath

    def parse_bsxST_header(self, headerpath):
        """Parse a bsxST header file. Contains stnid and starttime."""
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
            else:
                contents = yaml.load(hf)
                observer = contents['Observer']
                project = contents['Project']
                datatype = contents['DataType']
                stnid = contents['StationID']
                starttime = contents['StartTime']
                beamctl_line = contents['BeamctlCmds']
                rspctl_lines = contents['RspctlCmds'].split('\n')
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
        # Allocate object attributes
        self.observer = observer
        self.project = project
        self.datatype = datatype
        self.stnid = stnid
        self.starttime = starttime
        self.beamctl_cmd = beamctl_cmd
        self.rspctl_cmd = rspctl_cmd
        self.septonconf = septonconf
        return starttime, stnid, beamctl_cmd

    def isLOFARdatatype(self, obsdatatype):
        """Test if a string 'obsdatatype' is one of iLiSA's recognized LOFAR data types"""
        if (obsdatatype == 'bst' or
            obsdatatype == 'bst-357' or
            obsdatatype == 'sst' or
            obsdatatype == 'xst' or
            obsdatatype == 'xst-SEPTON' or
            obsdatatype == 'bfs'):
            return True
        else:
            return False

    def create_LOFARst_header(self, datapath):
        """Create a header file for LOFAR standalone observation."""
        LOFARstTYPE = self.LOFARdatTYPE
        LOFARstObsEpoch = self.datetime
        if type(self.beamctl_cmd) is not list:
            beamctl_CMD = [self.beamctl_cmd]
        else:
            beamctl_CMD = self.beamctl_cmd
        beamctl_CMD = '\n'.join(beamctl_CMD)
        rspctl_CMD = self.rspctl_cmd
        caltableInfo = self.caltabinfo
        septonconfig = self.septonconf
        def indenttext(txt):
            indentstr = "  "
            return indentstr+txt.replace("\n","\n"+indentstr)
        headerversion = "3"
        if not self.isLOFARdatatype(LOFARstTYPE):
            raise ValueError("Unknown LOFAR statistic type {}.\
                              ".format(LOFARstTYPE))
        LOFARstHeaderFile = LOFARstObsEpoch+"_"+LOFARstTYPE+".h"
        f = open(os.path.join(datapath, LOFARstHeaderFile), "w")
        f.write("# HeaderType: bsxSTdata (YAML)\n")
        f.write("# Header version {}\n".format(headerversion))
        f.write("DataType: {}\n".format(LOFARstTYPE))
        starttime = LOFARstObsEpoch[0:4]+'-'+LOFARstObsEpoch[4:6]+'-'\
                        + LOFARstObsEpoch[6:8]+'T'+LOFARstObsEpoch[9:11]+':'\
                        + LOFARstObsEpoch[11:13]+':'+LOFARstObsEpoch[13:15]
        f.write("StartTime: "+starttime+"\n")
        if septonconfig is not "":
            f.write("SEPTONconfig: {}\n".format(septonconfig))
        f.write("BeamctlCmds: |-\n")
        f.write(indenttext(beamctl_CMD)+"\n")
        # f.write(rspsetup_CMD+"\n")
        # FIX separation of beamctl and rspsetup
        # (Currently rspsetup is in beamctl)
        f.write("RspctlCmds: |-\n")
        f.write(indenttext(rspctl_CMD)+"\n")
        if LOFARstTYPE == 'bst' or LOFARstTYPE == 'bfs':
            f.write("CalTabInfo: |-\n")
            #f.write(indenttext(caltableInfo))
            f.write(str(caltableInfo))
        f.close()


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
        raise ValueError, "Filename not in bst_ext format."
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
        nz = ilisa.observations.modeparms.rcumode2NyquistZone(rcumode)
        freqlo = modeparms.sb2freq(sblo, nz)
        freqhi = modeparms.sb2freq(sbhi, nz)
        obsfileinfo['frequencies'] = numpy.append(obsfileinfo['frequencies'],
                                                  numpy.linspace(freqlo, freqhi, nrsbs))
        totnrsbs += nrsbs

    # When the beamlets allocated is less than the maximum (given by bit depth) the
    # RSPs fill the remaining ones regardless. Hence we have to account for them:
    if totnrsbs > 244:
        maxnrsbs = 2*244
    else:
        maxnrsbs = 244
    missing_sbs = maxnrsbs - totnrsbs
    if  missing_sbs > 0:
        nrsbs = missing_sbs
        sblo = sbhi + 1
        sbhi = sblo + nrsbs
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
        raise ValueError, "Folder name not in sst_ext format."
    return obsfolderinfo

def parse_sstfilename(SSTfilepath):
    SSTfilename = os.path.basename(SSTfilepath)
    obsfileinfo = {}
    try:
        (Ymd, HMS, sststr, rcudatstr) = SSTfilename.split('_')
        obsfileinfo['datetime'] = datetime.datetime.strptime(Ymd+'T'+HMS,'%Y%m%dT%H%M%S')
        #obsfolderinfo['sst'] = sststr
        (rcu, datext) = rcudatstr[3:].split('.')
        obsfileinfo['rcu'] = int(rcu)
    except:
        raise ValueError, "File name not in sst format."
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
    data: list of array_like
        Each item in list corresponds to one CVC file and is the actual covariance matrix
        cube with shape cvcdim0 x cvcdim1 x cvcdim2.
    fileobstimes: list of str
        A list of datetimes as specified by the CVC filename.
    samptimes: list of datetimes
        The datetime of the visibility matrix sample.

    """
    def __init__(self, datapath):
        self.data = []
        self.fileobstimes = []
        self.samptimes = []
        datapath = os.path.abspath(datapath)
        if os.path.isdir(datapath):
            self.filefolder = datapath
            self._readcvcfolder()
        elif os.path.isfile(datapath):
            self._readcvcfile(datapath)
        else:
            raise ValueError('Path does not exist')

    def _parse_cvcfile(self, cvcfilepath):
        """Parse the cvc file name. Sets fileobstimes and sets the dimensions of the
        visibility cube, which most generally is: dimtimes * dimrcu0 * dimrcu1.

        The file name should have the format:
            Ymd_HMS_xst.dat for xst file
            Ymd_HMS_acc_nrsbsxnrrcus0xnrrcus1.dat for acc file

        :param cvcfilepath: str
        :return: filenamedatetime, cvcdim0, cvcdim1, cvcdim2
        """
        cvcfilename = os.path.basename(cvcfilepath)
        (Ymd, HMS, cvcextrest) = cvcfilename.split('_',2)
        cvcext, restdat = cvcextrest[0:3], cvcextrest[3:]
        (rest, datstr) = restdat.split('.')
        if cvcext == 'acc':
            rest = rest.lstrip('_')
            (nrsbs, nrrcus0, nrrcus1) = map(int, rest.split('x'))
            self.cvcdim0=nrsbs
        else:
            self.cvcdim0 = 0
            nrrcus0 = 192
            nrrcus1 = 192
        self.cvcdim1 = nrrcus0
        self.cvcdim2 = nrrcus1
        filenamedatetime = datetime.datetime.strptime(Ymd + 'T' + HMS, '%Y%m%dT%H%M%S')
        # Save the observation datetime of file and the cvc dims.
        self.fileobstimes.append(filenamedatetime)
        return filenamedatetime, self.cvcdim0, self.cvcdim1, self.cvcdim2

    def _parse_cvcfolder(self, cvcfolderpath):
        """Parse the cvc filefolder.

        The filefolder should have the format:
            stnid_Ymd_HMS_rcumode_subband_integration_duration_pointing_cvc

        :param cvcfilepath: str
        :return: obsfolderinfo
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
            d0 = datetime.datetime(int(obsdirinfo['year']),
                                   int(obsdirinfo['month']),
                                   int(obsdirinfo['day']),
                                   int(obsdirinfo['hour']),
                                   int(obsdirinfo['minute']),
                                   int(obsdirinfo['second']))
            obsfolderinfo['datetime'] = d0
            obsfolderinfo['rcumode'] = obsdirinfo['rcumode']
            obsfolderinfo['subband'] = '0:511'
            obsfolderinfo['integration'] = 1.0
            obsfolderinfo['duration'] = int(obsdirinfo['duration'])
            obsfolderinfo['calsrc'] = obsdirinfo['calsrc']
            obsfolderinfo['pointing'] = modeparms.stdPointings(obsfolderinfo['calsrc'])
            obsfolderinfo['stnid'] = obsdirinfo['stnid']
        else:
            raise(ValueError, "Folder not expected xst or acc format.")
        obsfolderinfo['cvc-type'] = cvcextstr
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
        self.obsinfos = []
        try:
            self.obsfolderinfo = self._parse_cvcfolder(self.filefolder)
        except ValueError:
            self.obsfolderinfo = None
        cvcdirls = os.listdir(self.filefolder)
        # Select only data files in folder
        cvcfiles = [f for f in cvcdirls if f.endswith('.dat')]
        cvcfiles.sort()  # This enforces chronological order
        self.filenames = []
        for cvcfile in cvcfiles:
            self.filenames.append(cvcfile)
            print("Reading cvcfile: {}".format(cvcfile))
            self._readcvcfile(os.path.join(self.filefolder,cvcfile))
            if self.obsfolderinfo['cvc-type'] != 'acc':
                try:
                    (d,t, _rest) = cvcfile.split('_', 2)
                    hfilename = '{}_{}_{}.h'.format(d, t, self.obsfolderinfo['cvc-type'])
                    hfilepath = os.path.join(self.filefolder, hfilename)
                    obsinfo = ObsInfo()
                    obsinfo.parse_bsxST_header(hfilepath)
                    self.obsinfos.append(obsinfo)
                except:
                    print("Couldn't find a header file for {}".format(cvcfile))

    def _readcvcfile(self, cvcfilepath):
        """Reads in a single acc or xst data file by filepath and creates
        corresponding sample times.

        The contents of the data file is appended to the object attribute list
        `data`.

        Parameters
        ----------
        cvcfilepath : str
        """
        filenamedatetime, cvcdim0, cvcdim1, cvcdim2 = self._parse_cvcfile(cvcfilepath)
        integration = self.obsfolderinfo['integration']  # FIXME get this from st header
        t_end = filenamedatetime
        # Get cvc data from file.
        cvc_dtype = numpy.dtype(('c16', (192,192)))
        with open(cvcfilepath, "rb") as fin:
            datafromfile = numpy.fromfile(fin, dtype=cvc_dtype)
        self.data.append(datafromfile)

        # Compute time of each autocovariance matrix sample per subband
        obsaccdates = [None] * cvcdim0
        for sb in range(cvcdim0):
            # Previously forgot to add 1 in this formula:
            sbobstimedelta = datetime.timedelta(
                seconds=(sb - cvcdim0 + 1) * integration
            )
            obsaccdates[sb] = t_end + sbobstimedelta
        self.samptimes.append(obsaccdates)

    def getnrfiles(self):
        """Return number of data files in this filefolder."""
        return len(self.filenames)

    def getdata(self, filenr=None):
        """Return the data payload of the filefolder. For ACC each file is a sweep through
        512 frequency. For XST they represent another observation."""
        if filenr is None:
            return self.data
        else:
            return self.data[filenr]

    def getobsfolderinfo(self):
        if self.obsinfos:
            return self.obsinfos
        elif self.obsfolderinfo:
            return self.obsfolderinfo
        else:
            raise RuntimeError("No metadata available")


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


def saveacc2bst((bstXX, bstXY, bstYY), filestarttimes, calrunstarttime,
                calrunduration, rcumode, calsrc, calibmeta, stnid, used_autocorr,
                saveformat = "hdf5"):
    """Save acc2bst data to file. Dataformat can be hdf or numpy."""
    version = '3'  # Version of this dataformat
    calrundurationstr = str(int(calrunduration.total_seconds()))
    caltabID = calibmeta['Date']
    # Calculate start of ACC run.
    # Form self describing filename.
    dtlabel = 'acc2bst'
    acc2bstbase = stnid+'_'+calrunstarttime.strftime("%Y%m%dT%H%M%S")\
                  +'_rcu'+rcumode +'_'+calsrc+'_dur'+calrundurationstr\
                  +'_ct'+caltabID+'_v'+version+'_'+dtlabel
    #acc2bstsuffix = '.dat'
    pntstr = modeparms.stdPointings(calsrc)
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
