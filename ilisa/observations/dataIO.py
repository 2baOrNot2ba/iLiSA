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


import sys
import os
import numpy
import re
import datetime
import h5py
import yaml
import ilisa.observations.stationcontrol as stationcontrol
import ilisa.observations.observing as observing


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


def parse_bsxST_header(headerpath):
    """Parse a bsxST header file. Contains stnid and starttime."""
    # TODO extract beam_ctl commands and CalTable info.
    files = os.listdir(headerpath)
    headerfiles = [f for f in files if f.endswith('.h')]
    headerfile = os.path.join(headerpath,headerfiles.pop())
    stnid = None
    starttime = None
    with open(headerfile,'r') as hf:
        for hline in hf:
            if "Header version" in hline:
                headerversion = hline.split()[-1]
    with open(headerfile, 'r') as hf:
        if headerversion == '1':
            for line in hf:
                if "StationID" in line:
                    label, stnid = line.split('=')
                    stnid = stnid.strip()
                if "StartTime" in line:
                    label, starttime = line.split('=')
                    starttime = starttime.strip()
                if "beamctl" in line:
                    # HACK
                    beamctl_line = line
        else:
            contents = yaml.load(hf)
            stnid = contents['StationID']
            starttime = contents['StartTime']
            beamctl_line = contents['BeamctlCmds']
    multishellcmds = beamctl_line.split('&')
    beamctl_cmd = multishellcmds[0]
    (antennaset, rcus, rcumode, beamlets, subbands, anadir, digdir) \
        = stationcontrol.parse_beamctl_args(beamctl_cmd)
    beamctl = {'antennaset': antennaset,
               'rcus': rcus,
               'rcumode': rcumode,
               'beamlets': beamlets,
               'subbands': subbands,
               'anadir': anadir,
               'digdir': digdir}
    return starttime, stnid, beamctl


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
        obsfileinfo['subbands'] =     sbstr[2:]
        obsfileinfo['integration'] = int(intstr[3:])
        obsfileinfo['duration'] =    int(durstr[3:])
        obsfileinfo['pointing'] =    dirstr[3:].split(',')
        #polextstr
        #(obsfileinfo['pol'], datext) = polextstr[2:].split('.')
    except:
        raise ValueError, "Filename not in bst_ext format."
    return obsfileinfo


def readbstfolder(BSTfilefolder):
    obsfileinfo = parse_bstfolder(BSTfilefolder)
    (sblo,sbhi) = obsfileinfo['subbands'].split(':')
    (sblo, sbhi) = (int(sblo), int(sbhi))
    obsfileinfo['sblo'] = sblo
    obsfileinfo['sbhi'] = sbhi
    nrsbs = sbhi-sblo+1
    BSTdirls = os.listdir(BSTfilefolder)
    BSTfiles = [ f for f in BSTdirls if f.endswith('.dat')]

    # Now read the BST pol data
    BST_dtype = numpy.dtype(('f8', (nrsbs,)))
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

class XSTdata(object):
    """Provides functionality for XST data."""
    def __init__(self, datapath):
        if os.path.isdir(datapath):
            self.readxstfolder(datapath)
        elif os.path.isfile(datapath):
            self.readxstfile(datapath)
        else:
            raise ValueError('Path does not exist')

    def parse_xstfolder(self, XSTfilepath):
        """Parse the xst filefolder.

        The filefolder should have the format:
            Ymd_HMS_rcumode_subband_integration_duration_pointing_xst

        :param XSTfilepath: str
        :return: obsfileinfo
        """
        XSTfilename = os.path.basename(XSTfilepath)
        obsfileinfo = {}
        try:
            (Ymd, HMS, rcustr, sbstr, intstr, durstr, dirstr, xstextstr
            ) = XSTfilename.split('_')
            obsfileinfo['datetime'] = datetime.datetime.strptime(Ymd+'T'+HMS,'%Y%m%dT%H%M%S')
            obsfileinfo['rcumode'] =     rcustr[3:]
            obsfileinfo['subband'] =     int(sbstr[2:])
            obsfileinfo['integration'] = float(intstr[3:])
            obsfileinfo['duration'] =    float(durstr[3:])
            obsfileinfo['pointing'] =    dirstr[3:].split(',')
        except:
            raise ValueError, "Filename not in xst_ext format."
        return obsfileinfo

    def readxstfolder(self, XSTfilefolder):
        """Read in XST data from an filefolder.

        The filefolder name should have the format as specified in the parse_xstfolder() method.
        The contents of the data file is stored in the class attribute:
           XSTdata : (N,192,192)
        where N is the number of time samples.

        Parameters
        ----------
        XSTfilefolder : str
            The name of the XST filefolder.
        """
        self.obsfileinfo = self.parse_xstfolder(XSTfilefolder)
        XSTdirls = os.listdir(XSTfilefolder)
        XSTfiles = [ f for f in XSTdirls if f.endswith('.dat')]
        XSTfile = XSTfiles[0]    # FIX probably should warn if more than 1 xst file
        self.readxstfile(os.path.join(XSTfilefolder,XSTfile))

    def readxstfile(self, XSTfilepath):
        """Reads in a single xst data file by filepath.

        The contents of the data file is stored in the class attribute:
           XSTdata : (N,192,192)
        where N is the number of time samples.

        Parameters
        ----------
        XSTfilepath : str
        """
        XST_dtype = numpy.dtype(('c16', (192,192)))
        with open(XSTfilepath, "rb") as fin:
            self.XSTdata = numpy.fromfile(fin, dtype=XST_dtype)

    def getdata(self):
        return self.XSTdata

    def getobsfileinfo(self):
        return self.obsfileinfo

    def xst2XY(self, xst):
        """Return polarized components of flat XST data.

        Parameters
        ----------
        xst : (N,M) array of complex

        Returns
        -------
        XX, YY, XY, YX: (N/2,M/2) array of complex
        """
        XX = xst[::2, ::2]
        YY = xst[1::2,1::2]
        XY = xst[::2,1::2]
        YX = xst[1::2, ::2]
        return XX, YY, XY, YX


# BEGIN ACC related code
def parse_accfolder(caldumpdir):
    """Parse an ACC calibration folder name for obs parameters."""
    caldumpdir = os.path.basename(os.path.abspath(caldumpdir))
    dirpat = re.compile(regex_ACCfolder)
    obsdirinfo_m = dirpat.match(caldumpdir)
    if obsdirinfo_m is None:
        raise ValueError, "Calibration directory does not have correct syntax."
    obsdirinfo = obsdirinfo_m.groupdict()
    d0 = datetime.datetime(int(obsdirinfo['year']),
                           int(obsdirinfo['month']),
                           int(obsdirinfo['day']),
                           int(obsdirinfo['hour']),
                           int(obsdirinfo['minute']),
                           int(obsdirinfo['second']))
    return d0, obsdirinfo['rcumode'], obsdirinfo['calsrc'], \
           int(obsdirinfo['duration']), obsdirinfo['stnid']


def parse_accfilename(filepath):
    filename = os.path.basename(filepath)
    caldumpdir = os.path.dirname(os.path.abspath(filepath))
    d0, rcumode, calsrc, duration, stnid = parse_accfolder(caldumpdir)
    accextpat = re.compile(regex_ACCfilename)
    obsfileinfo_m = accextpat.match(filename)
    if obsfileinfo_m is None:
        raise ValueError, "Filename not in nominal ACC format."
    obsfileinfo = obsfileinfo_m.groupdict()
    t_end = datetime.datetime(int(obsfileinfo['year']),
                           int(obsfileinfo['month']),
                           int(obsfileinfo['day']),
                           int(obsfileinfo['hour']),
                           int(obsfileinfo['minute']),
                           int(obsfileinfo['second']))
    return t_end, rcumode, calsrc, int(obsfileinfo['totnrsb']), \
           int(obsfileinfo['nrrcu0']), int(obsfileinfo['nrrcu1']), stnid


def readacc(ACCfilepath):
    """Read an ACC data file."""
    integtime = 1 # integration time of an ACC sb is 1second.
    # Parse the ACC filename for datetime of observation.
    t_end, rcumode, calsrc, totnrsb, nrrcu0, nrrcu1, stnid =\
                                                  parse_accfilename(ACCfilepath)
    # Now read the ACC data
    ACC_dtype = numpy.dtype(('c16',
                      (nrrcu0,nrrcu1)))
    fin = open(ACCfilepath, "rb")
    ACCdata = numpy.fromfile(fin, dtype=ACC_dtype, count=totnrsb)
    fin.close()

    # Compute time of each autocovariance matrix sample per subband
    obsaccdates = [None] * totnrsb
    for sb in range(totnrsb):
        # Previously forgot to add 1 in this formula:
        sbobstimedelta = datetime.timedelta(
                           seconds=(sb-totnrsb+1)*integtime
                         )
        obsaccdates[sb] = t_end+sbobstimedelta
    return ACCdata, obsaccdates


def readaccfolder(ACCfolderpath):
    """Read an ACC folder."""
    obsfolderinfo = parse_accfolder(ACCfolderpath)
    ACCfiles = os.listdir(ACCfolderpath)
    ACCfiles.sort()
    nrACCfiles = len(ACCfiles)
    ACCdata = [None]*nrACCfiles
    obsaccdates = [None]*nrACCfiles
    for n, accfile in enumerate(ACCfiles):
        ACCdata[n], obsaccdates[n] = readacc(os.path.join(ACCfolderpath, accfile))
    return ACCdata, obsaccdates, obsfolderinfo


def readacc2bst(anacc2bstfilepath, datformat = 'hdf'):
    """Read an acc2bst file. The fileformat can be either hdf or numpy."""
    anacc2bstfilepath = os.path.abspath(anacc2bstfilepath)
    acc2bstfiledir = os.path.dirname(anacc2bstfilepath)
    anacc2bstfilename = os.path.basename(anacc2bstfilepath)
    (stnid, beginUTCstr, rcuarg, calsrc, durarg, caltabdate, acc2bst
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
                calrunduration, rcumode, calsrc, calibmeta, stnid,
                saveformat = "hdf5"):
    """Save acc2bst data to file. Dataformat can be hdf or numpy."""
    calrundurationstr = str(int(calrunduration.total_seconds()))
    caltabID = calibmeta['Date']
    # Calculate start of ACC run.
    # Form self describing filename.
    dtlabel = 'acc2bst'
    acc2bstbase = stnid+'_'+calrunstarttime.strftime("%Y%m%dT%H%M%S")\
                  +'_rcu'+rcumode +'_'+calsrc+'_dur'+calrundurationstr\
                  +'_ct'+caltabID+'_'+dtlabel
    #acc2bstsuffix = '.dat'
    pntstr = observing.stdPointings(calsrc)
    # Write out the data.
    if saveformat == 'hdf5':
        hf = h5py.File(acc2bstbase+".hdf5", "w")
        freqs = stationcontrol.rcumode2sbfreqs(rcumode)
        hf.attrs['DataDescription'] = 'LOFAR acc2bst data'
        hf.attrs['StationID'] = stnid
        hf.attrs['calibrationSource'] = calsrc
        hf.attrs['pointing'] = pntstr
        hf.attrs['ObservationStart'] = calrunstarttime.isoformat()
        hf.attrs['ObservationDuration'] = calrundurationstr
        hf.attrs['calibrationTableDate'] = caltabID
        hf.attrs['version'] = '2.0'
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
# END ACC related code




if __name__=="__main__":
    ACCdata = readacc(sys.argv[1])
