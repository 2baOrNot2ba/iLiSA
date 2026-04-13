#!/usr/bin/python3
"""A python module to process beamformed LOFAR data."""
# TobiaC 2022-12-04 (2015-04-14)
import os
import sys
import struct
import math
from multiprocessing import Pool
from itertools import repeat
from datetime import datetime, timedelta, timezone
import argparse
import warnings

import numpy as np
import matplotlib.pyplot as plt

import ilisa.operations
from ilisa.operations.filefolder import obsinfo2filefolder, filefolder2obsinfo
from ilisa.operations.modeparms import (MAX_NRLANES, NRBEAMLETSBYBITS,
                                        nrbits2bitmode)
from ilisa.pipelines.bfbackend import rawfilesinfolder

# BF data/header format constants:
NRTIMS_PACKET = 16  # Number of sample times in packet
FFTSIZE = 1024
NRPOLS = 2  # Number of polarization channels
NRCMPLX = 2  # Every data sample is complex i.e. real & imag so size is doubled
PACKETHDRSZ = 16  # Size of packet header in bytes
BytesBFPacket = NRCMPLX * NRPOLS * 2 * NRTIMS_PACKET * 61 + PACKETHDRSZ  # =7824
PacketH_struct_fmt = 'BBHHBBII'
PacketH_struct_len = struct.calcsize(PacketH_struct_fmt)
cint16 = np.dtype([('re', '<i2'), ('im', '<i2')])
cint8 = np.dtype([('re', '<i1'), ('im', '<i1')])

# PCAP file format data:
BytesPcapHeader = 42  # 42
BytesPcapFooter = 16
BytesPcapFileHeader = 40  # 40
BytesPerPcapPacket = BytesPcapHeader + BytesBFPacket + BytesPcapFooter


def get_samprate(is200mhz):
    """Get samprate"""
    if is200mhz == 1:
        samprate = 200.0e6  # sample rate per second
    else:
        samprate = 160.0e6  # sample rate per second
    return samprate


def bffmtparams(header, checkdrbits=False):
    """\
    Beamform stream format parameters

    Parameters
    ----------
    header : dict
        Packet header as dict.
    checkdrbits :  bool
        Whether to check dynamic range bits or not. Default: False means do not
        check what the bit depth is.

    Returns
    -------
    packetData_dtype :
        dtype of packet data payload.
    cint_smp :
        Complex integer data sample datatype.
    """
    if not hasattr(bffmtparams, 'drbits') or checkdrbits:
        # Get size, in bits, of data samples
        # Seems like bitmode bm sometimes is not set,
        # so use nrbeamlets instead
        #if header['bm'] == 0:
        if header['nrbeamlets'] == 1*61:
            bffmtparams.drbits = 16
        #elif header['bm'] == 1:
        elif header['nrbeamlets'] == 2*61:
            bffmtparams.drbits = 8
        else:
            raise RuntimeError("Only 8 and 16 bit mode supported")
        xrxiyryi16_dtype = np.dtype([('xr', '<i2'), ('xi', '<i2'),
                                        ('yr', '<i2'), ('yi', '<i2')])
        xrxiyryi8_dtype = np.dtype([('xr', '<i1'), ('xi', '<i1'),
                                       ('yr', '<i1'), ('yi', '<i1')])

        if bffmtparams.drbits == 16:
            nrbeamletsperlane = 61  # Number of beamlets in packet
            nrdatasampsBFpacket = NRTIMS_PACKET * header['nrbeamlets']
            # PacketD_struct_fmt = str(2 * 4 * nrdatasampsBFpacket) + 's'
            bffmtparams.cint_smp = cint16
            xrxiyryi_dtype = xrxiyryi16_dtype
        elif bffmtparams.drbits == 8:
            nrbeamletsperlane = 2 * 61  # Number of beamlets in packet
            nrdatasampsBFpacket = NRTIMS_PACKET * header['nrbeamlets']
            # PacketD_struct_fmt = str(2 * 4 * nrdatasampsBFpacket) + 's'
            bffmtparams.cint_smp = cint8
            xrxiyryi_dtype = xrxiyryi8_dtype
        bffmtparams.packetData_dtype = np.dtype((xrxiyryi_dtype,
                                                    (nrbeamletsperlane,
                                                     NRTIMS_PACKET)))
    return bffmtparams.packetData_dtype, bffmtparams.cint_smp


def nrpackets_in_file(bfs_filepath):
    nrpackets = int(np.floor(os.path.getsize(bfs_filepath) / BytesBFPacket))
    return nrpackets


def read_bf_packet(filepointer, keepstruct=False, pcapfile=False,
                   checkdrbits=False):
    """Read a Beamlet packet"""
    pcapheader = None
    if pcapfile:
        pcapheader = filepointer.read(BytesPcapHeader)
    # print(hex(ord(filepointer.read(1))))
    buf = filepointer.read(PacketH_struct_len)
    if buf == b'':
        return '', '', ''

    # Unpack BL header
    header = {}
    try:
        (version, sourceinfo, configuration, station, nrbeamlets, nrblocks,
         timestamps, blocksequencenumber
         ) = struct.unpack(PacketH_struct_fmt, buf)
    except struct.error:
        raise RuntimeError("Could not unpack header: ", buf)
    header['version'] = version
    header['_sourceinfo'] = sourceinfo
    # Parse sourceinfo for addition info
    lanenr = sourceinfo & int('11111', 2)
    error = sourceinfo >> 6 & 1
    is200mhz = sourceinfo >> 7 & 1
    bm = sourceinfo >> 8 & 3
    header['lanenr'] = lanenr
    header['error'] = error
    header['is200mhz'] = is200mhz
    header['bm'] = bm
    header['configuration'] = configuration
    header['station'] = station
    header['nrbeamlets'] = nrbeamlets
    header['nrblocks'] = nrblocks
    header['_timestamps'] = timestamps
    header['datetime'] = datetime.fromtimestamp(timestamps, timezone.utc)
    header['_blocksequencenumber'] = blocksequencenumber
    _microsecflt = (blocksequencenumber * 0.005 * FFTSIZE if is200mhz
                    else blocksequencenumber * 0.00625 * FFTSIZE)
    _microsecfrac, microsecs = math.modf(_microsecflt)
    _nanoseconds = _microsecfrac * 1000
    header['datetime'] = (header['datetime']
                          + timedelta(microseconds=int(microsecs)))
    # Since 0.2 samp/ns doesn't go out evenly for packets with 16*FFTSIZE
    # samples on full second timestamps, I add half an FFTSIZE sample time
    # to the nanosecs calculate for odd timestamps: (Not sure about this)
    header['nanosecs'] = (_nanoseconds
                          + (timestamps % 2)*0.5*FFTSIZE/(0.16+is200mhz*0.04))
    header['datetime64'] = (np.datetime64(header['datetime'])
                            + np.timedelta64(round(header['nanosecs']), 'ns'))
    #header['packno'] =
    # ((timestamps*1000000*(160+is200mhz*40)+oe*512)/1024+blocksequencenumber)/16
    # Get packet format parameters:
    packetData_dtype, cint_smp = bffmtparams(header, checkdrbits)

    # Unpack BL data
    alt = 2
    if alt == 1:
        LofarElemSamps = np.fromfile(filepointer, dtype=packetData_dtype,
                                     count=1)
    else:
        buf = filepointer.read(packetData_dtype.itemsize)
        LofarElemSamps = np.frombuffer(buf, dtype=packetData_dtype)
    xr = (LofarElemSamps[0, :, ]['xr']).squeeze()
    xi = (LofarElemSamps[0, :, ]['xi']).squeeze()
    yr = (LofarElemSamps[0, :, ]['yr']).squeeze()
    yi = (LofarElemSamps[0, :, ]['yi']).squeeze()
    if keepstruct:
        # Keep packet as integer values
        x = np.zeros((nrbeamlets, NRTIMS_PACKET), dtype=cint_smp)
        y = np.zeros((nrbeamlets, NRTIMS_PACKET), dtype=cint_smp)
        x[:, :]['re'] = xr
        x[:, :]['im'] = xi
        y[:, :]['re'] = yr
        y[:, :]['im'] = yi
    else:
        x = xr + 1j * xi
        y = yr + 1j * yi

    pcapfooter = None
    if pcapfile:
        pcapfooter = filepointer.read(BytesPcapFooter)

    return header, x, y


def integrate_int_samples(integration_req, samprate):
    """\
    Return number of samples matching desired integration at given samplerate

    Parameters
    ----------
    integration_req: float
        Requested integration time in seconds.
    samprate: float
        Sample rate in samples per second or Hz.

    Returns
    -------
    integration: float
        Integration time corresponding to whole number of samples integrated
        over.
    integrate_samps: int
        Number of samples intgrated over.
    """
    integrate_samps = (math.floor(integration_req * samprate
                                  / (FFTSIZE * NRTIMS_PACKET))
                       * NRTIMS_PACKET)
    integration = integrate_samps / samprate * FFTSIZE
    return integration, integrate_samps


def _missing_in_sequence(seqdif):
    """\
    Compute number of missing packets from sequence increment

    Parameters
    ----------
    seqdif: str
        Difference in block sequence number since last one.

    Returns
    -------
    missedpackets: int
        Number of missed packets.
    """
    missedpackets = 0
    if seqdif != 16 and seqdif != -195297 and seqdif != -195296:
        if seqdif < 0:
            if (195312 + seqdif) % 16 == 0:
                seqdif += 195312
            elif (195313 + seqdif) % 16 == 0:
                seqdif += 195313
            else:
                raise RuntimeError('Anomalous seqdif: ' + str(seqdif))
        missedpackets = seqdif // 16 - 1
    return missedpackets


def next_bfpacket(bfs_filename, keepstruct=True, padmissing=True,
                  firstpacket=0, checkdrbits=False):
    if bfs_filename.endswith('.pcap'):
        pcapfile = True
        bytesperpacket = BytesPerPcapPacket
        startskip = BytesPcapFileHeader
        endpad = -BytesPcapFooter
    else:
        pcapfile = False
        bytesperpacket = BytesBFPacket
        startskip = 0
        endpad = 0

    # nrPackets =
    # float(os.stat(bfs_filename).st_size - startskip + endpad) / bytesperpacket
    EOF = False
    packetnr = firstpacket
    missed_pkts_tot = 0
    _fil_open = open(bfs_filename, "rb")
    import mmap
    fin = mmap.mmap(_fil_open.fileno(), 0, access=mmap.ACCESS_READ)
    fin = _fil_open
    fin.seek(startskip + packetnr*BytesBFPacket)
    while not EOF:
        header, x, y = read_bf_packet(fin, keepstruct, pcapfile, checkdrbits)
        if header == '':
            EOF = True
            break
        if padmissing:
            if packetnr == firstpacket:
                seqprev = header['_blocksequencenumber']
                seqdif = 16
            else:
                seqdif = header['_blocksequencenumber'] - seqprev
                seqprev = header['_blocksequencenumber']
            missingpackets = _missing_in_sequence(seqdif)
            if missingpackets:
                missed_pkts_tot += missingpackets
                for blkpktidx in range(missingpackets):
                    _header_ = dict(header)  # Todo: update header
                    _x = np.zeros_like(x)
                    _y = np.zeros_like(y)
                    yield _header_, _x, _y
        packetnr += 1
        yield header, x, y
    if padmissing:
        print('Missed packets total:', missed_pkts_tot, packetnr)
    _fil_open.close()


class BFSmeta:
    lcufftlen = 1024
    lcusampfreq = 200e6
    sbsampprd = lcufftlen / lcusampfreq

    def __init__(self, bfs_filename, skip=0):
        self.bfs_filename =  bfs_filename
        # Lookup first and last packet in file to set things up:
        with open(self.bfs_filename, "rb") as fin:
            fin.seek(0, os.SEEK_SET)
            header_first, _x, _y = read_bf_packet(fin)
            fin.seek(0, os.SEEK_END)
            self.nrpkts_actual = fin.tell()//BytesBFPacket
            fin.seek(-BytesBFPacket, os.SEEK_END)
            header_last, _x, _y = read_bf_packet(fin, keepstruct=True)
        self.nrbeamlets = header_first['nrbeamlets']
        self.samprate = get_samprate(header_first['is200mhz'])
        _pkt_dur = np.timedelta64(round(
            NRTIMS_PACKET * FFTSIZE / self.samprate * 1e9), 'ns')
        self.start_time = header_first['datetime64']
        self.stop_time = header_last['datetime64'] + _pkt_dur
        self.dur = (self.stop_time - self.start_time).item() / 1e9
        self.nrpkts_nominal = round(self.dur * self.samprate
                                    / (NRTIMS_PACKET * FFTSIZE))
        self.xy_dtype = _x.dtype
        self.lane = header_first['lanenr']
        obsinfo = filefolder2obsinfo(os.path.dirname(self.bfs_filename))
        beamlet_0 = self.lane*self.nrbeamlets
        beamlet_n = beamlet_0 + self.nrbeamlets
        self.freqs = obsinfo['frequencies'][beamlet_0:beamlet_n]
        self.stnid = obsinfo['station_id']
        bffmtparams(header_first)
        self.bits = bffmtparams.drbits


    def timeaxis(self):
        sampperiod = np.timedelta64(round(
            NRTIMS_PACKET*FFTSIZE/self.samprate * 1e9), 'ns')
        #return np.arange(self.start_time, self.stop_time, sampperiod)
        return np.arange(np.timedelta64(0,'ns'),
                         np.timedelta64(round(self.dur*1e9), 'ns'), sampperiod)

    def freqaxis(self):
        return self.freqs

    def save(self):
        np.savez(self.bfs_filename+'_meta', times=self.timeaxis(),
                 freqs=self.freqaxis(), start_datetime=self.start_time,
                 stnid=self.stnid, bits=self.bits)

    def __str__(self):
        ret = (f"start: {self.start_time}\n"
               f"stop : {self.stop_time}\n"
               f"dur: {self.dur}\n"
               f"pkts_act: {self.nrpkts_actual}\n"
               f"pkts_nom: {self.nrpkts_nominal}\n"
               f"lane: {self.lane}\n"
               f"bits: {self.bits}")
        return ret


def readbfsfile(bfs_filename, bfsmeta, savenpy=True):
    """Read BFS file"""
    _paylodshp = (bfsmeta.nrbeamlets, NRTIMS_PACKET * (bfsmeta.nrpkts_nominal))
    if savenpy:
        x_strm = np.lib.format.open_memmap(bfs_filename + '_X.npy',
                                           dtype=bfsmeta.xy_dtype,
                                           mode='w+',
                                           shape=_paylodshp)
        y_strm = np.lib.format.open_memmap(bfs_filename + '_Y.npy',
                                           dtype=bfsmeta.xy_dtype,
                                           mode='w+',
                                           shape=_paylodshp)
    else:
        x_strm = np.empty(shape=_paylodshp, dtype=complex)
        y_strm = np.empty_like(x_strm)
    # Now read in all the packets
    pktnr = 0
    _leftpc_prev = 0
    for header, x, y in next_bfpacket(bfs_filename, True, True):
        if header == '':
            break
        pktnr_pe = pktnr*NRTIMS_PACKET + NRTIMS_PACKET
        x_strm[:, NRTIMS_PACKET*pktnr:pktnr_pe] = x
        y_strm[:, NRTIMS_PACKET*pktnr:pktnr_pe] = y
        leftpc = (pktnr+1)/bfsmeta.nrpkts_actual * 100
        if leftpc - _leftpc_prev > 1.0:
            _leftpc_prev = leftpc
            print('\r'+str(round(leftpc))+'% completed', end='')
        pktnr += 1
    print()
    return x_strm, y_strm


def parse_bfs_filename(bfs_filepath):
    """\
    Parse the name of a BFS file

    Parameters
    ----------
    bfs_filepath: str
        Path to bfs file.

    Returns
    -------
    prefix: str
        String prefixing BFS filename.
    stnid: str
        Station ID of observation.
    port: int
        Port number used to rec UDP stream.
    druname: str
        Name of Data-Recording Unit.
    fstart_dt: datetime
        Start datetime based on file name.
    rec_chunk_nr: int
        Ordinal of recording chunk.
    """
    bfs_filename = os.path.basename(bfs_filepath)
    _dest, druname, fnametime, rec_chunk_nr = bfs_filename.split('.', 3)
    prefix, stnid, port = _dest.rsplit('_', 2)
    if prefix != 'udp':
        warnings.warn('BFS file prefix is "'+prefix+'" rather than nominal "udp"')
    fstart_dt = (datetime.strptime(fnametime, ilisa.operations.DATETIMESTRFMT)
                 .replace(tzinfo=timezone.utc))
    port = int(port)
    return prefix, stnid, port, druname, fstart_dt, rec_chunk_nr


def format_bfs_filename(prefix, stnid, port, druname, fstart_dt, rec_chunk_nr):
    _dest = '_'.join([prefix, stnid, str(port)])
    bfs_filename = '.'.join([_dest, druname, fstart_dt, rec_chunk_nr])
    return bfs_filename


def parse_bfs_ff(bfs_filefolder):
    """\
    Parse BFS filefolder to get its root, name, bfs files and lane port numbers

    Parameters
    ----------
    bfs_filefolder: str
        Path to BFS filefolder

    Returns
    -------
    bfs_root: str
        Root path containing BFS filefolder
    bfs_ff_name: str
        Name of BFS filefolder
    bfs_files: list of str
        Name of BFS files contained in filefolder
    laneports: list
        Port numbers of lanes associated with the BFS files
    """
    bfs_filefolder = bfs_filefolder.rstrip('/')
    if bfs_filefolder.endswith('_bfs'):
        (bfs_root, bfs_ff_name) = os.path.split(bfs_filefolder)
    else:
        raise RuntimeError('Not BFS filefolder: '+bfs_filefolder)
    udp_files = filter(lambda _f: _f.startswith('udp_'),
                       os.listdir(bfs_filefolder))
    bfs_files = filter(lambda _f: not _f.endswith('.zst') and
                       not _f.endswith('.bin') and not _f.endswith('.npy')
                       and not _f.endswith('.npz'),
                       udp_files)
    bfs_files = list(bfs_files)
    laneports = []
    for _p in bfs_files:
        _, _, port, _, _, _ = parse_bfs_filename(_p)
        laneports.append(port)
    return bfs_root, bfs_ff_name, bfs_files, laneports


def _firstwholesecond(bfs_filename, filestart=None):
    """\
    Find first whole second in BFS file

    Parameters
    ----------
    bfs_filename: str
        Name of BFS file
    filestart: datetime
        Datetime of filestart. If None (default) then determine from filename.

    Returns
    -------
    packetnr: int
        Packet number.
    start_dt: datetime
        Datetime of first whole second.
    """
    packetnr = 0
    print('Searching for first whole second in {}...'
          .format(os.path.basename(bfs_filename)))
    if filestart is None:
        _, _, _, _, filestart, _ = parse_bfs_filename(bfs_filename)
    for header, x, y in next_bfpacket(bfs_filename, padmissing=False,
                                      checkdrbits=True):
        # BFS can have left over packets at beginning of file;
        # they have datetime < filestart
        print('pktnr',header['datetime'],  filestart)
        if header['datetime'] >= filestart:
            break
        packetnr += 1
    print('endpkt', packetnr)
    start_dt = header['datetime']
    return packetnr, start_dt


def fix_filestart_bfsff(bfs_file_path):
    """\
    Fix start of BFS file in filefolder

    Remove leftover packets at start of a BFS file.

    Parameters
    ----------
    bfs_file_path: str
        Path to a BFS file.
    """
    nrpackets = nrpackets_in_file(bfs_file_path)
    pre, stnid, port, drunm, fstart_dt, chunknr = parse_bfs_filename(bfs_file_path)
    startpacket, start_dt = _firstwholesecond(bfs_file_path)
    fnametime = start_dt.strftime(ilisa.operations.DATETIMESTRFMT)
    totnrpackets_new = nrpackets - startpacket + 0
    newfp = format_bfs_filename(pre, stnid, port, drunm, fnametime, chunknr)
    newfp = '__'+newfp
    bfs_ff_path = os.path.dirname(bfs_file_path)
    newfp = os.path.join(bfs_ff_path, newfp)
    print("Copying rest of file...")
    with open(bfs_file_path, 'rb') as fp, open(newfp, 'wb') as nfp:
        fp.seek(startpacket * BytesBFPacket)
        #bufsiz = 1024*1024
        for _p in range(totnrpackets_new):
            _buf = fp.read(BytesBFPacket)
            nfp.write(_buf)


def convert2binary(bfs_filepath):
    """\
    Convert BFS file to binary data file

    Output:
    t=0: Xr[b=0],Xi[b=0],Yr[b=0],Yi[b=0],...,Xr[b=60],Xi[b=60],Yr[b=60],Yi[b=60]
    ...
    t=16: Xr[b=0],Xi[b=0],Yr[b=0],Yi[b=0],...,Xr[b=60],Xi[b=60],Yr[b=60],Yi[b=60]
    (packet boundary)
    t=17: Xr[b=0],Xi[b=0],Yr[b=0],Yi[b=0],...,Xr[b=60],Xi[b=60],Yr[b=60],Yi[b=60]
    ...
    """
    bin_file = bfs_filepath + ".bin"
    fout = open(bin_file, "wb")
    xy = None
    for header, x, y in next_bfpacket(bfs_filepath, True):
        if xy is None:
            packetData_dtype, cint_smp = bffmtparams(header)
            xy = np.zeros((2, header['nrbeamlets'], NRTIMS_PACKET),
                             dtype=cint_smp)
        xy[0, :, :] = x
        xy[1, :, :] = y
        (xy.flatten('F')).tofile(fout)
    fout.close()


def correlate_bfs(bfs_filepath_skip, bst_abspath, integration_req=1.0):
    """\
    Correlate BFS file

    Parameters
    ----------
    bfs_filepath_skip: tuple
        A 2-tuple with bfs filepath and number packet number to skip to.
    bst_abspath: str
        Path to output BST filefolder
    integration_req: int
        Integration length in float seconds
    """
    (bfs_filepath, skip) = bfs_filepath_skip
    (bfs_ff, bfs_filenamebase) = os.path.split(os.path.abspath(bfs_filepath))
    print('Converting', bfs_filepath)
    bst_pathbase = os.path.join(bst_abspath, bfs_filenamebase)
    foutXX = open(bst_pathbase + "_bst_XX.dat_", "wb")
    foutYY = open(bst_pathbase + "_bst_YY.dat_", "wb")
    foutXY = open(bst_pathbase + "_bst_XY.dat_", "wb")

    nrpackets = nrpackets_in_file(bfs_filepath) - skip
    totnrtimsamp = 0
    initialize = True
    for header, x, y in next_bfpacket(bfs_filepath, False,
                                      firstpacket=skip):
        if initialize:
            # Set packet format parameters using first packet
            #    Get samprate
            samprate = get_samprate(header['is200mhz'])
            #    Normalize to Volts^2/second/channel
            nrbeamlets = header['nrbeamlets']
            xx = np.zeros((nrbeamlets,), dtype=float)
            yy = np.zeros((nrbeamlets,), dtype=float)
            xy = np.zeros((nrbeamlets,), dtype=complex)
            integration, integrate_samps = integrate_int_samples(
                integration_req, samprate)
            print("Requested {}s integration, but will instead use {}s".format(
                integration_req, integration))
            normfac = 1.0 / integrate_samps * (samprate / FFTSIZE)
            initialize = False
        xx += np.real(np.sum(x * np.conj(x), axis=-1))
        yy += np.real(np.sum(y * np.conj(y), axis=-1))
        xy += np.sum(x * np.conj(y), axis=-1)
        totnrtimsamp += NRTIMS_PACKET
        if (totnrtimsamp % integrate_samps) == 0:
            print("Integrated {}/{} time samples from:".format(
                totnrtimsamp // integrate_samps,
                nrpackets * NRTIMS_PACKET / integrate_samps), bfs_filenamebase)
            # Normalize to covariance
            xx = xx * normfac
            yy = yy * normfac
            xy = xy * normfac
            # Serial format: t=0...16:  X*X^H[b=0], X*X^H[b=1],...X*X^H[b=60]
            (xx.flatten('F')).tofile(foutXX)
            (yy.flatten('F')).tofile(foutYY)
            (xy.flatten('F')).tofile(foutXY)
            xx = np.zeros((nrbeamlets,), dtype=float)
            yy = np.zeros((nrbeamlets,), dtype=float)
            xy = np.zeros((nrbeamlets,), dtype=complex)
    foutXX.close()
    foutYY.close()
    foutXY.close()


def convert2bst(bfs_filefolder, integration_req=1.0):
    """Recreate bst style data based on the streamed beamformed data

    It also adds the novel combination of X*Y data.
    """
    bfs_root, bfs_ff_name, filen_pl, bmlt_pl = laneinfo_bfs_file(bfs_filefolder)
    obsinfo_bfs = filefolder2obsinfo(bfs_ff_name)

    bfs_filepath_skips = []
    bfs_files = list(filter(None, filen_pl))
    for _filen in bfs_files:
        _filep = os.path.join(bfs_root, bfs_ff_name, _filen)
        _, _, _, _, fstart_dt, _ = parse_bfs_filename(_filep)
        packetstart = 0
        bfs_filepath_skips.append((_filep, packetstart))
    # Create BST filefolder
    obsinfo_bst = dict(obsinfo_bfs)
    obsinfo_bst['integration'] = integration_req
    obsinfo_bst['filenametime'] = datetime.strftime(fstart_dt, '%Y%m%d_%H%M%S')
    lanes_nr = len(bfs_files)
    nrblsperfile = sum(bmlt_pl) // lanes_nr
    _bitsbynrbeamlets = {v//MAX_NRLANES: k for k, v in NRBEAMLETSBYBITS.items()}
    bm = nrbits2bitmode(_bitsbynrbeamlets[nrblsperfile])
    bmlt_pl_bin = int(''.join([str(e//nrblsperfile) for e in bmlt_pl]), 2)
    suffix_head = f'{bm}{bmlt_pl_bin:x}'
    obsinfo_bst['ldat_type'] = suffix_head + 'bst'
    bst_ff_name = obsinfo2filefolder(obsinfo_bst)
    bst_abspath = os.path.join(bfs_root, bst_ff_name)
    os.makedirs(bst_abspath, exist_ok=True)

    with Pool(lanes_nr) as p:
        p.starmap(correlate_bfs, zip(bfs_filepath_skips, repeat(bst_abspath),
                                     repeat(integration_req)))
    # Non multiprocess version:
    # list(starmap(correlate_bfs, zip(bfs_filepaths, repeat(bst_abspath), repeat(integration_req))))
    # Transpose and concatenate the <=4 lanes
    ls_bst = os.listdir(bst_abspath)
    bstdat = {}
    for corr in ['XX', 'YY', 'XY']:
        bstdat[corr] = []
        bst_dtype = np.dtype(('f8', (nrblsperfile)))
        if corr == 'XY':
            bst_dtype = np.dtype(('c16', (nrblsperfile)))
        for _f in sorted(ls_bst):
            if _f.endswith(corr+'.dat_'):
                _bstdat = np.fromfile(os.path.join(bst_abspath, _f),
                                     dtype=bst_dtype)
                bstdat[corr].append(_bstdat.T)
    # The nr of int samps may be different in the different lanes,
    # so pad the missing samps at the ends:
    nrintsmps_max = np.max([bstdat['XX'][lane].shape[-1] for lane in range(lanes_nr)])
    for lane in range(lanes_nr):
        for corr in ['XX', 'YY', 'XY']:
            nrintsmps_cur = bstdat[corr][lane].shape[-1]
            nrintsmps_xtr = nrintsmps_max - nrintsmps_cur
            bstdat[corr][lane] = np.pad(bstdat[corr][lane],
                                        [(0, 0), (0, nrintsmps_xtr)], 'constant',
                                        constant_values=0.)
    # Transpose concatenated bst data and save
    for corr in ['XX', 'YY', 'XY']:
        bstdat[corr] = np.concatenate(bstdat[corr])
        bstcorrlbl = corr[0]
        if corr == 'XY':
            bstcorrlbl += 'Y'
        outfile = os.path.join(bst_abspath, obsinfo_bst['filenametime']
                               + '_bst_00' + bstcorrlbl+'.dat')
        bstdat[corr].T.tofile(outfile)
    # Remove the intermediate untransposed and unconcatenated file
    for _f in os.listdir(bst_abspath):
        if _f.endswith('.dat_'):
            os.remove(os.path.join(bst_abspath, _f))
    return bst_abspath


def convert2npy(bfs_filepath):
    """\
    Convert BFS file to np npy file
    """
    bfsmeta = BFSmeta(bfs_filepath)
    readbfsfile(bfs_filepath, bfsmeta, savenpy=True)
    bfsmeta.save()


def get_packet_h_x_y_fromfile(bfs_filename, packetnr=0):
    """\
    Get a packet from a bfs file (normal or via pcap)
    """
    if bfs_filename.endswith('.pcap'):
        pcapfile = True
        pcktlen = BytesBFPacket + BytesPcapHeader + BytesPcapFooter
        startskip = BytesPcapFileHeader
    else:
        pcapfile = False
        pcktlen = BytesBFPacket
        startskip = 0
    f = open(bfs_filename, "rb")
    if packetnr >= 0:
        f.seek(packetnr * (pcktlen) + startskip, os.SEEK_SET)
    else:
        f.seek(packetnr * pcktlen, os.SEEK_END)
    header, x, y = read_bf_packet(f, pcapfile=pcapfile)
    f.close()
    return header, x, y


def print_packet(header, x, y, print_names=True):
    """\
    Print out a packet
    """
    # Print packet header and X,Y data
    # print header:
    if print_names:
        print("version:", header['version'])
        print("sourceinfo:", header['_sourceinfo'])
        print("  lanenr:", header['lanenr'])
        print("  error:", header['error'])
        print("  is200mhz:", header['is200mhz'])
        print("  bm:", header['bm'])
        print("configuration:", header['configuration'])
        print("station:", header['station'])
        print("nrbeamlets:", header['nrbeamlets'])
        print("nrblocks:", header['nrblocks'])
        #print("datetime:", header['datetime'].isoformat())
        #print("nanosecs:", header['nanosecs'])
        print("datetime64:", header['datetime64'])
        # print("blocksequencenumber:", header['_blocksequencenumber'])
    # print X, Y beamlet data
    for beamletnr in range(header['nrbeamlets']):
        if print_names:
            print("beamlet_{}_X:".format(beamletnr), end=" ")
        print(*x[beamletnr,:])
        if print_names:
            print("beamlet_{}_Y:".format(beamletnr), end=" ")
        print(*y[beamletnr,:])


def plot_packet(header, x, y, print_names=True):
    """\
    Plot a packet
    """
    # print X, Y beamlet data
    plt.subplot(2,2,1)
    plt.pcolormesh(np.abs(x))
    plt.ylabel('Beamlet [#]')
    plt.colorbar()
    plt.title('abs(X)')

    plt.subplot(2,2,2)
    plt.pcolormesh(np.angle(x))
    plt.colorbar()
    plt.title('arg(X)')

    plt.subplot(2,2,3)
    plt.pcolormesh(np.abs(y))
    plt.colorbar()
    plt.xlabel('Time bin [#]')
    plt.ylabel('Beamlet [#]')
    plt.title('abs(Y)')

    plt.subplot(2,2,4)
    plt.pcolormesh(np.angle(y))
    plt.colorbar()
    plt.xlabel('Time bin [#]')
    plt.title('arg(Y)')

    plt.suptitle('BFS packet lane {} @ {}'.format(
        header['lanenr'], header['datetime'].isoformat()))
    plt.show()


def laneinfo_bfs_file(bfs_filefolder):
    """\
    Info about BFS files in folder

    Parameters
    ----------
    bfs_filefolder: str
        Name of BFS file

    Returns
    -------
    bfs_root: str
        BFS root directory.
    bfs_ff_name: str
        BFS filefolder name.
    filen_perlane:

    bmlts_perlane
    """
    bfs_root, bfs_ff_name, bfs_files, laneports = parse_bfs_ff(bfs_filefolder)

    bmlts_perlane = [0  for _ln in range(MAX_NRLANES)]
    filen_perlane = ['' for _ln in range(MAX_NRLANES)]
    for _f in bfs_files:
        _fp = os.path.join(bfs_root, bfs_ff_name, _f)
        startpacketnr, start_dt = _firstwholesecond(_fp)
        print(_f, startpacketnr, start_dt, end=' ')
        _, stnid, port, druname, fstart, chunk_nr = parse_bfs_filename(_fp)
        print(stnid, port, druname, fstart, chunk_nr)
        needs_start_fix = (start_dt-fstart).total_seconds()>=1.0

        if needs_start_fix:
            print('Need to fix start of', end=' ')
        else:
            print('Do not need to fix start of', end=' ')
        print(_fp)
        h,x,y = get_packet_h_x_y_fromfile(_fp)
        bmlts_perlane[h['lanenr']] = h['nrbeamlets']
        filen_perlane[h['lanenr']] = _f
    return bfs_root, bfs_ff_name, filen_perlane, bmlts_perlane

def check_packets(bfs_filename):
    """\
    Check the packets of a BFS file
    """
    # Report on if file start is OK
    startpacketnr, start_dt = _firstwholesecond(bfs_filename)
    print('Start packetnr:', startpacketnr, ' @ ', start_dt)
    # Go through all packets from start
    packetnr = startpacketnr
    missedpackets = 0
    for header, x, y in next_bfpacket(bfs_filename, padmissing=False):
        seq = header['_blocksequencenumber']
        if header['error'] == 1:
            print("Error in packet")
        if packetnr != startpacketnr:
            seqdif = seq - seqprev
        else:
            seqdif = 16
        missedpackets += _missing_in_sequence(seqdif)
        seqprev = seq
        packetnr += 1
    print("Missed packets: ", missedpackets, "/", packetnr, " ",
          100 * missedpackets / float(packetnr), '%')


class BFS_dataset:

    def __init__(self, bfs_ff, segdur=None):
        """Load bfs dataset

        Format: x[<chanbins>, <t_segment_maj>, <t_segment_min>]
        """
        bfsfiles = rawfilesinfolder(bfs_ff)
        self.attrs_ln = []
        for bfsfile in bfsfiles:
            if bfsfile:
                self.attrs_ln.append(np.load(bfsfile + '_meta.npz'))
        nrbeamlets = len(self.attrs_ln[0]['freqs'])
        self.xs = []
        self.ys = []
        for bfsfile in bfsfiles:
            self.xs.append(np.lib.format.open_memmap(bfsfile+'_X.npy',
                                                     mode='r',
                                                     shape=(nrbeamlets, )))
            self.ys.append(np.lib.format.open_memmap(bfsfile+'_Y.npy',
                                                     mode='r',
                                                     shape=(nrbeamlets, )))
        self.seglen = 1
        if segdur is not None:
            self.set_seglen(segdur)

    def set_seglen(self, segdur=1):
        # Set up segment
        if type(segdur) == float:
            self.seglen = int(segdur / BFSmeta.sbsampprd)
        elif type(segdur) == int:
            self.seglen = segdur

    def _sel_freq(self, freq_cntr, ts=(0.0, None)):
        freqs = np.array(
            [self.attrs_ln[_lane]['freqs']
             for _lane in range(len(self.attrs_ln))])
        lane_fr, self.fftbin = np.unravel_index(
            np.argmin((freqs - freq_cntr) ** 2), freqs.shape)
        xseg, yseg = self._get_segmented(lane_fr, ts)
        if self.attrs['bits'] == 8:
            xsegbin = xseg[self.fftbin, ...].view(np.int8).astype(np.float32).view(
                np.complex64)
            ysegbin = yseg[self.fftbin, ...].view(np.int8).astype(np.float32).view(
                np.complex64)
        else:
            xsegbin = xseg[self.fftbin, ...].view(np.int16).astype(np.float32).view(
                np.complex64)
            ysegbin = yseg[self.fftbin, ...].view(np.int16).astype(np.float32).view(
                np.complex64)
        return xsegbin, ysegbin

    def _get_segmented(self, lane_fr, ts):
        """Segment time series"""
        x = self.xs[lane_fr]
        y = self.ys[lane_fr]
        ### Time start
        t_smp0, dur = ts
        if not dur:
            dur = x.shape[1] * BFSmeta.sbsampprd
        seg0 = int(t_smp0 / BFSmeta.sbsampprd) // self.seglen
        segN = int((t_smp0 + dur) / BFSmeta.sbsampprd) // self.seglen
        ### Time end
        nrsbs = x.shape[0]
        ofs = x.shape[1] % self.seglen
        endofs = -ofs if ofs > 0 else None
        # Next lines change to F_CONTIGUOUS
        _xseg = x[:, :endofs].reshape((nrsbs, -1, self.seglen))
        _yseg = y[:, :endofs].reshape((nrsbs, -1, self.seglen))
        del x, y
        xseg = _xseg[:, seg0:segN, :].squeeze()
        yseg = _yseg[:, seg0:segN, :].squeeze()
        del _xseg, _yseg
        nrsegs = xseg.shape[1]
        segdur = self.seglen * BFSmeta.lcufftlen / BFSmeta.lcusampfreq
        self.ts = np.arange(nrsegs) * segdur + t_smp0
        toffset = seg0 * segdur
        self.attrs = self.attrs_ln[lane_fr]
        self.fftfreqs = np.fft.fftshift(np.fft.fftfreq(self.seglen,
                                                       d=BFSmeta.sbsampprd))
        start_dt = datetime.fromtimestamp(
            (self.attrs['start_datetime'] - np.datetime64(0, 's')) \
             / np.timedelta64(1, 's'), timezone.utc)
        # tstart = obsinfo['start_datetime'] + timedelta(seconds=toffset)
        return xseg, yseg

    def sel(self, freq=None, ts=None):
        if freq is not None:
            if type(freq) == float:
                return self._sel_freq(freq, ts)
            lanenr = int(freq)
        xseg, yseg = self._get_segmented(lanenr, ts)
        xseg = np.ascontiguousarray(xseg)
        yseg = np.ascontiguousarray(yseg)
        if self.attrs['bits'] == 8:
            xseg = xseg.view(np.int8)
            yseg = yseg.view(np.int8)
        else:
            xseg = xseg.view(np.int16)
            yseg = yseg.view(np.int16)
        xseg = xseg.astype(np.float32).view(np.complex64)
        yseg = yseg.astype(np.float32).view(np.complex64)
        return xseg, yseg


def __read_use_np(bfsfile):
    """Stub testing alternative way to read BFS data files using numpy"""

    def __pktsinfile(bfsfile):
        """Stub testing alternative way to read BFS data files using numpy"""
        header_dtype = ('header', 'B', PacketH_struct_len)
        xrxiyryi8_dtype = np.dtype([('xr', '<i1'), ('xi', '<i1'),
                                    ('yr', '<i1'), ('yi', '<i1')])
        payload8_dtype = ('payload', xrxiyryi8_dtype, (2 * 61, NRTIMS_PACKET))
        lofpacket8_dtype = np.dtype([header_dtype, payload8_dtype])
        data = np.fromfile(bfsfile, dtype=lofpacket8_dtype)
        datseg = data[:-190].reshape((-1, 512))
        for p in datseg:
            yield p['payload']

    with open('xx.bst', 'wb') as fx, open('yy.bst', 'wb') as fy,\
            open('xy.bst', 'wb') as fxy:
        nr = 0
        for p in __pktsinfile(bfsfile):
            _p=np.ascontiguousarray(p).view('int8').astype(np.float32).view(np.complex64)
            _pp = _p * np.conj(_p)
            ppx = np.real(np.mean(_pp[..., 0::2], axis=(0, -1)))
            ppy = np.real(np.mean(_pp[..., 1::2], axis=(0, -1)))
            ppxy = np.mean(_p[..., 0::2]*_p[..., 1::2], axis=(0,-1))
            ppx.tofile(fx)
            ppy.tofile(fy)
            ppxy.tofile(fxy)
            nr += 1
    print(nr)


def main_cli():
    cmdln_prsr = argparse.ArgumentParser()
    subparsers = cmdln_prsr.add_subparsers(help='sub-command help')

    parser_check = subparsers.add_parser('check', help='Check BFS file')
    parser_check.set_defaults(func=check_packets)
    parser_check.add_argument('bfs_filename', help="BFS filename")

    parser_show = subparsers.add_parser('show', help='Show packet from BFS file')
    parser_show.set_defaults(func=print_packet)
    parser_show.add_argument('-p', '--plot', action='store_true',
                             help="Plot packet")
    parser_show.add_argument('bfs_filename', help="BFS filename")
    parser_show.add_argument('packetnr', type=int, default=0,
                             help="Packet number")

    parser_bin = subparsers.add_parser('bin', help='Convert BFS to binary')
    parser_bin.set_defaults(func=convert2binary)
    parser_bin.add_argument('bfs_filename', help="BFS filename")

    parser_bst = subparsers.add_parser('bst', help='Convert to BST files')
    parser_bst.set_defaults(func=convert2bst)
    parser_bst.add_argument('bfs_ff', help="BFS filefolder path")
    parser_bst.add_argument('integration', type=float,
                            help="Integration in float seconds")

    # Use Numpy for parsing
    parser_bin = subparsers.add_parser('num', help='NOT IMPLEMENTED YET (testing)')
    parser_bin.set_defaults(func=__read_use_np)
    parser_bin.add_argument('bfs_filename', help="BFS filename")

    parser_npy = subparsers.add_parser('npy', help='Convert BFS to np npy')
    parser_npy.set_defaults(func=convert2npy)
    parser_npy.add_argument('bfs_filename', help="BFS filename")

    parser_meta = subparsers.add_parser('meta', help='BFS metadata')
    parser_meta.set_defaults(func=BFSmeta)
    parser_meta.add_argument('bfs_filename', help="BFS filename")

    parser_strip = subparsers.add_parser('fix', help='Strip BFS')
    parser_strip.set_defaults(func=fix_filestart_bfsff)
    parser_strip.add_argument('bfs_ff', help="BFS filefolder")

    parser_info = subparsers.add_parser('info', help='Info on BFS filefolder')
    parser_info.set_defaults(func=laneinfo_bfs_file)
    parser_info.add_argument('bfs_ff', help="BFS filefolder")

    args = cmdln_prsr.parse_args()

    try:
        bfs_filename = args.bfs_filename
    except:
        bfs_filename = ''
    if bfs_filename.endswith('zst'):
        print("Please uncompress the BFS file first")
        sys.exit()

    if args.func == print_packet:
        packetnr =  args.packetnr
        while True:
            header, x, y = get_packet_h_x_y_fromfile(args.bfs_filename,
                                                     packetnr)
            if not header:
                print("Passed end of datafile, nothing to show")
                break
            print("packetnr:", packetnr ,'/', nrpackets_in_file(args.bfs_filename))
            if args.plot:
                plot_packet(header, x, y)
            else:
                print_packet(header, x, y)
                try:
                    input("Press return for next packet, ctrl-c to stop ")
                except KeyboardInterrupt:
                    print('\n')
                    break
            packetnr += 1
    elif args.func == check_packets:
        check_packets(args.bfs_filename)
    elif args.func == convert2binary:
        convert2binary(args.bfs_filename)
    elif args.func == convert2bst:
        convert2bst(args.bfs_ff, args.integration)
    elif args.func == __read_use_np:
        __read_use_np(args.bfs_filename)
    elif args.func == convert2npy:
        convert2npy(args.bfs_filename)
    elif args.func == BFSmeta:
        print(BFSmeta(args.bfs_filename))
    elif args.func == fix_filestart_bfsff:
        bfs_root, bfs_ff_name, bfs_files, laneports = parse_bfs_ff(args.bfs_ff)
        bfs_fps = []
        for _fn in bfs_files:
            bfs_fps.append(os.path.join(bfs_root, bfs_ff_name, _fn))
        with Pool() as pool:
            pool.map(fix_filestart_bfsff, bfs_fps)
    elif args.func == laneinfo_bfs_file:
        _, _, filen_perlane, bmlts_perlane = laneinfo_bfs_file(args.bfs_ff)
        for _ln_idx, _ln in enumerate(filen_perlane):
            print('lane', _ln, ':', bmlts_perlane[_ln_idx])


if __name__ == "__main__":
    main_cli()
