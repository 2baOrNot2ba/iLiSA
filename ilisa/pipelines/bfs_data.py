#!/usr/bin/python3
"""A python module to process beamformed LOFAR data."""
# TobiaC 2022-12-04 (2015-04-14)
import os
import sys
import struct
import math

import datetime
import argparse
import numpy as np
import matplotlib.pyplot as plt

import ilisa.operations.data_io as dio
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


def bffmtparams(header):
    if not hasattr(bffmtparams, 'drbits'):
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


def read_bf_packet(filepointer, keepstruct=False, pcapfile=False):
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
    header['datetime'] = datetime.datetime.utcfromtimestamp(timestamps)
    header['_blocksequencenumber'] = blocksequencenumber
    _microsecflt = (blocksequencenumber * 0.005 * FFTSIZE if is200mhz
                    else blocksequencenumber * 0.00625 * FFTSIZE)
    _microsecfrac, microsecs = math.modf(_microsecflt)
    _nanoseconds = _microsecfrac * 1000
    header['datetime'] = (header['datetime']
                          + datetime.timedelta(microseconds=int(microsecs)))
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
    packetData_dtype, cint_smp = bffmtparams(header)

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


def next_bfpacket(bfs_filename, keepstruct=True, padmissing=True):
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
    fin = open(bfs_filename, "rb")
    EOF = False
    packetnr = 0
    missed_pkts_tot = 0
    fin.seek(startskip)
    while not EOF:
        # print(str(packetNr)+" / "+str(nrPackets))
        header, x, y = read_bf_packet(fin, keepstruct, pcapfile)
        if header == '':
            EOF = True
            break
        if padmissing:
            if packetnr == 0:
                seqprev = header['_blocksequencenumber']
                seqdif = 16
            else:
                seqdif = header['_blocksequencenumber'] - seqprev
                seqprev = header['_blocksequencenumber']
            if seqdif != 16 and seqdif != -195297 and seqdif != -195296:
                missingpackets = seqdif//16 - 1
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
    fin.close()


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
        obsinfo = dio.filefolder2obsinfo(os.path.dirname(self.bfs_filename))
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


def convert2bst(bfs_filepath, integration_req=1.0):
    """Recreate bst style data based on the streamed beamformed data

    It also adds the novel combination of X*Y data.
    """
    (bfs_ff, bfs_filenamebase) = os.path.split(os.path.abspath(bfs_filepath))
    if bfs_ff.endswith('_bfs'):
        (bfs_path, bfs_ff_name) = os.path.split(bfs_ff)
        obsinfo_bsf = dio.filefolder2obsinfo(bfs_ff_name)
        obsinfo_bst = dict(obsinfo_bsf)
        obsinfo_bst['ldat_type'] = 'bst'
        obsinfo_bst['integration'] = integration_req
        bst_ff_name = dio.obsinfo2filefolder(obsinfo_bst)
        bst_abspath = os.path.join(bfs_path, bst_ff_name)
        print(obsinfo_bst['max_nr_bls']//4)
        os.mkdir(bst_abspath)
    bfs_pathbase = os.path.join(bst_abspath, bfs_filenamebase)
    foutXX = open(bfs_pathbase + "_bst_XX.dat", "wb")
    foutYY = open(bfs_pathbase + "_bst_YY.dat", "wb")
    foutXY = open(bfs_pathbase + "_bst_XY.dat", "wb")

    nrpackets = os.path.getsize(bfs_filepath)/BytesBFPacket
    totnrtimsamp = 0
    initialize = True
    for header, x, y in next_bfpacket(bfs_filepath, False):
        if initialize:
            # Set packet format parameters using first packet
            #    Get samprate
            samprate = get_samprate(header['is200mhz'])
            #    Normalize to Volts^2/second/channel
            nrbeamlets = header['nrbeamlets']
            xx = np.zeros((nrbeamlets,), dtype=float)
            yy = np.zeros((nrbeamlets,), dtype=float)
            xy = np.zeros((nrbeamlets,), dtype=complex)
            integrate_samps = (math.floor(integration_req * samprate
                                          / (FFTSIZE * NRTIMS_PACKET))
                               * NRTIMS_PACKET)
            integration = integrate_samps / samprate * FFTSIZE
            print("Requested {}s integration, but will instead use {}s".format(
                integration_req, integration))
            normfac = 1.0 / integrate_samps * (samprate / FFTSIZE)
            initialize = False
        xx += np.real(np.sum(x * np.conj(x), axis=-1))
        yy += np.real(np.sum(y * np.conj(y), axis=-1))
        xy += np.sum(x * np.conj(y), axis=-1)
        totnrtimsamp += NRTIMS_PACKET
        if (totnrtimsamp % integrate_samps) == 0:
            print("Integrated {}/{} time samples".format(
                totnrtimsamp//integrate_samps,
                nrpackets*NRTIMS_PACKET/integrate_samps))
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
    f.seek(packetnr * (pcktlen) + startskip)
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
        print("datetime:", header['datetime'].isoformat())
        print("nanosecs:", header['nanosecs'])
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


def check_packets(bfs_filename):
    """\
    Check the packets of a BFS file
    """
    packetnr = 0
    missedpackets = 0
    for header, x, y in next_bfpacket(bfs_filename, padmissing=False):
        seq = header['_blocksequencenumber']
        if header['error'] == 1:
            print("Error in packet")
        if packetnr != 0:
            seqdif = seq - seqprev
        else:
            seqdif = 16
        if seqdif != 16 and seqdif != -195297 and seqdif != -195296:
            missedpackets += seqdif//16-1
            if False:
                print(seqdif)
        seqprev = seq
        packetnr += 1
    print("Missed packets: ", missedpackets, "/", packetnr, " ",
          100 * missedpackets / float(packetnr), '%')


class BFS_dataset:

    def __init__(self, bfs_ff):
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

    def segment(self, segdur=None):
        self.seglen = 1
        if segdur:
            if type(segdur) == float:
                self.seglen = int(segdur / BFSmeta.sbsampprd)
            elif type(segdur) == int:
                self.seglen = segdur
        return self

    def sel(self, freq_cntr, ts=(0.0, None)):
        freqs = np.array(
            [self.attrs_ln[_lane]['freqs']
             for _lane in range(len(self.attrs_ln))])
        lane_fr, self.fftbin = np.unravel_index(
            np.argmin((freqs - freq_cntr) ** 2), freqs.shape)
        x = self.xs[lane_fr]
        y = self.ys[lane_fr]
        ### Time start
        t_smp0, dur = ts
        if not dur:
            dur  = x.shape[1] * BFSmeta.sbsampprd
        seg0 = int(t_smp0 / BFSmeta.sbsampprd) // self.seglen
        segN = int((t_smp0 + dur) / BFSmeta.sbsampprd) // self.seglen
        ### Time end
        nrsbs = x.shape[0]
        ofs = x.shape[1] % self.seglen
        endofs = -ofs if ofs > 0 else None
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
        start_dt = datetime.datetime.utcfromtimestamp(
            (self.attrs['start_datetime'] - np.datetime64(0, 's')) / np.timedelta64(1, 's'))
        # tstart = obsinfo['start_datetime']+datetime.timedelta(seconds=toffset)
        if self.attrs['bits'] == 8:
            xsegbin = xseg[self.fftbin, ...].view(np.int8).astype(np.float32).view(
                np.complex64)
            ysegbin = yseg[self.fftbin, ...].view(np.int8).astype(np.float32).view(
                np.complex64)
        else:
            xsegbin = xseg[self.fftbin, ...].view(np.int16).astype(np.float32).view(
                np.complex64)
            ysegbin = xseg[self.fftbin, ...].view(np.int16).astype(np.float32).view(
                np.complex64)
        return xsegbin, ysegbin


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
    parser_bst.add_argument('bfs_filename', help="BFS filename")
    parser_bst.add_argument('integration', type=float,
                            help="Integration in float seconds")

    parser_npy = subparsers.add_parser('npy', help='Convert BFS to np npy')
    parser_npy.set_defaults(func=convert2npy)
    parser_npy.add_argument('bfs_filename', help="BFS filename")

    parser_meta = subparsers.add_parser('meta', help='BFS metadata')
    parser_meta.set_defaults(func=BFSmeta)
    parser_meta.add_argument('bfs_filename', help="BFS filename")

    args = cmdln_prsr.parse_args()

    if args.bfs_filename.endswith('zst'):
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
            print("packetnr:", packetnr)
            if args.plot:
                plot_packet(header, x, y)
            else:
                print_packet(header, x, y)
                input("Press return for next packet, ctrl-c to stop ")
            packetnr += 1
    elif args.func == check_packets:
        check_packets(args.bfs_filename)
    elif args.func == convert2binary:
        convert2binary(args.bfs_filename)
    elif args.func == convert2bst:
        integration = 1.0
        convert2bst(args.bfs_filename, args.integration)
    elif args.func == convert2npy:
        convert2npy(args.bfs_filename)
    elif args.func == BFSmeta:
        print(BFSmeta(args.bfs_filename))


if __name__ == "__main__":
    main_cli()
