#!/usr/bin/python3
"""A python module to process beamformed LOFAR data."""
# TobiaC 2022-12-04 (2015-04-14)
import os
import os.path
import sys
import struct
import math

import matplotlib.pyplot as plt
import numpy
import datetime
import argparse

# BF data/header format constants:
import numpy as np
import ilisa.operations.data_io as dio

NRTIMS_PACKET = 16  # Number of sample times in packet
FFTSIZE = 1024
NRPOLS = 2  # Number of polarization channels
NRCMPLX = 2  # Every data sample is complex i.e. real & imag so size is doubled
PACKETHDRSZ = 16  # Size of packet header in bytes
BytesBFPacket = NRCMPLX * NRPOLS * 2 * NRTIMS_PACKET * 61 + PACKETHDRSZ  # =7824
PacketH_struct_fmt = 'BBHHBBII'
PacketH_struct_len = struct.calcsize(PacketH_struct_fmt)
cint16 = numpy.dtype([('re', '<i2'), ('im', '<i2')])
cint8 = numpy.dtype([('re', '<i1'), ('im', '<i1')])

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
        xrxiyryi16_dtype = numpy.dtype([('xr', '<i2'), ('xi', '<i2'),
                                        ('yr', '<i2'), ('yi', '<i2')])
        xrxiyryi8_dtype = numpy.dtype([('xr', '<i1'), ('xi', '<i1'),
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
        bffmtparams.packetData_dtype = numpy.dtype((xrxiyryi_dtype,
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
    #header['packno'] = ((timestamps*1000000*(160 + is200mhz*40) +oe*512)/1024+blocksequencenumber)/16
    # Get packet format parameters:
    packetData_dtype, cint_smp = bffmtparams(header)

    # Unpack BL data
    alt = 2
    if alt == 1:
        LofarElemSamps = numpy.fromfile(filepointer, dtype=packetData_dtype,
                                        count=1)
    else:
        buf = filepointer.read(packetData_dtype.itemsize)
        LofarElemSamps = numpy.frombuffer(buf, dtype=packetData_dtype)
    xr = (LofarElemSamps[0, :, ]['xr']).squeeze()
    xi = (LofarElemSamps[0, :, ]['xi']).squeeze()
    yr = (LofarElemSamps[0, :, ]['yr']).squeeze()
    yi = (LofarElemSamps[0, :, ]['yi']).squeeze()
    if keepstruct:
        # Keep packet as integer values
        x = numpy.zeros((nrbeamlets, NRTIMS_PACKET), dtype=cint_smp)
        y = numpy.zeros((nrbeamlets, NRTIMS_PACKET), dtype=cint_smp)
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

    # nrPackets = float(os.stat(bfs_filename).st_size - startskip + endpad) / bytesperpacket
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
                    _x = numpy.zeros_like(x)
                    _y = numpy.zeros_like(y)
                    yield _header_, _x, _y
        packetnr += 1
        yield header, x, y
    if padmissing:
        print('Missed packets total:', missed_pkts_tot, packetnr)
    fin.close()


class BFSmeta:
    def __init__(self, bfs_filename, sb=0):
        # Lookup first and last packet in file to set things up:
        with open(bfs_filename, "rb") as fin:
            fin.seek(0, os.SEEK_SET)
            header_first, _x, _y = read_bf_packet(fin)
            fin.seek(0, os.SEEK_END)
            self.nrpkts_actual = fin.tell()//BytesBFPacket
            fin.seek(-BytesBFPacket, os.SEEK_END)
            header_last, _x, _y = read_bf_packet(fin, keepstruct=True)
        self.nrbeamlets = header_first['nrbeamlets']
        self.start_time = header_first['datetime64']
        self.stop_time = header_last['datetime64']
        self.samprate = get_samprate(header_first['is200mhz'])
        self.dur = (self.stop_time - self.start_time).item() /1e9
        self.nrpkts_nominal = round(self.dur * self.samprate
                                    / (NRTIMS_PACKET * FFTSIZE)) + 1
        self.xy_dtype = _x.dtype
        self.lane = header_first['lanenr']

    def __str__(self):
        ret = (f"start: {self.start_time}\n"
               f"stop : {self.stop_time}\n"
               f"dur: {self.dur}\n"
               f"pkts_act: {self.nrpkts_actual}\n"
               f"pkts_nom: {self.nrpkts_nominal}\n"
               f"lane: {self.lane}")
        return ret

def readbfsfile(bfs_filename, savenpy=True):
    """Read BFS file"""
    bfsmeta = BFSmeta(bfs_filename)
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
    pkt_abs_nrs = []
    # Now read in all the packets
    pktnr = 0
    _leftpc_prev = 0
    for header, x, y in next_bfpacket(bfs_filename, True, True):
        if header == '':
            break
        pkt_abs_nr = \
            round((header['datetime64'] - bfsmeta.start_time).item()
                  / (NRTIMS_PACKET * FFTSIZE / bfsmeta.samprate * 1e9))
        pktnr_pe = pktnr*NRTIMS_PACKET + NRTIMS_PACKET
        x_strm[:, NRTIMS_PACKET*pktnr:pktnr_pe] = x
        y_strm[:, NRTIMS_PACKET*pktnr:pktnr_pe] = y
        leftpc = (pkt_abs_nr+1)/bfsmeta.nrpkts_nominal * 100
        if leftpc - _leftpc_prev > 1.0:
            _leftpc_prev = leftpc
            print('\r'+str(round(leftpc))+'% completed', end='')
        pkt_abs_nrs.append(pkt_abs_nr)
        pktnr += 1
    print()
    return x_strm, y_strm, bfsmeta.start_time, pkt_abs_nrs


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
            xy = numpy.zeros((2, header['nrbeamlets'], NRTIMS_PACKET),
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
        exit()
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
            xx = numpy.zeros((nrbeamlets,), dtype=float)
            yy = numpy.zeros((nrbeamlets,), dtype=float)
            xy = numpy.zeros((nrbeamlets,), dtype=complex)
            integrate_samps = (math.floor(integration_req * samprate
                                          / (FFTSIZE * NRTIMS_PACKET))
                               * NRTIMS_PACKET)
            integration = integrate_samps / samprate * FFTSIZE
            print("Requested {}s integration, but will instead use {}s".format(
                integration_req, integration))
            normfac = 1.0 / integrate_samps * (samprate / FFTSIZE)
            initialize = False
        for bi in range(nrbeamlets):
            xx[bi] += numpy.real(numpy.vdot(x[bi, :], x[bi, :]))
            yy[bi] += numpy.real(numpy.vdot(y[bi, :], y[bi, :]))
            xy[bi] += numpy.vdot(x[bi, :], y[bi, :])
            # yx[bi]=numpy.vdot(y[bi,:],x[bi,:])
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
            xx = numpy.zeros((nrbeamlets,), dtype=float)
            yy = numpy.zeros((nrbeamlets,), dtype=float)
            xy = numpy.zeros((nrbeamlets,), dtype=complex)
    foutXX.close()
    foutYY.close()
    foutXY.close()


def convert2npy(bfs_filepath):
    """\
    Convert BFS file to numpy npy file
    """
    readbfsfile(bfs_filepath, savenpy=True)


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
    plt.pcolormesh(numpy.abs(x))
    plt.ylabel('Beamlet [#]')
    plt.colorbar()
    plt.title('abs(X)')

    plt.subplot(2,2,2)
    plt.pcolormesh(numpy.angle(x))
    plt.colorbar()
    plt.title('arg(X)')

    plt.subplot(2,2,3)
    plt.pcolormesh(numpy.abs(y))
    plt.colorbar()
    plt.xlabel('Time bin [#]')
    plt.ylabel('Beamlet [#]')
    plt.title('abs(Y)')

    plt.subplot(2,2,4)
    plt.pcolormesh(numpy.angle(y))
    plt.colorbar()
    plt.xlabel('Time bin [#]')
    plt.title('arg(Y)')

    plt.suptitle('BFS packet lane {} @ {}'.format(header['lanenr'],
                                                  header['datetime'].isoformat()))
    plt.show()


def check_packets(bfs_filename):
    """\
    Check the packets of a BFS file
    """
    packetnr = 0
    missedpackets = 0
    prev_dt =None
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
        if prev_dt:
            print(header['datetime64'] -prev_dt)
        prev_dt = header['datetime64']

    print("Missed packets: ", missedpackets, "/", packetnr, " ",
          100 * missedpackets / float(packetnr), '%')


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

    parser_npy = subparsers.add_parser('npy', help='Convert BFS to numpy npy')
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
