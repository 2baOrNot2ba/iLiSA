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
NRTIMS_PACKET = 16  # Number of sample times in packet
FFTSIZE = 1024
NRPOLS = 2  # Number of polarization channels
NRCMPLX = 2  # Every data sample is complex i.e. real & imag so size is doubled
PACKETHDRSZ = 16  # Size of packet header in bytes
BytesBFPacket = NRCMPLX * NRPOLS * 2 * NRTIMS_PACKET * 61 + PACKETHDRSZ
PacketH_struct_fmt = 'BBHHBBII'
PacketH_struct_len = struct.calcsize(PacketH_struct_fmt)
cint16 = numpy.dtype([('r', '<i2'), ('i', '<i2')])
cint8 = numpy.dtype([('r', '<i1'), ('i', '<i1')])

# PCAP file format data:
BytesPcapHeader = 42  # 42
BytesPcapFooter = 16
BytesPcapFileHeader = 40  # 40
BytesPerPcapPacket = BytesPcapHeader + BytesBFPacket + BytesPcapFooter


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
    header['datetime'] = header['datetime']+datetime.timedelta(microseconds=int(microsecs))
    header['nanosecs'] = round(_microsecfrac * 1000)
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
        x[:, :]['r'] = xr
        x[:, :]['i'] = xi
        y[:, :]['r'] = yr
        y[:, :]['i'] = yi
    else:
        x = xr + 1j * xi
        y = yr + 1j * yi

    pcapfooter = None
    if pcapfile:
        pcapfooter = filepointer.read(BytesPcapFooter)

    return header, x, y


# def plotPacket(sbdata):
#    im = plt.imshow(numpy.abs(sbdata))
#    plt.show()


def next_bfpacket(bfs_filename, keepstruct=True):
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
    packetNr = 0
    fin.seek(startskip)
    while not EOF:
        # print(str(packetNr)+" / "+str(nrPackets))
        header, x, y = read_bf_packet(fin, keepstruct, pcapfile)
        if header == '':
            EOF = True
            break
        packetNr += 1
        yield header, x, y
    fin.close()


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
    bfs_filenamebase = os.path.basename(bfs_filepath)
    nrpackets = os.path.getsize(bfs_filepath)/BytesBFPacket
    foutXX = open(bfs_filenamebase + ".XX.bst", "wb")
    foutYY = open(bfs_filenamebase + ".YY.bst", "wb")
    foutXY = open(bfs_filenamebase + ".XY.bst", "wb")

    totnrtimsamp = 0
    initialize = True
    for header, x, y in next_bfpacket(bfs_filepath, False):
        if initialize:
            # Set packet format parameters using first packet
            #    Get samprate
            if header['is200mhz'] == 1:
                samprate = 200.0e6  # sample rate per second
            else:
                samprate = 160.0e6  # sample rate per second
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
    for header, x, y in next_bfpacket(bfs_filename, False):
        seq = header['blocksequencenumber']
        if header['error'] == 1:
            print("Error in packet")
        if packetnr != 0:
            seqdif = seq - seqprev
        else:
            seqdif = 16
        if seqdif != 16 and seqdif != -195297 and seqdif != -195296:
            missedpackets += 1
            if False:
                print(seqdif)
        seqprev = seq
        packetnr += 1
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

    args = cmdln_prsr.parse_args()

    if args.bfs_filename.endswith('zst'):
        print("Please uncompress the BFS file first")
        sys.exit()

    if args.func == print_packet:
        packetnr =  args.packetnr
        while True:
            header, x, y = get_packet_h_x_y_fromfile(args.bfs_filename,
                                                     packetnr)
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


if __name__ == "__main__":
    main_cli()
