#!/usr/bin/env python

import sys
import numpy
import datetime
import argparse
import ilisa.antennameta.antennafieldlib as antennafieldlib
import ilisa.antennameta.calibrationtables as calibrationtables
import ilisa.observations.dataIO as dataIO
import ilisa.observations.imaging as imaging
import ilisa.observations.modeparms


def main(accrunfolder, desiredsrc, use_autocorr=False):
    """Convert the ACC folder into a hdf or numpy file.

    The output file will be created in the current working path.

    Parameters
    ----------
    accrunfolder : str
        Path to the folder that contains a set of temporally contiguous ACC files.

    desiredsrc : str
        The source to point at.
    """
    accffobj = dataIO.CVCfiles(accrunfolder)
    obsfolderinfo = accffobj.getobsfolderinfo()
    obsdate = obsfolderinfo['datetime']
    rcumode = obsfolderinfo['rcumode']
    calsrc = obsfolderinfo['calsrc']
    # dur = obsfolderinfo['duration']
    stnid = obsfolderinfo['stnid']
    if calsrc is None and desiredsrc is None:
        raise ValueError("No calibration source specified")
    elif int(rcumode) == 3 and desiredsrc is not None:
        calsrc = desiredsrc
    elif desiredsrc != calsrc:
        print("Warning: calibration source {} is not the desired one {}."
              .format(calsrc, desiredsrc))
    pntstr = ilisa.observations.modeparms.stdPointings(calsrc)
    pointing = pntstr.split(',')
    bandarr = ilisa.observations.modeparms.rcumode2antset(rcumode).split("_")[0]
    stn_pos, stn_rot, antpos, stn_intile_pos = antennafieldlib.getArrayBandParams(stnid,
                                                                                  bandarr)
    freqs = ilisa.observations.modeparms.rcumode2sbfreqs(rcumode)
    accfolderdatestr = obsdate.strftime("%Y%m%d") + "T120000"  # FIX put appropriate time
    accfiles = accffobj.filenames
    nraccfiles = len(accfiles)
    (caltab, ctheader) = calibrationtables.getcaltab(rcumode, stnid, accfolderdatestr)

    sb = None
    if sb is None:
        (bstXX, bstXY, bstYY) = (
            numpy.zeros((nraccfiles, ilisa.observations.modeparms.TotNrOfsb)),
            numpy.zeros((nraccfiles, ilisa.observations.modeparms.TotNrOfsb),
                        dtype=complex),
            numpy.zeros((nraccfiles, ilisa.observations.modeparms.TotNrOfsb)))
    else:
        (bstXX, bstXY, bstYY) = (numpy.zeros((nraccfiles, )),
                                 numpy.zeros((nraccfiles, ), dtype=complex),
                                 numpy.zeros((nraccfiles, )))
    filestarttimes = numpy.zeros((nraccfiles,), dtype='datetime64[s]')

    for tidx, ACCfile in enumerate(accfiles):
        accunc = accffobj.getdata(tidx)
        sbobstimes = accffobj.samptimes[tidx]
        filestarttimes[tidx] = numpy.datetime64(sbobstimes[0])
        if tidx == 0:
            calrunstarttime = sbobstimes[0]
        else:
            ftd = numpy.asscalar((filestarttimes[tidx]-filestarttimes[tidx-1])
                                 .astype('int'))
            if ftd != 519:
                print("Time delta between ACC files is {}s.".format(ftd))
                raise RuntimeError("""Time between ACC files nr {}-{} not the \
                                   nominal 519s.""".format(tidx, tidx-1))
        # Calibrate
        acc = calibrationtables.calibrateACCwithtab(accunc, caltab)
        accpol = dataIO.cvc2cvpol(acc)
        sys.stdout.write('({}/{})\n'.format(tidx+1, nraccfiles))
        if sb is None:
            (bstXX[tidx, :], bstXY[tidx, :], bstYY[tidx, :]
             ) = imaging.accpol2bst(accpol, sbobstimes, freqs, stn_pos, antpos, pointing,
                                    use_autocorr=use_autocorr)
        else:
            (bstXX[tidx], bstXY[tidx], bstYY[tidx]
             ) = imaging.xst2bst(accpol[:, :, sb, ...].squeeze(), sbobstimes[sb],
                                 freqs[sb], stn_pos, antpos, pointing)
    calrunendtime = sbobstimes[-1]
    accsbsampleduration = datetime.timedelta(seconds=1)
    calrunduration = calrunendtime - calrunstarttime + accsbsampleduration
    acc2bstname = dataIO.saveacc2bst((bstXX, bstXY, bstYY), filestarttimes,
                                     calrunstarttime, calrunduration, rcumode,
                                     calsrc, ctheader['Calibration'], stnid,
                                     use_autocorr)
    print("Created file: {}".format(acc2bstname))


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("accrunfolder")
    parser.add_argument("desiredsrc", type=str, nargs='?', default=None)
    args = parser.parse_args()
    main(args.accrunfolder, args.desiredsrc)
