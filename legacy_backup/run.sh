#!/bin/bash


#!/bin/bash

    # The inputs are as follows:
    # arg 1  =  Filename                    The file to be processed
    # arg 2  =  Processing coordinates      These are the coordinates that we process the data in
    # arg 3  =  Flux Channel                Processing channel
    # arg 4  =  Calibration Method          Which calibration data to be used
    # arg 5  =  Background Scale            In beamwidths
    # arg 6  =  RFI Scale                   In beamwidths
    # arg 7  =  Weight Scale                  ...
    # arg 8  =  Time Shifting               1 or 0 corresponds to on/off
    # arg 9  =  Raw                         1 or 0 corresponds to on/off (Generate a raw map of the originial datafile)
    # arg 10 =  LSS                         1 or 0 corresponds to large scale structures (LSS)/small scale structures (SSS)
    # arg 11 =  Photometry on               1 or 0 corresponds to on/off
    # arg 12 =  Aperture Radius             The radius around the source
    # arg 13 =  Annulus Radius              The outer radius for background measurements
    # arg 14 =  Centroid Finding Method     Look for center, brightest, or coordinates
    # arg 15 =  Coordinates                 Location of sources
    # arg 16 =  Trim Size                   Trims the turning edges
    # arg 17 =  M10+ Criteria               1 or 0 corresponds to true/on or false/off

            # Filename                                    Coord       Channel     Cali            Bg   RFI  Wgt  TS  Raw LSS Pho Apt   Ann   Finding     Loc Trm M10

# ./radio-cartographer   Skynet_58592_cyg_a_38879_46208.cyb.fits       equatorial  composite   interpolated    6.0  0.8  0.6  1   0   0   0   1.25  5.0   center      0   0   0


# oj 287
# ./radio-cartographer   Skynet_58709_oj287_8-14_200_41200_49713.cyb.fits       equatorial  composite   interpolated    6.0  0.8  0.6  1   0   0   0   1.25  5.0   center      0   0   0

# ./radio-cartographer   Skynet_58709_oj287_8-14_230_41202_49714.cyb.fits       equatorial  composite   interpolated    6.0  0.8  0.6  1   0   0   0   1.25  5.0   center      0   0   0

          # Filename                                    Coord       Channel     Cali            Bg   RFI   Wgt  TS  Raw LSS Pho Apt   Ann   Finding     Loc Trm M10
./radio-cartographer "40963_oj287_jul_30_1230.fits"                equatorial  left        interpolated    6.0  0.35  0.6  1   0   0   1   1.25  5.0   center      0   0   0
