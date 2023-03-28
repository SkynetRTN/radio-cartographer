#!/bin/bash

#             # Filename                                    Coord       Channel     Cali            Bg   RFI  Wgt  TS  Raw LSS Pho Apt   Ann   Finding     Loc Trm M10

# ./RDP.out   ../OJ287/oj287_025_39082.cyb.fits             equatorial  left        interpolated    6.0  0.7  0.33 1   0   0   0   1.25  5.0   center      0   0   0
# cp SSS_main.txt oj287_025_39082.cyb.fits.txt

./RDP.out   ../sxu/Skynet_58570_CRAB_38560_45388.cyb.fits   equatorial  left        interpolated    6.0  0.7  0.33 1   0   0   0   1.25  5.0   center      0   0   0
cp SSS_main.txt Skynet_58570_CRAB_38560_45388.cyb.fits.txt
