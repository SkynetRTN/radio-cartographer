import os

def main():
    # the directory containing all images to be processed
    dirr = '../OJ287/'

    # The arguments are as follows:
    # arg 1  =  Filename                    The file to be processed
    # arg 2  =  Processing coordinates      These are the coordinates that we process the data in
    # arg 3  =  Flux Channel                Processing channel
    # arg 4  =  Calibration Method          Which calibration data to be used
    # arg 5  =  Background Scale            In beamwidths
    # arg 6  =  RFI Scale                   In beamwidths
    # arg 7  =  Weight Scale                ...
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

    #        Coord       Channel     Cali            Bg   RFI  Wgt  TS  Raw LSS Pho Apt   Ann   Finding     Loc Trm M10
    args = ' equatorial  left        interpolated    6.0  0.7  0.33 1   0   0   0   1.25  5.0   center      0   0   0'

    # get all files in the directory 'dirr'
    for filename in os.listdir(dirr):

        # only process fits files.
        if filename[-4:] != 'fits':
            continue

        # generate command to be run by bash shell
        cmd = './RDP.out ' + filename + args

        # run the command
        os.system(cmd)

        # rename the result files
        os.rename('SSS_main.txt', filename + '.txt')

if __name__ == '__main__':
    main()
