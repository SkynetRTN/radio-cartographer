#!/bin/bash

# Working directory
DIR=$(pwd)
# "$( cd "$( dirname "${BASH_SOURCE[0]}" )" >/dev/null 2>&1 && pwd )"
echo pwd
# Files to be processed in the directory
# FILES=testing/good-test-file.cyb.fits        # arg 1  =  Filename                    The file(s) to be processed
# *_2.txt this was good when we tested it using .txt

# The following are the processing parameters
Coord=equatorial    # arg 2  =  Processing coordinates      These are the coordinates that we process the data in
Channel=left        # arg 3  =  Flux Channel                Processing channel
Cali=interpolated   # arg 4  =  Calibration Method          Which calibration data to be used
Bg=6.0              # arg 5  =  Background Scale            In beamwidths
RFI=0.05            # arg 6  =  RFI Scale                   In beamwidths
Wgt=0.33            # arg 7  =  Weight Scale                ...
TS=0                # arg 8  =  Time Shifting               1 or 0 corresponds to on/off
Raw=1               # arg 9  =  Raw                         1 or 0 corresponds to on/off (Generate a raw map of the originial datafile)
LSS=0               # arg 10 =  LSS                         1 or 0 corresponds to large scale structures (LSS)/small scale structures (SSS)
Phot=0              # arg 11 =  Photometry on               1 or 0 corresponds to on/off
Apt=1.25            # arg 12 =  Aperture Radius             The radius around the source
Ann=3.0             # arg 13 =  Annulus Radius              The outer radius for background measurements
Fd=center           # arg 14 =  Centroid Finding Method     Look for center, brightest, or coordinates
Loc=0               # arg 15 =  Coordinates                 Location of sources
Trm=0               # arg 16 =  Trim Size                   Trims the turning edges
M10=0               # arg 17 =  M10+ Criteria               1 or 0 corresponds to true/on or false/off

# The location and filename for the log file
Log="${DIR%/}/run_sh.log"

# This function sets the directory name for each processed image
set_output_directory () {
    filename="${1##*/}"         # Remove the path from the beginning
    filename="${filename%.*}"   # Remove the extension from the end
    echo outdir
    echo -n "${DIR%/}/${filename}_phot"
}

# This function outputs processing info, and write the info to log file as well
info () {
    echo -e "$1" | expand -t 4 | tee -a "$Log"
}

# Output processing info, and write the info to log file as well
info "$(date)"
info "Working directory: ${DIR%/}"
info "Parameters:"
info "\tCoord\t= $Coord"
info "\tChannel\t= $Channel"
info "\tCali\t= $Cali"
info "\tBg\t\t= $Bg"
info "\tRFT\t\t= $RFI"
info "\tWgt\t\t= $Wgt"
info "\tTS\t\t= $TS"
info "\tRaw\t\t= $Raw"
info "\tLSS\t\t= $LSS"
info "\tPhot\t= $Phot"
info "\tApt\t\t= $Apt"
info "\tAnn\t\t= $Ann"
info "\tFd\t\t= $Fd"
info "\tLoc\t\t= $Loc"
info "\tTrm\t\t= $Trm"
info "\tM10\t\t= $M10"

# The ${FILES} should not be quote enclosed because we need to treat the expansion as an array
#   of multiple items (instead of a string with spaces in it)
# for file in "${DIR%/}/"${FILES}
shopt -s nullglob
files_array=(testing/good-test-file.cyb.fits)

for ((i=0; i<${#files_array[@]}; i++)); do
  files_array[$i]="${DIR}/${files_array[$i]}"
done
# echo ${files_array[@]}

    file=${files_array[0]}
    echo $file
    info "Processing $file"

    ./radio-cartographer $file $Coord $Channel $Cali $Bg $RFI $Wgt $TS $Raw $LSS $Phot $Apt $Ann $Fd $Loc $Trm $M10 ${files_array[@]}

    if [ $? -eq 0 ]
    then
        # If the image was successfully processed, output info and copy processing results
        Out=$(set_output_directory "$file")

        info "Copying to ${Out}"

        mkdir "$Out"
        cp *.txt "${Out}/"

        info "$file has been processed successfully"
    else
        # If the image failed to process, output info only
        info "$file failed to process and is not finished" 
    fi

    info "---------------"
    
# done

# ------- The old-fashioned way --------

#             # Filename                                      Coord       Channel     Cali            Bg   RFI  Wgt  TS  Raw LSS Pho Apt   Ann   Finding     Loc Trm M10

# ./radio-cartographer   ../OJ287/oj287_025_39082.cyb.fits               equatorial  left        interpolated    6.0  0.7  0.33 1   0   0   0   1.25  5.0   center      0   0   0
# cp SSS_main.txt oj287_025_39082.cyb.fits.txt

# ./radio-cartographer   ../OJ287/oj287_025_39081.cyb.fits               equatorial  left        interpolated    6.0  0.7  0.33 1   0   0   0   1.25  5.0   center      0   0   0
# cp SSS_main.txt oj287_025_39081.cyb.fits.txt
