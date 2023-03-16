# Single-Dish Radio Mapping Algorithm
This software was developed for use with Skynet's 20-m radio telescope located in Green Bank, West 
Virginia. It has been integrated into the Skynet and automatically processes radio observations. 

For a full description of the algorithm, see our systems paper: 
[SKYNET ALGORITHM FOR SINGLE-DISH RADIO MAPPING](https://arxiv.org/abs/1808.06128).

# Installation
The Radio Development Project (RDP) will run on both Windows and Linux. It's possible that it will run 
on MacOS, but we have never attempted that internally.

Most of these instructions are condensed from the install docs for each library. If any of these steps
fail, I recommend going to the following links for additional details and alternative methods:
[CMake](https://cmake.org/), 
[ZLIB](https://zlib.net/), 
[CFITSIO](https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/README.win), 
[CCFITS](https://heasarc.gsfc.nasa.gov/fitsio/CCfits/html/installation.html),
[FFTW](http://www.fftw.org/install/windows.html).

## Windows

### 1. Install CMake
CMake cane be freely downloaded from http://www.cmake.org. It is recommended that you choose the 
"Add CMake to the system PATH for current user" option during the installation setup process for 
convenience when running CMake later on.

### 2. Install ZLIB
ZLIB can be freely downloaded from https://zlib.net. Once downloaded, unpack it, then build and 
install it from a parallel directory, for example:

```
mkdir zlib.build
cd zlib.build
cmake ..\zlib-1.2.11 -DCMAKE_INSTALL_PREFIX=c:\Users\myname\zlib
cmake --build . --config Release
cmake --install .
```
The cmake comands below will use the path "c:\Users\myname\zlib" as an example for the installed 
zlib location.

### 3. Install CFITSIO
CFITSIO can be freely downloaded from http://heasarc.gsfc.nasa.gov/fitsio. CFITSIO is a dependency for
another library (CCFITS) that we use to read in the data files.

Unzip the CFITSIO .zip file. This will create a new `\cfitsio-4.2.0` subdirectory that contains the 
source code and documentation files. Version `4.2.0` is the most recent release as of writing, but 
future releases should work OK.

Open the Visual Studio Command Prompt window and `cd` into the parent directory that is one level
above the directory containing the CFITSIO source files that was created in the previous step.

Create a new subdirectory, and cd into it with the following commands:
```
mkdir cfitsio.build
cd cfitsio.build
```

Create the CMake files necessary for the build by running (replace the path with the location to ZLIB):
```
cmake ..\cfitsio-4.0.0 -DCMAKE_PREFIX_PATH=c:\Users\myname\zlib
```

Execute the following command to build the CFITSIO library:

```
cmake --build . --config Release
cmake install .
```

 The "." argument following "--build" here tells CMake to build the files in the current directory 
 (i.e., in "cfitsio.build"). If this process completes successfully, you should find the CFITSIO 
 library files in the "cfitsio.build\Release" subdirectory.

If anything went wrong, try checking the CFITSIO [README](https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/README.win).
for additional details on this process.
 
### 4. Install CCFits
CCFITS can be freely downloaded from https://heasarc.gsfc.nasa.gov/fitsio/CCfits/. After unzipping and 
untarring the CCfits source code tarball, the source code will appear in a new `\CCfits` subdirectory.

Open the Visual Studio Developer Command Prompt window and create a directory named `CCfits.build` 
parallel to this CCfits source code directory. Run the following command:

```
mkdir CCfits.build cd CCfits.build
```

Create the CMake files necessary for the build, bukd, then install by runnning:
```
cmake.exe -G"NMake Makefiles" -DCMAKE_PREFIX_PATH=C:\path\to\your\CFITSIO ..\CCfits
nmake
namke install
```

### 5. Install Fastest Fourier Transform in the West (FFTW3)
FFTW3 can be freely downloaded from http://www.fftw.org/install/windows.html. The developers provide
pre-compiled Windows DLLs for you to use. You should be able to call them from any compiler. In order 
to link to them from Visual C++, you will need to create .lib "import libraries" using the lib.exe 
program included with VC++. Run:
```
lib /def:libfftw3-3.def
lib /def:libfftw3f-3.def
lib /def:libfftw3l-3.def
```
On Visual Studio 2008 in 64-bit mode, and possibly in other cases, you may need to specify the machine
explicitly:
```
lib /machine:x64 /def:libfftw3l-3.def
```

## Building
### With Visual Studio
Needs updating

### Without Visual Studio
We have a (crude) compile script called `compileStandard.sh` that will compile the project. I often
dream of writing a `makefile`, but until that day, we must recompile the entire project with the script
even if we make a one-line change to a single file.

## Running
RDP takes an ever growing number of parameters via the command line.
### With Visual Studio
https://stackoverflow.com/questions/298708/debugging-with-command-line-parameters-in-visual-studio

### Without Visual Studio
We have written a run script called `run2.sh` that describes every input parameter and provides an
interface to modify them. Simply specify the path to the input `.FITS` file, modify your inputs,
and then run the script.