# Single-Dish Radio Mapping Algorithm

This software was developed for use with Skynet's 20-m radio telescope located in Green Bank, West
Virginia. It has been integrated into the Skynet and automatically processes radio observations.

For a full description of the algorithm, see our systems paper:
[SKYNET ALGORITHM FOR SINGLE-DISH RADIO MAPPING](https://arxiv.org/abs/1808.06128).

# Local Installation

The Radio Development Project (RDP) will run on both Windows and Linux. It's possible that it will run
on MacOS, but we have never attempted that internally.

Most of these instructions are condensed from the install docs for each library. If any of these steps
fail, it is recommended to go to the following links for additional details and alternative methods:
[CMake](https://cmake.org/),
[ZLIB](https://zlib.net/),
[CFITSIO](https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/README.win),
[CCFITS](https://heasarc.gsfc.nasa.gov/fitsio/CCfits/html/installation.html),
[FFTW](http://www.fftw.org/install/windows.html).

## Windows

### 1. Install CMake

CMake can be freely downloaded from http://www.cmake.org. It is recommended that you choose the
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

If anything went wrong, try checking the
CFITSIO [README](https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/README.win).
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

You can generate a Visual Studio solution with CMake. From a Developer Command Prompt, create a
build directory (for example, `build`) parallel to the source tree and run:

```
cmake -S . -B build -G "Visual Studio 17 2022" -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release
```

If your dependencies are installed outside the default search paths, add them to
`CMAKE_PREFIX_PATH`, e.g. `-DCMAKE_PREFIX_PATH="C:/path/to/cfitsio;C:/path/to/fftw"`.
The generated solution builds the `radio-cartographer` executable and the
`radio_cartographer` library target.

### Without Visual Studio

Use the CMake-based build instead of the old `compileStandard.sh` script. A typical build looks like:

```
mkdir -p build
cmake -S . -B build \
  -DCMAKE_BUILD_TYPE=Release \
  -DCMAKE_PREFIX_PATH="/path/to/fftw;/path/to/cfitsio" \
  -DBUILD_PYTHON_BINDINGS=OFF
cmake --build build
```

Key notes:

- The project depends on FFTW3, CFITSIO, CCFITS, pthreads/Win32 threads, and nlohmann_json
  (automatically fetched if not already available).
- Set `BUILD_PYTHON_BINDINGS=ON` to build the shared library with position-independent code and
  default symbol hiding for Python wrappers.
- Run `cmake --install build --prefix /your/prefix` to install the executable, library, headers,
  and CMake package config under the chosen prefix (`/usr/local` by default).

## Running

RDP takes an ever growing number of parameters via the command line.

### With Visual Studio

https://stackoverflow.com/questions/298708/debugging-with-command-line-parameters-in-visual-studio

### Without Visual Studio

We have written a run script called `run2.sh` that describes every input parameter and provides an
interface to modify them. Simply specify the path to the input `.FITS` file, modify your inputs,
and then run the script.

# Docker + CLion + CMake all in one dev solution
_added by Reed Fu, May 13 2024_
## Install CLion

Download and install CLion from [here].(https://www.jetbrains.com/clion/download/)
It's free for students or researchers at academic institutions.

Open the project in CLion.

## Docker

Install Docker Desktop from [here].(https://www.docker.com/products/docker-desktop)

In the terminal (either in CLion or your terminal of choice), under the project root directory, run:

```bash
docker build -t radio-cartographer .
```

This will build the docker image, and name it as `radio-cartographer`.

You could run a container from the image interactively inside docker app by running:

```bash
docker run -d -it --name rc -v "$(pwd)":/skynet/radio-cartographer radio-cartographer:latest
```

However, it's way easier to run and debug the code in CLion.

## CLion

First we want to set up the remote toolchain in CLion.

Open settings, then go to `Build, Execution, Deployment` -> `Toolchains`.

Click the `+` button, and select `Docker`.
Select the image as `radio-cartographer:latest`, and open the container setting popup.
We don't need to configure port bindings,
but we need to bind your local project directory to the container's `/skynet/radio-cartographer` directory.

Click `OK` to save the settings.

In the top right corner, you can select the tool chain and the run configuration.

There should be a `debug` configuration already set up for you, click the hammer icon to build the project.
Then click the green triangle to run the project. You can also set break points and debug the code.

### Modify Run Configuration

The working directory is set up to be `./testing/tmp`, all output txt will be ignored by git.

You can modify which file to run by changing the arguments in the configuration.
Just be aware that since everything is running inside the container,
it's the easiest to copy your file to /testing, and only modify the filename in the arguments.

You can also modify all the other RC arguments. Good luck figure out what each of them does. :)

## Antigravity Setup

To configure Antigravity to use the Docker environment, use the provided `docker-shell.sh` script. This script acts as a bridge between the local filesystem and the Docker container's build environment.

1. **Ensure the script is executable:**
   ```bash
   chmod +x docker-shell.sh
   ```

2. **Configure Antigravity Terminal Profile:**
   Since Antigravity is based on VSCode, you can configure the integrated terminal to use the Docker environment automatically.

   1. Open the Command Palette (`Cmd+Shift+P` on Mac, `Ctrl+Shift+P` on Windows/Linux).
   2. Type `Preferences: Open Settings (JSON)` and select it to open your `settings.json` file.
   3. Add the following configuration to the file.

   *Note: Change `osx` to `linux` or `windows` in the property keys depending on your operating system.*

   ```json
   "terminal.integrated.profiles.osx": {
     "Docker Sandbox": {
       "path": "${workspaceFolder}/docker-shell.sh",
       "icon": "container"
     }
   },
   "terminal.integrated.defaultProfile.osx": "Docker Sandbox"
   ```
