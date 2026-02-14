# docker build -t radio-cartographer .

FROM ubuntu:latest

#set up file system
RUN mkdir /skynet
WORKDIR /skynet

# Install dependencies
# Install dependencies
RUN apt-get update && apt-get -y install \
    python3 \
    python3-pip \
    cmake \
    gdb \
    curl \
    unzip \
    p7zip-full \
    tar \
    g++ \
    && rm -rf /var/lib/apt/lists/*

## zlib
RUN curl https://zlib.net/zlib131.zip -o /skynet/zlib131.zip
RUN unzip zlib131.zip
RUN rm zlib131.zip
RUN mkdir /skynet/zlib.build
WORKDIR /skynet/zlib.build
RUN cmake ../zlib-1.3.1 -DCMAKE_INSTALL_PREFIX=/skynet/zlib
RUN cmake --build . --config Release
RUN cmake --install .
WORKDIR /skynet


RUN export CXX=g++

## CFITSIO
RUN curl https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/c/cfitsio-4.6.3.tar.gz -o /skynet/cfitsio-4.6.3.tar.gz --retry 10
RUN tar -xf cfitsio-4.6.3.tar.gz
WORKDIR /skynet/cfitsio-4.6.3
# Set flags so configure can find the custom zlib headers and libs
ENV CPPFLAGS="-I/skynet/zlib/include"
ENV LDFLAGS="-L/skynet/zlib/lib"
RUN ./configure --prefix=/usr/local
RUN make
RUN make install
WORKDIR /skynet


## CCFits
RUN curl https://heasarc.gsfc.nasa.gov/FTP/software/fitsio/ccfits/CCfits.tar.gz -o /skynet/CCfits-2.7.tar.gz --retry 10

RUN tar -xf CCfits-2.7.tar.gz
WORKDIR /skynet/CCfits-2.7
RUN ./configure --with-cfitsio=/usr/local/cfitsio --prefix=/usr/local
RUN gmake
RUN make install
WORKDIR /skynet

## FFTW3
RUN curl http://www.fftw.org/fftw-3.3.10.tar.gz -o /skynet/fftw-3.3.10.tar.gz
RUN tar -xf fftw-3.3.10.tar.gz
WORKDIR /skynet/fftw-3.3.10
RUN ./configure
RUN make
RUN make install
WORKDIR /skynet

## Set environment
ENV LD_LIBRARY_PATH=/usr/local/lib
RUN ldconfig

# Setup RC
# Download testing standards (To safe static location)
RUN mkdir -p /skynet/test_standards_static
RUN curl -L "https://www.dropbox.com/scl/fo/1jd5vl7gta57yngrxjlyr/ABl5GBK2wnjpBBS1tLAF3x4?rlkey=y3gv086i7ab06ncn0ydx0ijjb&st=b2wgwi0a&dl=1" -o /skynet/test_standards.zip
RUN 7z x /skynet/test_standards.zip -o/skynet/test_standards_static
RUN rm /skynet/test_standards.zip

## Transfer files
RUN mkdir /skynet/radio-cartographer
COPY . /skynet/radio-cartographer/
RUN pip3 install -r /skynet/radio-cartographer/pyrc/requirements.txt --break-system-packages
# ENV PYTHONPATH="${PYTHONPATH}:/skynet/radio-cartographer/pyrc"


## Build
WORKDIR /skynet/radio-cartographer
RUN cmake .
RUN cmake --build .
WORKDIR /skynet

