#!/bin/bash
set -e

# Create build directory
mkdir -p build
cd build

# Run CMake and Make
cmake ..
make

# Copy executable to root
cp radio-cartographer ..

echo "Build successful. Executable 'radio-cartographer' is ready."
