#!/bin/bash
set -e

# Create build directory
mkdir -p build
cd build

# Run CMake and Make with Release build type
cmake -DCMAKE_BUILD_TYPE=Release ..
make

# Copy executable to root
cp radio-cartographer ..

echo "Production build successful. Executable 'radio-cartographer' is ready."
