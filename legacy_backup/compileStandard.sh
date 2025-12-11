#!/bin/bash

# Build the project using CMake
mkdir -p build
cd build
cmake ..
make

# Copy the executable to the root directory for backward compatibility with scripts
if [ -f "radio-cartographer" ]; then
    cp radio-cartographer ..
    echo "Build successful. Executable 'radio-cartographer' created."
else
    echo "Build failed."
    exit 1
fi