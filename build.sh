#!/bin/bash
set -e

# Detect number of cores (macOS/Linux)
num_cores=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 1)

# Create separate build directory for Debug builds
mkdir -p build/debug
cd build/debug

# Run CMake for Debug build
cmake -DCMAKE_BUILD_TYPE=Debug ../..

# Parallel compilation
make -j"$num_cores"

# Copy executable to project root
cp radio-cartographer ../..

echo "Done. Executable 'radio-cartographer' (Debug) is ready."
