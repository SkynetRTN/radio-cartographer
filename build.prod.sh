#!/bin/bash
set -e

# Detect number of cores (macOS/Linux)
num_cores=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 1)
echo "Building in Release mode with native optimizations using $num_cores cores..."

# Create separate build directory for Release builds
mkdir -p build/release
cd build/release

# Run CMake for Release build with optimizations
cmake -DCMAKE_BUILD_TYPE=Release -DENABLE_NATIVE_OPTIMIZATION=ON ../..

# Parallel compilation
make -j"$num_cores"

# Copy executable to project root
cp radio-cartographer ../..

echo "Done. Executable 'radio-cartographer' (Release, Optimized) is ready."
