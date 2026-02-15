#!/bin/bash
set -e

PROJECT_ROOT="$(pwd)"
BUILD_DIR="build/debug"

# Argument parsing
for arg in "$@"
do
    case $arg in
        --docker)
        BUILD_DIR="/skynet/tmp/build/debug"
        shift
        ;;
    esac
done

# Detect number of cores (macOS/Linux)
num_cores=$(sysctl -n hw.ncpu 2>/dev/null || nproc 2>/dev/null || echo 1)

# Create build directory
mkdir -p "$BUILD_DIR"
cd "$BUILD_DIR"

# Run CMake for Debug build
# We need to point back to PROJECT_ROOT.
cmake -DCMAKE_BUILD_TYPE=Debug "$PROJECT_ROOT"

# Parallel compilation
make -j"$num_cores"

# Copy executable to project root
cp radio-cartographer "$PROJECT_ROOT"

echo "Done. Executable 'radio-cartographer' (Debug) is ready."
