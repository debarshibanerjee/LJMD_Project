#!/bin/bash

echo "Cleaning old files..."
rm -f ljmd*.x
rm -rf build/

echo "Preparing to compile with CMake..."
mkdir -p build
cmake -S. -B ./build -D ENABLE_TESTING=yes
cmake --build build/

echo "Compilation Done. Go to examples/ and do 'make check' to run the program for certain pre-defined sample sizes"
echo "Running Tests..."
./test_verification.sh

