#!/bin/bash

echo "Cleaning old files..."
rm -f ljmd*.x
rm -rf build/

echo "Preparing to compile with CMake..."
mkdir -p build
cmake -S. -B ./build -D ENABLE_TESTING=yes
cmake --build build/

./test_verification.sh

./example_verification.sh


