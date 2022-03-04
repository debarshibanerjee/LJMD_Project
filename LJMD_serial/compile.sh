#!/bin/bash
echo "Cleaning old files..."
rm -f ljmd*.x
rm -rf build/

echo "Preparing to compile with CMake..."
mkdir -p build
cmake -S. -B ./build -DCMAKE_EXPORT_COMPILE_COMMANDS=1
cmake --build build/

