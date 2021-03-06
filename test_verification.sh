#!/bin/bash

echo "Compilation Done. Go to examples/ and do 'make check' to run the program for certain pre-defined sample sizes"
echo "Running Tests..."
cd build/tests
export OMP_DYNAMIC=FALSE
ctest -V
cd ../..
