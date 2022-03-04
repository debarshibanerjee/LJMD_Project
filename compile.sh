#!/bin/bash
if [[ $# -eq 0 ]] ; then
	echo 'Run with one argument: 1 or 0. If 1, OpenMP is enabled. If 0, OpenMP is disabled.'
	echo 'USAGE: ./compile.sh 1'
	echo 'For OpenMP threads, do: export OMP_NUM_THREADS=#threads'
	exit 0
fi

echo "Cleaning old files..."
rm -f ljmd*.x
rm -rf build/

echo "Preparing to compile with CMake..."
mkdir -p build
cmake -S. -B ./build -D ENABLE_TESTING=no -DCMAKE_EXPORT_COMPILE_COMMANDS=1 -DENABLE_OPENMP=$1
cmake --build build/

export OMP_DYNAMIC=FALSE
export OMP_NUM_THREADS=2
#./test_verification.sh

export OMP_DYNAMIC=FALSE
export OMP_NUM_THREADS=2
./result_verification.sh

