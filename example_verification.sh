#!/bin/bash

echo "Running basic make check in examples dir..."
cd examples/
make check
diff argon_108.dat ../reference/argon_108.dat
cd ..
