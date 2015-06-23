#!/bin/bash
mkdir build
cd build/
cmake ../
make
./fdtd_1D < ../in > ../out
cd ..
./scriptOctave
