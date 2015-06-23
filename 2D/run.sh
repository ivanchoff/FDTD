#!/bin/bash
let "N  = $1 + 1"
mkdir build
cd build
cmake ../
make
./fdtd_2D < ../in2d > ../out
cd ..
echo  $1
echo  $N
head -n $1 out > out0
tail -n $N out > outT
