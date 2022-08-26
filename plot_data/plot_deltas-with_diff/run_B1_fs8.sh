#!/bin/sh
g++ ../../code/delta_gen_friedmann_Cs_woD.cpp -lm
./a.out bimetric 0.3 0.0
./a.out bimetric 0.3 0.2
./a.out bimetric 0.3 0.5
./a.out bimetric 0.3 1.0
./a.out bimetric 0.3 2.0

./a.out bimetric 0.3 -0.2
./a.out bimetric 0.3 -0.5
./a.out bimetric 0.3 -1.0
./a.out bimetric 0.3 -2.0
#./a.out bimetric 0.3 3.0
#./a.out bimetric 0.3 6.0


