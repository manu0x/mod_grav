#!/bin/sh
g++ ../delta_gen_friedmann_Cs.cpp -lm
./a.out bimetric 0.3 0.0
./a.out bimetric 0.3 0.5
./a.out bimetric 0.3 1.0
./a.out bimetric 0.3 3.0
./a.out bimetric 0.3 6.0



./a.out lcdm 0.3 0.0
./a.out lcdm 0.4 0.0
./a.out lcdm 0.7 0.0
./a.out lcdm 1.0 0.0