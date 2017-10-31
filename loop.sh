#!/bin/bash


for b in `seq 0.8 0.01 1.2`; do 
   u1twist -Nx 12 -Ny 12 -Nz 16 -Nt 12 -beta $b -w 0.1 -warm 1500 -trajecs 1000 -meas 10 >> out.twistL12z16
done 
