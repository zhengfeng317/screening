#!/bin/bash


for k in   Mn V 
do
    ./pd_example.py $k "B" "C"
    cp mp_stable_ene.dat ${k}_stable.dat
done
