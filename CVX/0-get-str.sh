#!/bin/bash


for k in "Ce" "Li" "Ni"
do
    ./pd_example.py $k "B" "C"
    ./get-mp-struc.py
    cp mp_stable_ene.dat ${k}_stable.dat
done
