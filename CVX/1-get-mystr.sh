#!/bin/bash

rm -f mystr.dat

for out in ../FM/Ni/OUTCAR ../FM/Ce/OUTCAR ../FM/W/OUTCAR 
do
    name=${out%"/OUTCAR"}
    echo $name
    echo -n "$name " >> mystr.dat
    cp $name/CONTCAR.3 $name/CONTCAR
    vasp_out.py -p $name >> mystr.dat
done
