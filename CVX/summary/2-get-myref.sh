#!/bin/bash

rm -f ref.dat

for out in ../*/OUTCAR
do
    name=${out%"/OUTCAR"}
    echo $name
    echo -n "$name " >> ref.dat
    cp $name/CONTCAR.3 $name/CONTCAR
    vasp_out.py -p $name >> ref.dat
done
