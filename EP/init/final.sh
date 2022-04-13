#!/bin/bash

for sys in scr unscr
do
    cd $sys

    ndoc=`ls POSCAR-* |wc -l`
    
    files=""
    for ((k=1; k<=$ndoc ; k++))
    do
        f=`printf "%03d" $k`
        files+="$f/vasprun.xml "
    done
    
    echo $files
    
    phonopy -f $files # {001..002}/vasprun.xml
    
    cat >band.conf<<EOF
    ATOM_NAME = Mg B 
    DIM = 1 1 1 
    BAND = AUTO
    EIGENVECTORS = .TRUE.
EOF
    
    phonopy -p band.conf

    cd ..
done
