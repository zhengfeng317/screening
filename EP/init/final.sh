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
DIM = 1 1 1 
EIGENVECTORS = .TRUE.
EOF
    
    phonopy band.conf --pa '1 0 0 0 1 0 0 0 1' --band '0 0 0 0.5 0.5 0.5'
    cd ..
done
