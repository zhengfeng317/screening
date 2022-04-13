#!/bin/bash

main=`pwd`

for k in La W Cu Li Al Mg K V Ca Ba Sr Hf Pr Nd Sm Gd Tb \
    Dy Er Tm Yb Y Sc Ho Lu 
do
    echo -n "$k     "

    cd $k/scr 
    tot=`ls POSCAR-* | wc -l`
    finished1=`grep F= */vasp.scf | wc -l`
    cd $main

    cd $k/unscr
    finished2=`grep F= */vasp.nscf | wc -l`
    cd $main
    
    echo "$finished1/$tot   $finished2/$tot"
done

