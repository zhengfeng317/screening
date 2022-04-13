#!/bin/bash

main=`pwd`

for k in *_*
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

