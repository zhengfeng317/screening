#!/bin/bash

main=`pwd`

for k in Ba Sr Hf Pr Nd Sm Gd Tb \
    Dy Er Tm Yb Y Sc Ho #Lu 
do
    echo "$k     "

    cd $k/scr 
    tot=`ls POSCAR-* | wc -l`
    for j in POSCAR-*
    do
        num=${j#"POSCAR-"}
        name=$(printf "%03d" $num)
        if ! grep -Fq "F=" $name/vasp.scf
        then
            echo -n " sub-scr-$name "
            cd $name
            cp ~/cont-scr.pbs .
            sbatch --job-name=$k.$num.s cont-scr.pbs
            cd ..
        fi
    done
    finished1=`grep F= */vasp.scf | wc -l`
    cd $main

    cd $k/unscr
    for j in POSCAR-*
    do
        num=${j#"POSCAR-"}
        name=$(printf "%03d" $num)
        if ! grep -Fq "F=" $name/vasp.nscf
        then
            echo -n " sub-unscr-$name "
            cd $name
            cp ~/cont-unscr.pbs .
            sbatch --job-name=$k.$num.u cont-unscr.pbs
            cd ..
        fi
    done
    finished2=`grep F= */vasp.nscf | wc -l`
    cd $main
    
    echo "$finished1/$tot   $finished2/$tot"
done
