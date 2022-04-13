#!/bin/bash

for name in La W Cu Li Al Mg K V Ca Ba Sr Hf Pr Nd Sm Gd Tb \
    Dy Er Tm Yb Y Sc Ho Lu
do
    echo $name

    cp -r init $name
    cp ../NM/$name/CONTCAR.3 $name/POSCAR 
    cp ../NM/$name/POTCAR $name/
    cd $name
    ./init.sh
    sbatch --job-name=s.$name scr.pbs
    sbatch --job-name=u.$name unscr.pbs
    cd ..
done
