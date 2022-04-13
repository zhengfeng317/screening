#!/bin/bash

for sys in structures/*vasp
do
    n1=${sys#"structures/"}
    name=${n1%".vasp"}
    echo $name

    ele=`sed '6q;d' $sys`
    IFS=' ' read -ra ADDR <<< "$ele"
    A=${ADDR[0]}
    B=${ADDR[1]}
    
    #if [ ! -f pot/${A}.pot ] ; then
    #    echo ${A} not exist
    #fi

    cp -r init $name
    cat pot/${A}.pot pot/${B}.pot > $name/POTCAR
    cp $sys $name/POSCAR
    cd $name
    ./init.sh
    sbatch --job-name=$name scr.pbs
    sbatch --job-name=$name unscr.pbs
    cd ..
done
