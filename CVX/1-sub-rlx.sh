#!/bin/bash

POT="/jet/home/ys3339/works/POT"
INIT="../FM/init-rlx"

for sys in POSCAR.* 
do
    n1=${sys#"POSCAR.mp-"}
    name=${n1} #%".vasp"}
    echo $name

    ele=`sed '6q;d' $sys`
    IFS=' ' read -ra ADDR <<< "$ele"
    A=${ADDR[0]}
    B=${ADDR[1]}
    C=${ADDR[2]}

    #if [ ! -f pot/${A}.pot ] ; then
    #    echo ${A} not exist
    #fi

    cp -r $INIT $name
    
    for ee in ${ADDR[@]} 
    do
         cat $POT/${ee}.pot >> $name/POTCAR
    done
    #cat pot/${A}.pot pot/${B}.pot pot/${C}.pot  > $name/POTCAR
    cp $sys $name/POSCAR
    cd $name
    sbatch --job-name=$name job.pbs
    cd ..
done
