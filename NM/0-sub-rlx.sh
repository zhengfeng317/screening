#!/bin/bash

POT=/jet/home/ys3339/works/POT

for sys in ../str/*.vasp
do
    n1=${sys#"../str/"}
    name=${n1%".vasp"}
    echo $name

    ele=`sed '6q;d' $sys`
    IFS=' ' read -ra ADDR <<< "$ele"
    A=${ADDR[0]}
    B=${ADDR[1]}
    C=${ADDR[2]}

    #if [ ! -f pot/${A}.pot ] ; then
    #    echo ${A} not exist
    #fi

    cp -r init-rlx $name
    cat $POT/${A}.pot $POT/${B}.pot $POT/${C}.pot  > $name/POTCAR
    cp $sys $name/POSCAR
    cd $name
    sbatch --job-name=$name job.pbs
    cd ..
done
