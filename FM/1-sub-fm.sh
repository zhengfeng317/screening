#!/bin/bash

for sys in ../str/*.vasp
do
    n1=${sys#"../str/"}
    name=${n1%".vasp"}
    echo $name

    cp -r init-rlx $name
    cp ../NM/$name/CONTCAR.3 $name/POSCAR
    cp ../NM/$name/POTCAR $name/
    cd $name
    sbatch --job-name=fm.$name job.pbs
    cd ..
done
