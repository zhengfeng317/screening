#!/bin/bash

for name in Ni Ce
do
    echo $name

    cp -r init-rlx $name
    cp ../NM/$name/CONTCAR.3 $name/POSCAR
    cp ../NM/$name/POTCAR $name/
    cd $name
    sbatch --job-name=fm.$name job.pbs
    cd ..
done
