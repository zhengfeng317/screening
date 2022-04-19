#!/bin/bash

for k in Ef_*dat
do
    ./phdia.py $k
    echo $k
done
