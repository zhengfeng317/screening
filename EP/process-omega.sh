#!/bin/bash

rm -f lambda.dat
echo "#sys lmd_max lmd_sum" >> lambda.dat

for doc in La W Cu Li Al Mg K V Ca Ba Sr Hf Pr Nd Sm Gd Tb \
    Dy Er Tm Yb Y Sc Ho 
do

   echo -n "$doc " >> lambda.dat
   cd $doc
   echo $doc
   ana-lmd-phonopy.py
   ana-lmd-plot.py "$doc" >> ../lambda.dat
   mv lmd.png $doc.png
   cd ..
done
