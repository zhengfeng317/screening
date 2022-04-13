#!/usr/bin/env python3

RE=['La', 'Fe', 'Mn', 'Ni', 'Co', 'W', 'Cu', 'Li', 
    'Al', 'Mg', 'K', 'Cr', 'V', 'Ca', 'Ba', 'Sr', 
    'Hf', 'Pr', 'Nd', 'Sm', 'Gd', 'Tb', 'Dy', 'Er',
    'Tm', 'Yb', 'Y', 'Sc', 'Ho', 'Lu', 'Ce' ]

natom=20

for rr in RE:
    fin=open("template","r")
    fout=open(rr+".vasp","w+")
    print(rr, file=fout)
    fin.readline()
    for k in range(4):
        print(fin.readline().strip("\n"),file=fout)
    print(rr,"B C",file=fout)
    fin.readline()
    for k in range(natom+2):
        print(fin.readline().strip("\n"),file=fout)
    fin.close()
    fout.close()
