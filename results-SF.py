#!/usr/bin/env python3

import sys
import glob
import numpy as np
import argparse


def get_EM(infile, n_mag):
    J2eV= 1.602176634e-19
    fin=open(infile,"r")
    F=[];E=[]; press=[]; stress=[]; V=[]; Mag='NA'
    natom=0; vtag=0
    atom_mag=[0 for k in range(n_mag)]
    for line in fin:
        if("in kB" in line):
            ll=line.split()
            stress.append(ll[2]+" "+ll[3]+" "+ll[4]+" "+ll[5]+" "+ll[6]+" "+ll[7])
            try:
                press.append((float(ll[2])+float(ll[3])+float(ll[4]))/3)
            except:
                press.append(-99999)
        if(("volume of cell" in line) and vtag==0):
            vtag = 1
        elif(("volume of cell" in line) and vtag==1):
            vol_vasp=float(line.split()[4])
            fin.readline()
            vec=[]
            for k in range(3):
                ll=fin.readline().split()
                fl=[float(ll[s]) for s in range(3) ]
                vec.append(fl)
            vec=np.array(vec)
            vol=np.dot(np.cross(vec[0],vec[1]),vec[2])
            if(abs(vol-vol_vasp)>0.1): print("wrong volume!")
            V.append(vol)
        if("NIONS" in line):
            ll=line.split()
            natom=int(ll[len(ll)-1])
        if("energy  without entropy=" in line):
            ll=line.split()
            E.append(float(ll[-1])) #"sigma->0"
        if("free  energy   TOTEN" in line):
            ll=line.split()
            F.append(float(ll[-2]))
        if("NELM" in line and ("=" in line)):
            ll=line.split()
            nelm=int(ll[2].strip(";"))
        if("magnetization (x)" in line):
            tag=0; Mab=0 ; mab=0; count=0
            for line in fin:
                ll=line.split()
                if(len(ll)>0):
                    if(ll[0]=="tot"):
                        Mag=ll[-1] #float(ll[-1])
                        break
                    elif(tag==0 and "----" in line):
                        tag=1
                    elif(tag==1  and "----" in line):
                        tag=0
                        Mab=mab
                    elif(tag==1):
                        mab += abs(float(ll[-1]))
                        if(count<n_mag):
                            atom_mag[count]=float(ll[-1])
                            count+=1

    aE=np.array(E)
    aF=np.array(F)
    ap=np.array(press)
    aV=np.array(V)
    Ep=aE /natom
    Fp=aF /natom
    
    return Ep[-1], press[-1]*0.1, atom_mag

elem=['Ce','Ni']
conf=["FM","NM",'AFM']
natom_mag=2

file_e=open("SF-energy.dat","w+")
file_p=open("SF-press.dat","w+")
file_m=open("SF-mag.dat","w+")

print("%12s"%("#E(eV/atom)"),end=" ",file=file_e)
print("%12s"%("#P(GPa) "),end=" ",file=file_p)
print("%12s"%("#M(uB/atom)"),end=" ",file=file_m)
for ic in conf:
    print("%10s"%(ic),end=" ",file=file_e)
    print("%6s"%(ic),end=" ",file=file_p)
    print("%10s"%(ic),end=" ",file=file_m)
print(" ",file=file_e)
print(" ",file=file_p)
print(" ",file=file_m)

for el_i in elem:
    print("%12s"%(el_i),end=" ",file=file_e)
    print("%12s"%(el_i),end=" ",file=file_p)
    print("%12s"%(el_i),end=" ",file=file_m)
    for ic in conf:
        E,P,mag=get_EM(ic+"/"+el_i+"/OUTCAR",natom_mag)
        mag=np.array(mag)
        print("%10.6f"%(E),end=" ",file=file_e)
        print("%6.1f"%(P),end=" ",file=file_p)
        ss="%.1f,%.1f"%(sum(mag)/natom_mag,sum(abs(mag))/natom_mag)
        print("%10s"%(ss),end=" ",file=file_m)
    print(" ",file=file_e)
    print(" ",file=file_p)
    print(" ",file=file_m)
file_e.close()
file_p.close()
file_m.close()
