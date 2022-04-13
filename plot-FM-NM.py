#!/usr/bin/env python3

import numpy as np
import pylab as plt
plt.figure(figsize=(9,12))
plt.rcParams.update({'font.size': 14})

plt.subplot(2,1,1)
name=[]; dE=[]
with open('energy.dat','r') as f:
    f.readline()
    for line in f:
        ll=line.split()
        dE.append(float(ll[1]) - float(ll[2])) # E(FM)-E(NM)
        name.append(ll[0])
plt.xticks(range(len(name)), name, fontsize=12)
plt.plot(dE,'-o',label="E(FM)-E(NM)")
plt.ylabel('dE (eV/atom)')
plt.axhline(0,ls='--')
plt.legend()


plt.subplot(2,1,2)
name=[]; mag=[]
with open('mag.dat','r') as f:
    f.readline()
    for line in f:
        ll=line.split()
        mag_str=ll[1].split(",")
        mag_a=[float(k) for k in mag_str]
        mag.append(np.mean(mag_a))
        name.append(ll[0])
plt.xticks(range(len(name)), name, fontsize=12)
plt.plot(mag,'-o')
plt.ylabel('Mag. Mom. (uB/Metal)')
plt.axhline(1.0,ls='--')

plt.tight_layout()
plt.savefig('FM-NM.png')



