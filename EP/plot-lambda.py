#!/usr/bin/env python3

import numpy as np
import pylab as plt
plt.figure(figsize=(9,12))
plt.rcParams.update({'font.size': 16})

name=[]; lmd_max=[] ; lmd_sum=[]
with open('lambda.dat','r') as f:
    f.readline()
    for line in f:
        ll=line.split()
        lmd_max.append(float(ll[1]))
        lmd_sum.append(float(ll[2]))
        name.append(ll[0])

plt.subplot(2,1,1)
plt.xticks(range(len(name)), name, fontsize=12)
plt.plot(lmd_max,'-o',label="$\lambda$_max")
plt.ylabel('$\lambda_{\Gamma}$')
plt.ylim(-0.1,2)
plt.axhline(0,ls='--')
plt.legend()


plt.subplot(2,1,2)
plt.xticks(range(len(name)), name, fontsize=12)
plt.plot(lmd_sum,'-o',label="$\lambda$_sum")
plt.ylabel('$\lambda_{\Gamma}$')
plt.ylim(-0.1,2)
plt.axhline(0,ls='--')
plt.legend()

plt.tight_layout()
plt.savefig('lambda.png')



