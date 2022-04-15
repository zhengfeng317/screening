#!/usr/bin/env python3

import sys
import MD

fin=open(sys.argv[1],"r")
ll=fin.readline().split()
ele=[ll[k] for k in range(1,4)]
pts=[]
for line in fin:
    ll=line.split()
    n1=int(ll[1]); n2=int(ll[2]); n3=int(ll[3])
    Ef=float(ll[4])
    if(Ef>0.00001):
        print("Positive Ef! EXIT!!!")
        exit()
    pts.append([n1,n2,n3, Ef])

MD.draw_ternary_convex(pts,ele,ele[0]+"-"+ele[1])
