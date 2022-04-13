#!/usr/bin/env python3

import glob


def get_bad(files):
    lst=[] 
    for ifl in files:
        fin=open(ifl,"r")
        lines=fin.readlines()
        fin.close()
        if(len(lines)==0):
            print(ifl,"no data")
        elif("reached required" not in lines[-1]):
            sys=ifl.split("/")[1]
            print(ifl,'not converged')
            lst.append(sys)
        else:
            print(ifl,'done')
    return lst

files=glob.glob("FM/*/vasp.3")
lst_NM=get_bad(files)
