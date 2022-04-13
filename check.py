#!/usr/bin/env python3

import glob


def get_bad(files):
    lst=[] 
    for ifl in files:
        fin=open(ifl,"r")
        lines=fin.readlines()
        if(len(lines)==0):
            fin.close()
            ifl1=ifl.strip("2")+"1"
            fin1=open(ifl1,"r")
            lines=fin1.readlines()
            if(len(lines)==0):
                print(ifl,"no data")
            fin1.close()
        else:
            fin.close()
        if("reached required" not in lines[-1]):
            sys=ifl.split("/")[1]
            lst.append(sys)
    return lst

files=glob.glob("FM/*/vasp.3.2")
lst_FM=get_bad(files)

files=glob.glob("NM/*/vasp.3.2")
lst_NM=get_bad(files)

lst_comm=list(set(lst_FM).intersection(lst_NM))
lst_FM_only=[m for m in lst_FM if m not in lst_comm]
lst_NM_only=[m for m in lst_NM if m not in lst_comm]

print("common:")
for sys in sorted(lst_comm):
    print(sys)
print("")

print("FM only")
for sys in sorted(lst_FM_only):
    print(sys)
print("")

print("NM only")
for sys in sorted(lst_NM_only):
    print(sys)
