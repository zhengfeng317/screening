#!/usr/bin/env python3

def cal_Ef(pts):
    Ea=0; Eb=0; Ec=0
    for ipt in pts:
        if(ipt[1]==0 and ipt[2]==0):
            Ea=ipt[3]
        elif(ipt[0]==0 and ipt[2]==0):
            Eb=ipt[3]
        elif(ipt[0]==0 and ipt[1]==0):
            Ec=ipt[3]
    if(Ea==0 or Eb==0 or Ec==0):
        print("Error: cannot find formation energy reference!!! (Ea,Eb,Ec=", Ea, Eb, Ec,")")
        print(pts)
    pts_Ef=[]
    for ipt in pts:
        Ef=ipt[3]-(Ea*ipt[0]+Eb*ipt[1]+Ec*ipt[2])/(ipt[0]+ipt[1]+ipt[2])
        pts_Ef.append([ipt[0], ipt[1], ipt[2], Ef])
    return pts_Ef


def area(a,b,c):
    import numpy as np
    from numpy.linalg import norm
    if(type(a)==list):
        a=np.array(a,dtype=np.float32)
    if(type(b)==list):
        b=np.array(b,dtype=np.float32)
    if(type(c)==list):
        c=np.array(c,dtype=np.float32)
    return 0.5*norm(np.cross(b-a,c-a))

def get_plane(p1,p2,p3):
    import numpy as np
    if(type(p1)==list):
        p1=np.array(p1,dtype=np.float32)
    if(type(p2)==list):
        p2=np.array(p2,dtype=np.float32)
    if(type(p3)==list):
        p3=np.array(p3,dtype=np.float32)
    v1=p3-p1; v2=p2-p1
    cp=np.cross(v1,v2)
    a,b,c=cp
    d=np.dot(cp,p3)
    return a,b,c,d

def ternary_convex_Ed(pts):
    import numpy as np
    from scipy.spatial import ConvexHull
    
    # convert data to trianle set
    pts=np.array(pts)
    tpts = []
    for ipt in pts:
        comp=ipt[:3]
        comp=comp/sum(comp)
        x=comp[0] + comp[1]/2.
        y=comp[1]*np.sqrt(3)/2
        tpts.append([x,y,ipt[3]])

    # get convex hull
    # hll.vertices is the indice of stable compounds in pts,
    # hull.simplices is the connections among stable compounds
    hull=ConvexHull(tpts)
    stables=[]
    for iv in hull.vertices:
        name=ele[0]+str(int(pts[iv][0]))+ele[1]+str(int(pts[iv][1]))+ele[2]+str(int(pts[iv][2]))
        stables.append([tpts[iv][0],tpts[iv][1],name])  # still not sure how to plot names on the figure 06/24

    # plot data
    pdata=[]
    for pt in pts:
        mm=pt[:3]
        pdata.append(1.0*mm/sum(mm))
    
    #2 get meta-stable phases
    mstables=[]
    for i in range(len(pdata)):
        if(i not in hull.vertices):
            mstables.append(pdata[i])
    
    #3 plot meta-stable phases
    mstables=np.array(mstables)
    
    #4 find the distance to the convex hull
    #  4.1 get nearest 3 points
    Ed=999
    for k in range(1):   # only care the first structure
        if(k in hull.vertices):
            Ed=0
            continue # jump the stable ones
        x=tpts[k][:2] # metastable, as [x,y,Ef]
        heights=[]
        for isimp in hull.simplices: # loop the simplices
            A=tpts[isimp[0]][:2]; B=tpts[isimp[1]][:2]; C=tpts[isimp[2]][:2];
            # find if x in the A-B-C triangle
            area_ABC=area(A,B,C)
            sum_a = area(A,B,x)+area(A,C,x)+area(B,C,x)
            if( sum_a - area_ABC <=0.001):
                # in the ABC, get the ABC plane
                a,b,c,d=get_plane(tpts[isimp[0]],tpts[isimp[1]],tpts[isimp[2]])
                if(a==0 and b==0 and d==0):
                    continue
                if(c==0):
                    continue
                # get the cross point with ABC plane
                z = (d-a*x[0]-b*x[1])/c
                # height to convex hull
                h = tpts[k][2]-z
                heights.append(h)
        Ed=max(heights)
    return Ed

# main
fin=open("mystr.dat",'r')
print("%10s%16s"%("System", " Ed(eV/atom)"))
for line in fin:
    # info of one structure
    ll=line.split()
    if(len(ll)!=5):
        break
    syst=ll[0].split("/")[-1] # system name
    print("%10s"%(syst),end=" ") # system name here
    E=float(ll[1]) # energy
    chem={}
    for k in range(2,5):
        ss=ll[k].split("_")
        chem[ss[0]]=int(ss[1])
    #print(sys, E, chem)
    ele=list(chem.keys())

    pts=[]
    pts.append([chem[ele[0]], chem[ele[1]], chem[ele[2]], E])
    names=["new"]
    # find ref
    fr=open("../"+syst+"_stable.dat","r")
    for lr in fr:
        sid=lr.split()[1].split("-")[1]
        names.append(lr.split()[1])
        fr_my=open("ref.dat","r")
        ref_info=""
        for lrm in fr_my:
            sys=lrm.split()[0].split("/")[1].split(".")[0]
            if(sys==sid):
                ref_info=lrm
        fr_my.close()

        ref_info=ref_info.split()
        E=float(ref_info[1])
        chem=chem.fromkeys(chem,0)
        for k in range(2,len(ref_info)):
            ss=ref_info[k].split("_")
            chem[ss[0]]=int(ss[1])
        pts.append([chem[ele[0]], chem[ele[1]], chem[ele[2]], E])
    fr.close()
    
    # get Ef
    pts_Ef=cal_Ef(pts)
    if(pts_Ef[0][-1]>=0): # Ef must not be positive to get Ed!!!
        Ef_post=pts_Ef[0][-1]
        pts_Ef[0][-1]=-0.00001 # use almost ZERO to replace positive Ef
        Ed=ternary_convex_Ed(pts_Ef)
        Ed+=Ef_post
        pts_Ef[0][-1]=Ef_post
    else:
        Ed=ternary_convex_Ed(pts_Ef) 
    print("%16.6f"%(Ed)) # system name here

    # output final data
    fout=open("Ef_"+syst+".dat","w+")
    print("%16s"%("system"),end=" ",file=fout)
    for k in range(len(ele)):
        print("%6s"%(ele[k]),end=" ",file=fout)
    print("%12s"%("Ef(eV/atom)"),file=fout)
    for k in range(len(pts_Ef)):
        print("%16s"%(names[k]),end=" ",file=fout)
        for j in range(len(ele)):
            print("%6d"%(pts_Ef[k][j]),end=" ",file=fout)
        print("%12.6f"%(pts_Ef[k][-1]),file=fout)
    fout.close()
