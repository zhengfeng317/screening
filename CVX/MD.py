##--Python functions to analyze simulation results, with LAMMPS&VASP files---##
##--Yang Sun(yangsun017@gmail.com), OCT-16-2016 -----------------------------##
##--version 1.0.0----------------------------- ------------------------------##
##--Contents-----------------------------------------------------------------##
##---------xdat2movie(xfile,njump)-------------------------------------------##
##---------math functions----------------------------------------------------##
##----------------------AcrossB(a,b)-----------------------------------------##
##----------------------AdotB(a,b)-------------------------------------------##
##----------------------volume(a,b,c)----------------------------------------##
##----------------------R2K(x)-----------------------------------------------##
##----------------------rotation_matrix(x,y,z,ranx,rany,ranz)----------------##
##---------gr(xfile,njump,rrange,mesh)---------------------------------------##
##---------lamp2poscar(filein,step,atom_type)--------------------------------##
##---------dump2xdat(filein,atom_type)---------------------------------------##
##---------dump2datain(filein,step,mass)-------------------------------------##
##---------lmp_random_start(vector,mass,composition)-------------------------##
##---------BOO_Q(l,x,y,z)----------------------------------------------------##
##---------W3J(J1,J,J2,M1,M,M2)----------------------------------------------##
##---------BOO_WL(l,x,y,z)---------------------------------------------------##
##---------BOO_WL_bar(l,x,y,z)-----------------------------------------------##
##---------freq(x,xbin,xstart,xend)------------------------------------------##
##---------cry_elong(file_in,elongate)---------------------------------------##
##---------lammps2posave(filein,atom_type)-----------------------------------##
##---------cry_cluster(file_in, name_id, ntemp,center, neighbor)-------------##
##---------xdat2dump---------------------------------------------------------##
##---------lammps2xml--------------------------------------------------------##
##---------dumpave-----------------------------------------------------------##
##---------draw_ternary------------------------------------------------------##
##---------draw_ternary_convex(pts,ele,string)-------------------------------##
##---------draw_binary_convex(pts,ele,string)--------------------------------##
##---------vasp2spin---------------------------------------------------------##
##---------count_n_plot(list)------------------------------------------------##
##---------cry_cluster_cut---------------------------------------------------##
##---------RAS---------------------------------------------------------------##

#-------------------------------------XDATCAR to xyz fomat-----------------------------
def xdat2movie(xfile,njump,ntotal_step):    
#transfer XDATCAR(vasp output) to movie.dat
#xfile = input file(must XDATCAR format); njump = the number of withdrawed steps
   fin=open(xfile,'r')
   xlines=fin.readlines(700)
#atom number
   atom_type=xlines[5].split()
   ntype=len(atom_type)
   na=[]
   for i in xlines[6].split():
      na.append(int(i))
   natom=sum(na)
#steps
   nstep=ntotal_step-njump
#vectors
   v=[[0 for i in range(3)] for j in range(3)]
   line_v=[[0 for i in range(3)] for j in range(3)]
   for i in range(3):
      line_v[i]=xlines[i+2].split()
   for i in range(3):
      for j in range(3):
         v[i][j]=float(line_v[i][j])
#atom type
   fin.seek(0)
#output xyz data
   fout=open('movie.dat','w+')
   x_dir=[0 for i in range(3)]
   x_car=[0 for i in range(3)]
   for i in range(7+njump*(natom+1)):
      fin.readline()
   for i in range(nstep):
      istep=i+1 #int(line_step[len(line_step)-1])
      print('%10d'%(natom),file=fout)
      print('%10d'%(istep),file=fout)
      fin.readline()
      for itype in range(ntype):
         for j in range(na[itype]):
            line_xyz=fin.readline().split()
            for k in range(3):
               x_dir[k]=float(line_xyz[k])
               if(x_dir[k] > 1): x_dir[k] = x_dir[k] -1.0  
               if(x_dir[k] < 0): x_dir[k] = x_dir[k] +1.0  
   #         if(x_dir[k] > 0.5): x_dir[k] = x_dir[k] -1.0  # traditional treatment from Xiaowei
            for k in range(3):
               x_car[k]=x_dir[0]*v[k][0]+x_dir[1]*v[k][1]+x_dir[2]*v[k][2]
            print('%2s%16.6f%16.6f%16.6f'%(atom_type[itype],x_car[0],x_car[1],x_car[2]),file=fout)
   fin.close()
   fout.close()
   return nstep

#------------------------------math_functions--------------------------------------------------
# a,b,c are three-element list
def AcrossB(a,b): 
   c=[0.0 for i in range(3)]
   c[0]=a[1]*b[2]-a[2]*b[1]
   c[1]=a[2]*b[0]-a[0]*b[2]
   c[2]=a[0]*b[1]-a[1]*b[0]
   return c

def AdotB(a,b):
   c = a[0]*b[0]+ a[1]*b[1] + a[2]*b[2]
   return c

def volume(a,b,c):
   x=AcrossB(b,c)
   v=AdotB(a,x)
   return v

def R2K(x): 
# x = x[3][3] 
   y=[[0 for i in range(3)] for j in range(3)]
   v=volume(x[0],x[1],x[2])
   b1=AcrossB(x[1],x[2])
   b2=AcrossB(x[2],x[0])
   b3=AcrossB(x[0],x[1])
   for i in range(3):
      b1[i]=b1[i]/v
      b2[i]=b2[i]/v
      b3[i]=b3[i]/v
   y[0]=b1
   y[1]=b2
   y[2]=b3
   return y

def rotation_matrix(x,y,z,ranx,rany,ranz):
   import numpy as np
   thx=ranx*np.pi/180
   thy=ranx*np.pi/180
   thz=ranx*np.pi/180
   length=len(x)
   xn=[]
   for i in range(length):
      xr1=x[i]*np.cos(thz)+y[i]*np.sin(thz)
      yr1=-x[i]*np.sin(thz)+y[i]*np.cos(thz)
      zr1=z[i]
      xr2=xr1*np.cos(thy)-zr1*np.sin(thy)
      yr2=yr1
      zr2=xr1*np.sin(thy)+zr1*np.cos(thy)
      xr1=xr2
      yr1=yr2*np.cos(thx)+zr2*np.sin(thx)
      zr1=-yr2*np.sin(thx)+zr2*np.cos(thx)
      xn.append([xr1,yr1,zr1])
   return(xn)
      
#---------------------------------math function done-----------------------------------------


#---------------------------Pair Correlation Functions---------------------------------------
def gr(xfile,njump,rrange,mesh):
   from math import sqrt
# Calculate g(r) for XDATCAR file
# xfile = input file(MUST DATCAR format);  njump= withdrawed steps before calculating(if 0, calculte all the data);
# rrange = the r range; mesh = the mesh for the r range
   fin=open(xfile,'r')
   fout=open('gr-vasp.dat','w')
   xlines=fin.readlines()
#atom number
   line_num=xlines[6].split()
   n1=int(line_num[0])
   n2=int(line_num[1])
   natom=n1+n2
#steps
   nstep=int((len(xlines)-7)/(natom+1)-njump)
#vectors
   v=[[0 for i in range(3)] for j in range(3)]
   line_v=[[0 for i in range(3)] for j in range(3)]
   for i in range(3):
      line_v[i]=xlines[i+2].split()
   for i in range(3):
      for j in range(3):
         v[i][j]=float(line_v[i][j])
#atom type
   atom_type=xlines[5].split()
   fin.close()
#atom coordinates: x[nsteps][natoms]
   fin=open(xfile,'r')
   for i in range(7+njump*(natom+1)):
      fin.readline()
   x=[[0.0 for i in range(natom)] for j in range(nstep)]
   y=[[0.0 for i in range(natom)] for j in range(nstep)]
   z=[[0.0 for i in range(natom)] for j in range(nstep)]
   for i in range(nstep):
      fin.readline()
      for j in range(natom):
         line_xyz=fin.readline().split()
         x[i][j]=float(line_xyz[0])
         y[i][j]=float(line_xyz[1])
         z[i][j]=float(line_xyz[2])
   fin.close()

#gr calculation
   twopi=6.283185307
   rangr2=rrange*rrange
   scale=mesh/rrange
   delta=rrange/mesh
   vol=volume(v[0],v[1],v[2])
   grall=[0 for i in range(mesh+1)]
   gr11=[0 for i in range(mesh+1)]
   gr12=[0 for i in range(mesh+1)]
   gr22=[0 for i in range(mesh+1)]
   rr=[0 for i in range(3)]
   recip=R2K(v)
   i0max=int(rrange * sqrt(AdotB(recip[0],recip[0])) - 0.001 )
   i1max=int(rrange * sqrt(AdotB(recip[1],recip[1])) - 0.001 )
   i2max=int(rrange * sqrt(AdotB(recip[2],recip[2])) - 0.001 )
   for i in range(nstep):
      for i0 in range(-1*i0max-1,i0max+1):
         for i1 in range(-1*i1max-1,i1max+1):
            for i2 in range(-1*i2max-1,i2max+1):
               for ni in range(natom-1):
                  for nii in range(ni+1,natom):
                     d0=i0+(x[i][ni]-x[i][nii]+1)%1.0
                     d1=i1+(y[i][ni]-y[i][nii]+1)%1.0
                     d2=i2+(z[i][ni]-z[i][nii]+1)%1.0
                     for ii in range(3):
                        rr[ii]=d0*v[0][ii]+d1*v[1][ii]+d2*v[2][ii]
                     r2=AdotB(rr,rr)
                     if(r2< rangr2) :
                        index=int(sqrt(r2))
                        grall[index] += 1
                        if ( (ni < n1) and (nii < n1)): gr11[index] += 1
                        if ( (ni < n1) and (nii >= n1)): gr12[index] += 1
                        if ( (ni >= n1) and (nii >= n1)): gr22[index] += 1
      print(i)
   pcfak=3.0/twopi/natom**2*vol*scale**3
   pcfak11=3.0/twopi/n1**2*vol*scale**3
   pcfak12=3.0/twopi/(n1*n2)/2*vol*scale**3
   pcfak22=3.0/twopi/n2**2*vol*scale**3
   for i in range(mesh):
      rdist=i*delta
      grall[i]=pcfak*grall[i]/nstep/(3*i*(i+1)+1)
      gr11[i]=pcfak11*gr11[i]/nstep/(3*i*(i+1)+1)
      gr12[i]=pcfak12*gr12[i]/nstep/(3*i*(i+1)+1)
      gr22[i]=pcfak22*gr22[i]/nstep/(3*i*(i+1)+1)
      print('%10.4f%16.6f%16.6f%16.6f%16.6f'%(rdist,grall[i],gr11[i],gr12[i],gr22[i]),file=fout)
   return 0
#-----------------------End: pair correlation functions--------------------------------------

#---------------------------------lammps2poscar----------------------------------------------
# Transfer LAMMPS dump file to POSCAR
def lamp2poscar(filein,step,atom_type):
   fin=open(filein,'r')
   for i in range(step-1):
      for k in range(3):
         fin.readline()
      natom=int(fin.readline())
      for k in range(5+natom):
         fin.readline()
   fin.readline()
   atep=int(fin.readline().strip('\n'))
   fin.readline()
   natom=int(fin.readline())
   fin.readline()


#   if(step== 1):
#      fin.readline()
#      atep=int(fin.readline().strip('\n'))
#      fin.readline()
#      natom=int(fin.readline())
#      fin.readline()
#   else:
#      for i in range(3): fin.readline()      
#      natom=int(fin.readline())
#      for i in range((step-1)*(natom+9)-3): fin.readline()
#      atep=int(fin.readline().strip('\n'))
#      for i in range(3):
#          fin.readline()
   print(atep)
# vector
   v=[[0 for i in range(3)] for j in range(3)]
   vlow=[0 for i in range(3)]
   vhigh=[ 0 for i in range(3)]
   l1=fin.readline().split()
   l2=fin.readline().split()
   l3=fin.readline().split()
   vlow[0]=float(l1[0]); vlow[1]=float(l2[0]); vlow[2]=float(l3[0])
   vhigh[0]=float(l1[1]); vhigh[1]=float(l2[1]); vhigh[2]=float(l3[1])
   if(len(l1)==2) :
      v[0][0]=float(l1[1])-float(l1[0])
      v[1][1]=float(l2[1])-float(l2[0])
      v[2][2]=float(l3[1])-float(l3[0])
   else :
      # see details on https://lammps.sandia.gov/doc/Howto_triclinic.html
      xlo_b=float(l1[0]); xhi_b=float(l1[1]); xy=float(l1[2])
      ylo_b=float(l2[0]); yhi_b=float(l2[1]); xz=float(l2[2])
      zlo_b=float(l3[0]); zhi_b=float(l3[1]); yz=float(l3[2])
      xlo=xlo_b-min(0.0,xy,xz,xy+xz)
      xhi=xhi_b-max(0.0,xy,xz,xy+xz)
      ylo=ylo_b-min(0.0,yz)
      yhi=yhi_b-max(0.0,yz)
      zlo=zlo_b; zhi=zhi_b
      v[0][0]=xhi-xlo
      v[1][1]=yhi-ylo
      v[2][2]=zhi-zlo
      v[1][0]=xy
      v[2][0]=xz
      v[2][1]=yz
# xyz type ......
   arg_id=[];arg_type=[];arg_xs=[];arg_ys=[];arg_zs=[];arg_x=[];arg_y=[];arg_z=[]
   arg_fx=[];arg_fy=[];arg_fz=[];arg_cpa=[]; arg_vx=[]; arg_vy=[]; arg_vz=[];arg_cpa2=[];
   arg_cka=[]; arg_diff=[]; arg_nb=[]; arg_V=[]; arg_sro=[]; arg_rmsd=[]; arg_QConn=[];
   arg_xu=[];arg_yu=[];arg_zu=[]
 #  dic={arg_id:'id',arg_type:'type',arg_xs:'xs',arg_ys:'ys',arg_zs:'zs',arg_x:'x',arg_y:'y',arg_z:'z'}
   dic={'id':arg_id,'type':arg_type,'xs':arg_xs,'ys':arg_ys,'zs':arg_zs,'x':arg_x,'y':arg_y,'z':arg_z,
           'fx':arg_fx,'fy':arg_fy,'fz':arg_fz,'c_peratom':arg_cpa,'c_pperatom':arg_cpa, 'vx':arg_vx, 'vy':arg_vy, 'vz':arg_vz,
           'c_pea':arg_cpa2,'c_kea':arg_cka, 'diffusion':arg_diff, 'nb':arg_nb, 'V':arg_V, 'sro':arg_sro,
           'RMSD':arg_rmsd, 'c_QConn':arg_QConn, 'StructureType':arg_type, 'xu':arg_xu,'yu':arg_yu,'zu':arg_zu}
   l_arg=fin.readline().split()
   arg_num=len(l_arg)-2
   for i in range(natom):
      line=fin.readline().split()
      for j in range(arg_num):
         dic[l_arg[2+j]].append(line[j])
# output
   if(len(arg_type)==0): 
      print('Error: Include "type" in your dump commond!')
      return 0
 #  ntype=int(arg_type[natom-1])
   ntype=len(atom_type)
   n=[0 for i in range(ntype)]
   for i in range(ntype):  n[i]=arg_type.count(str(i+1))
   #for i in range(ntype):  n[i]=arg_type.count(str(i)) # w/ type==0
   fout=open('snapshot.vasp','w')
   print('Lammps2Vasp',file=fout)
   print('%5.3f'%(1.0),file=fout)
   # output vector
   for i in range(3):  print('%16.6f%16.6f%16.6f'%(v[i][0],v[i][1],v[i][2]),file=fout)
   # output atom species
   for i in range(ntype): print('%5s'%(atom_type[i]),end=' ',file=fout)
   print(' ',file=fout)
   # output atom numbers for each species
   for i in range(ntype): print('%5d'%(n[i]),end=' ',file=fout)
   print(' ',file=fout)
   # output atom coordinates
   if(len(arg_xs) != 0 ):
      print('Direct',file=fout)
      for k in range(ntype):
         for i in range(natom): 
            if(arg_type[i]==str(k+1)): print(arg_xs[i],arg_ys[i],arg_zs[i],arg_id[i],file=fout)
   elif(len(arg_x) !=0) :
      print('Direct',file=fout)
      for k in range(ntype):
         for i in range(natom): 
            if(arg_type[i]==str(k+1)): 
            #if(arg_type[i]==str(k)):  # w/ type==0
                dx= (float(arg_x[i])-vlow[0])/v[0][0]
                dy= (float(arg_y[i])-vlow[1])/v[1][1]
                dz= (float(arg_z[i])-vlow[2])/v[2][2]
                print(dx,dy,dz,arg_id[i],file=fout)
   elif(len(arg_xu) !=0) :
      print('Direct',file=fout)
      for k in range(ntype):
         for i in range(natom): 
            if(arg_type[i]==str(k+1)): 
                dx= (float(arg_xu[i])-vlow[0])/v[0][0]
                dy= (float(arg_yu[i])-vlow[1])/v[1][1]
                dz= (float(arg_zu[i])-vlow[2])/v[2][2]
                dx=dx-int(dx)
                dy=dy-int(dy)
                dz=dz-int(dz)
                print(dx,dy,dz,arg_id[i],file=fout)
   if(len(arg_vx) != 0):
      print("",file=fout)
      for k in range(ntype):
         for i in range(natom): 
            if(arg_type[i]==str(k+1)):
               print(float(arg_vx[i])*0.001,float(arg_vy[i])*0.001,
                     float(arg_vz[i])*0.001,arg_id[i],file=fout) # A/fs in vasp
#      print('Cartesian',file=fout)
#      for k in range(ntype):
#         for i in range(natom): 
#            if(arg_type[i]==str(k+1)): print(arg_x[i],arg_y[i],arg_z[i],i+1,file=fout)
      return 2
#--------------------End: Lmp_dump 2 POSCAR------------------------------------------------

#--------------------dump file to xdatcar type--------------------------------------------- 
def dump2xdat(filein,atom_type,nevery):
    fin=open(filein,'r')
    for i in range(3): fin.readline()
    natom=int(fin.readline())
    line=fin.readline()
# vector
    steps=0
    v=[[0 for i in range(3)] for j in range(3)]
    while(line !=''):
        l1=fin.readline().split()
        l2=fin.readline().split()
        l3=fin.readline().split()
        if(len(l1)==2) :
            v[0][0] += float(l1[1])-float(l1[0])
            v[1][1] += float(l2[1])-float(l2[0])
            v[2][2] += float(l3[1])-float(l3[0])
            steps += 1
        else :
            v[0][0] += float(l1[1])-float(l1[0])
            v[1][1] += float(l2[1])-float(l2[0])
            v[2][2] += float(l3[1])-float(l3[0])
            v[1][0] += float(l1[2])
            v[2][0] += float(l2[2])
            v[2][1] += float(l3[2])
            steps += 1
        for i in range(natom+6): line=fin.readline()
    if(len(l1)==2):
        v[0][0] /= steps
        v[1][1] /= steps
        v[2][2] /= steps
    else : 
        v[0][0] /= steps
        v[1][1] /= steps
        v[2][2] /= steps
        v[1][0] /= steps
        v[2][0] /= steps
        v[2][1] /= steps
    fin.close()
    print("total steps:", steps)
    fin=open(filein,'r')
    fout=open('XDATCAR','w')
    ntype=len(atom_type)
    print('Dumpfile -> XDATCAR',file=fout)
    print('%5.3f'%(1.0),file=fout)
   # output vector
    for i in range(3):  print('%16.6f%16.6f%16.6f'%(v[i][0],v[i][1],v[i][2]),file=fout)
   # output atom species
    for i in range(ntype): print('%5s'%(atom_type[i]),end=' ',file=fout)
    print(' ',file=fout)

    for i in range(8): fin.readline()
    # xyz type ......
    arg_id=[];arg_type=[];arg_xs=[];arg_ys=[];arg_zs=[];arg_x=[];arg_y=[];arg_z=[]
    arg_fx=[];arg_fy=[];arg_fz=[];arg_cpa=[];arg_vx=[];arg_vy=[];arg_vz=[];arg_conn=[]
 #  dic={arg_id:'id',arg_type:'type',arg_xs:'xs',arg_ys:'ys',arg_zs:'zs',arg_x:'x',arg_y:'y',arg_z:'z'}
    dic={'id':arg_id,'type':arg_type,'xs':arg_xs,'ys':arg_ys,'zs':arg_zs,'x':arg_x,'y':arg_y,'z':arg_z,
            'fx':arg_fx,'fy':arg_fy,'fz':arg_fz,'c_peratom':arg_cpa, 'c_pea':arg_cpa, 'vx':arg_vx, 'vy':arg_vy, 'vz':arg_vz,
            'c_QConn':arg_conn}
    l_arg=fin.readline().split()
    arg_num=len(l_arg)-2
    for i in range(natom):
       line=fin.readline().split()
       for j in range(arg_num):
          dic[l_arg[2+j]].append(line[j])
    n=[0 for i in range(ntype)]
    for i in range(ntype):  n[i]=arg_type.count(str(i+1))
   # output atom numbers for each species
    for i in range(ntype): print('%5d'%(n[i]),end=' ',file=fout)
    print(' ',file=fout)
    print('Direct configuration= 0',file=fout)
    if(len(arg_xs) != 0 ):
        for k in range(ntype):
            for i in range(natom): 
                if(arg_type[i]==str(k+1)): print(arg_xs[i],arg_ys[i],arg_zs[i],i+1,file=fout)
        for ii in range(steps-1):
            if(ii%nevery != 0):
                for jj in range(9+natom):
                    fin.readline()
            else:
                print(str(ii)+"/"+str(steps))
                print('Direct configuration= ',ii+1,file=fout)
                for jj in range(9):fin.readline()
                for i in range(natom):
                    line=fin.readline().split()
                    for j in range(arg_num):
                        dic[l_arg[2+j]][i]=line[j]
                for k in range(ntype):
                    for i in range(natom):
                        if(arg_type[i]==str(k+1)): print(arg_xs[i],arg_ys[i],arg_zs[i],i+1,file=fout)
        return 1
    else:
        return ('use xs in lammps input')
#--------------------End: Lmp_dump 2 XDATCAR------------------------------------------------




#--------------------dump to data.in-------------------------------------------------------
def dump2datain(filein,step,mass):
   fin=open(filein,'r')
   fin.readline()
   astep=fin.readline().strip('\n')
   fin.readline()
   natom=int(fin.readline())
   fin.readline()
   if( step != 1):
       for i in range((step-1)*(natom+9)-4): fin.readline()
       astep=fin.readline().strip('\n')
       for i in range(3): fin.readline()
# vector
   v=[[0 for i in range(3)] for j in range(3)]
   l1=fin.readline().split()
   l2=fin.readline().split()
   l3=fin.readline().split()
   if(len(l1)==2) :
      v[0][0]=float(l1[1])-float(l1[0])
      v[1][1]=float(l2[1])-float(l2[0])
      v[2][2]=float(l3[1])-float(l3[0])
   else :
      v[0][0]=float(l1[1])-float(l1[0])
      v[1][1]=float(l2[1])-float(l2[0])
      v[2][2]=float(l3[1])-float(l3[0])
      v[1][0]=float(l1[2])
      v[2][0]=float(l2[2])
      v[2][1]=float(l3[2])
   bl=[float(l1[0]),float(l2[0]),float(l3[0])]
   bh=[float(l1[1]),float(l2[1]),float(l3[1])]
   #tilt=[float(l1[2]),float(l2[2]),float(l3[2])]
# xyz type ......
   arg_id=[];arg_type=[];arg_xs=[];arg_ys=[];arg_zs=[];arg_x=[];arg_y=[];arg_z=[]
   arg_fx=[];arg_fy=[];arg_fz=[];arg_vx=[];arg_vy=[];arg_vz=[]; arg_cpa=[]; arg_cka=[];
   arg_cpa2=[]; arg_nconn=[] ; arg_vq6=[]; arg_c_cluster=[]
 #  dic={arg_id:'id',arg_type:'type',arg_xs:'xs',arg_ys:'ys',arg_zs:'zs',arg_x:'x',arg_y:'y',arg_z:'z'}
   dic={'id':arg_id,'type':arg_type,'xs':arg_xs,'ys':arg_ys,'zs':arg_zs,'x':arg_x,'y':arg_y,'z':arg_z,
           'fx':arg_fx,'fy':arg_fy,'fz':arg_fz, 'vx':arg_vx, 'vy':arg_vy, 'vz':arg_vz, 'c_peratom':arg_cpa,
           'c_kea':arg_cka, 'c_pea':arg_cpa2, 'c_QConn':arg_nconn,  'v_q6Coor':arg_vq6, 'c_cluster_index':arg_c_cluster }
   l_arg=fin.readline().split()
   arg_num=len(l_arg)-2
   for i in range(natom):
      line=fin.readline().split()
      for j in range(arg_num):
         dic[l_arg[2+j]].append(line[j])
# output
   if(len(arg_type)==0): 
      print('Error: Include "type" in your dump commond!')
      return 0
 #  ntype=int(arg_type[natom-1])
   ntype=len(mass)
   n=[0 for i in range(ntype)]
   for i in range(ntype):  n[i]=arg_type.count(str(i+1))
   fout=open('dump2data.in','w')
   print('Dump2datain,step:'+astep,file=fout)
   print(' ',file=fout)
   print('%10d%8s'%(natom,'atoms'),file=fout)
   print('%5d%14s'%(ntype,'atom types'),file=fout)
   print(' ',file=fout)
   print(l1[0],l1[1],'xlo xhi',file=fout)
   print(l2[0],l2[1],'ylo yhi',file=fout)
   print(l3[0],l3[1],'zlo zhi',file=fout)
   print(' ',file=fout)
   if( len(l1) == 3):
      print(l1[2],l2[2],l3[2],'xy xz yz',file=fout)
      print(' ',file=fout)
   print(' Masses',file=fout)
   print(' ',file=fout)
   for i in range(ntype):
       print('%5d%10.3f'%(i+1,mass[i]),file=fout)
   print(' ',file=fout)
   # output atom coordinates
   print(' Atoms',file=fout)
   print(' ',file=fout)
   if(len(arg_x) !=0 ) :
       for i in range(natom): 
           print("%s  %s  %12.6f  %12.6f  %12.6f"%(arg_id[i],arg_type[i],float(arg_x[i]),float(arg_y[i]),float(arg_z[i])),file=fout)
   elif(len(arg_xs) != 0 ):
       for i in range(natom): 
           xx=bl[0]+float(arg_xs[i])*v[0][0]+float(arg_ys[i])*v[1][0]+float(arg_zs[i])*v[2][0]
           yy=bl[1]+float(arg_xs[i])*v[0][1]+float(arg_ys[i])*v[1][1]+float(arg_zs[i])*v[2][1]
           zz=bl[2]+float(arg_xs[i])*v[0][2]+float(arg_ys[i])*v[1][2]+float(arg_zs[i])*v[2][2]
#           print(arg_id[i],arg_type[i],xx,yy,zz,file=fout)
           print("%s  %s  %12.6f  %12.6f  %12.6f"%(arg_id[i],arg_type[i],xx,yy,zz),file=fout)
   # atom velocity
   if(len(arg_vx) !=0):
       print('Velocities',file=fout)
       print(' ',file=fout)
       for i in range(natom):
           print(arg_id[i],arg_vx[i],arg_vy[i],arg_vz[i],file=fout)
       print(' ',file=fout)
   return 1
#--------------------End: Lmp_dump 2 LMP_data.in------------------------------------------------
       
#--------------------Lammps_random_configuration_generator--------------------------------------
def lmp_random_start(vector,mass,composition):
# vector=vector[[],[],[]], mass=mass[...],composition=[...]
    import random
#check
    if(len(mass)!= len(composition)):
        print('Input Error: Composition Wrong')
        return 0
#start
    ntype=len(mass)
    natom=0
    for i in range(ntype): natom += composition[i]
    fin=open('data_random.in','w')
    print('Randomly generated, composition:',end=' ',file=fin)
    for i in range(ntype):
        print(int(composition[i]),end=' ',file=fin)
    print(' ',file=fin)
    print(' ',file=fin)
    print(int(natom),'    atoms',file=fin)
    print(int(ntype),'    atom types',file=fin)
    print(' ',file=fin)
    print('%12.6f%12.6f%8s'%(0.0,vector[0][0],'xlo xhi'),file=fin)
    print('%12.6f%12.6f%8s'%(0.0,vector[1][1],'ylo yhi'),file=fin)
    print('%12.6f%12.6f%8s'%(0.0,vector[2][2],'zlo zhi'),file=fin)
    print(' ',file=fin)
    print('%12.6f%12.6f%12.6f%12s'%(vector[1][0],vector[2][0],vector[2][1],'xy xz yz'),file=fin)
    print(' ',file=fin)
    print(' Masses',file=fin)
    print(' ',file=fin)
    for i in range(ntype):
        print('%5d%10.3f'%(i+1,mass[i]),file=fin)
    print(' ',file=fin)
    print(' Atoms',file=fin)
    print(' ',file=fin)
    atype=[]
    for i in range(ntype):
        for k in range(composition[i]):
            atype.append(str(i+1))
    x=[];y=[];z=[]
    i=0
    while (i<natom):
        xs=random.random();ys=random.random();zs=random.random()
        index=0
        for j in range(i):
            dx=xs-x[j];dy=ys-y[j];dz=zs-z[j]
            xx=dx*vector[0][0]+dy*vector[1][0]+dz*vector[2][0]
            yy=dx*vector[0][1]+dy*vector[1][1]+dz*vector[2][1]
            zz=dx*vector[0][2]+dy*vector[1][2]+dz*vector[2][2]
            dr2=xx**2+yy**2+zz**2
            if(dr2 < 0) : #1.8**2
                index=1
                break
        if(index == 0):
            i += 1
            x.append(xs);y.append(ys),z.append(zs)
    for i in range(natom):
        xx=x[i]*vector[0][0]+y[i]*vector[1][0]+z[i]*vector[2][0]
        yy=x[i]*vector[0][1]+y[i]*vector[1][1]+z[i]*vector[2][1]
        zz=x[i]*vector[0][2]+y[i]*vector[1][2]+z[i]*vector[2][2]
        print(i+1,atype[i],xx,yy,zz,file=fin)
#        print('%5d%5s%16.6f%16.6f%16.6f'%(i+1,atype[i],xx,yy,zz),file=fin)
    return 1
#--------------------End: LRC_generator-----------------------------------------------------   

#--------------------Bond orientation order: qlm(m,l,theta,phi)
def BOO_QLM(m,l,phi,theta):
   import numpy as np
   import scipy.special
   nb=len(phi)
   qlm=complex(0)
   for i in range(nb):
      aa=scipy.special.sph_harm(m,l,phi[i],theta[i])
      qlm += scipy.special.sph_harm(m,l,phi[i],theta[i])
      #if(l==4 and m==-4):
       #  print(aa.real,aa.imag,qlm.real,qlm.imag,theta[i],phi[i])
   qlm /= nb
   return qlm
#-------------------END:  qlm--------------------------------------------------------------

#--------------------Bond orientation order: ql--------------------------------------------
def BOO_QL(l,x,y,z):
   import numpy as np
   import scipy.special
   nb=len(x)-1  # neighbour number
   m=[i for i in range(-l,l+1)]
   lm=len(m)
   theta=[0 for i in range(nb)]  #polar angle
   phi=[0 for i in range(nb)]    #azimuthal angle
   for i in range(1,nb+1):
      dx=x[i]-x[0];dy=y[i]-y[0];dz=z[i]-z[0]
      ii=i-1
      if( dz > 0 ) : theta[ii]=np.arctan((np.sqrt(dx**2+dy**2)/dz))
      if( dz == 0 ) : theta[ii]=np.pi/2
      if( dz < 0 ) : theta[ii]=np.arctan((np.sqrt(dx**2+dy**2)/dz)) + np.pi
      if( dx == 0 ) : 
         if( dy > 0 ) : phi[ii]= np.pi/2
         if( dy < 0 ) : phi[ii]= np.pi*3/2
      if( dx>0 and dy>=0 ): phi[ii]=np.arctan(dy/dx)
      if( dx<0 ): phi[ii]=np.arctan(dy/dx) + np.pi
      if( dx>0 and dy<0 ): phi[ii]=np.arctan(dy/dx) + 2*np.pi
   qlm=0
   for i in range(lm):
      #print(l,m[i],abs(BOO_QLM(m[i],l,phi,theta))**2)
      qlm += abs(BOO_QLM(m[i],l,phi,theta))**2
   qlm = np.sqrt(qlm*4*np.pi/(2*l+1))
   return qlm
#-----------------END:Bond orientation order_ql--------------------------------------------

#-----------------Wigner 3j symbol---------------------------------------------------------
def W3J(J1,J,J2,M1,M,M2):
   import numpy as np
   import scipy
   import math
       # set flags to determine allowed 3j symbols based on calling arguments
   if ((abs(J1-J) <= J2) and (J2 <= J1+J)):
      tri = True
   else:
      tri = False

   if ((M1 in np.arange(-J1, J1+1)) and (M in np.arange(-J, J+1)) and \
        (M2 in np.arange(-J2, J2+1))):
      mem = True
   else:
      mem = False

   if (np.floor(J1 + J + J2) == J1 + J + J2):
      perim = True
   else:
      perim = False

   if (M1 + M + M2 == 0):
      mag = True
   else:
      mag = False

    # now compute allowed 3j symbol, return 0 if not allowed

   flag = tri and mem and perim and mag

   if (flag):
      delta = np.sqrt(math.factorial(J1+J2-J)*math.factorial(J1-J2+J)* \
                     math.factorial(-J1+J2+J)/math.factorial(J1+J+J2+1))

      a = np.sqrt(math.factorial(J2+M2)*math.factorial(J2-M2)/math.factorial(J+M)/\
                 math.factorial(J-M)/math.factorial(J1-M-M2)/math.factorial(J1+M+M2))
      s0 = 0

      u  = -J1+J2+J
      z1 = max(-J1-M-M2,u-J2-M2,0)
      z2 = min(J+J2-M-M2,J2-M2,u)

      for z in np.arange(z1,z2+1):
         if ((0 <= J+J2-M-M2-z) and (0 <= J2-M2-z) and \
                (0 <= u-z) and (0 <= J2+M2-u+z)):
            stp_1=math.factorial(z)
            stp_2=math.factorial(J2-M2-z)
            stp_3=math.factorial(u-z)
            stp_4=math.factorial(J2+M2-u+z)
            stp_ll=math.factorial(J+J2-M-M2-z)*math.factorial(J1+M+M2+z)
            stp_ll=stp_ll/stp_1
            stp_ll=stp_ll/stp_2
            stp_ll=stp_ll/stp_3
            stp_ll=stp_ll/stp_4
            s0 = s0 + (-1)**(2*J-J1-M1+z)*stp_ll

    #        s0 = s0 + (-1)**(2*J-J1-M1+z) * \
    #                 math.factorial(J+J2-M-M2-z)*math.factorial(J1+M+M2+z)/\
    #                 math.factorial(z)/math.factorial(J2-M2-z)/math.factorial(u-z)/\
    #                 math.factorial(J2+M2-u+z)
            s = a*delta*s0
         else:
            s = 0
   return s,flag
#-----------------END: Wigner 3j symbol---------------------------------------------------------

#-----------------Bond orientation order: Wl----------------------------------------------------
def BOO_WL(l,x,y,z):
   import numpy as np
   import scipy.special
   nb=len(x)-1  # neighbour number
   m=[i for i in range(-l,l+1)]
   lm=len(m)
   theta=[0 for i in range(nb)]  #polar angle
   phi=[0 for i in range(nb)]    #azimuthal angle
   for i in range(1,nb+1):
      dx=x[i]-x[0];dy=y[i]-y[0];dz=z[i]-z[0]
      ii=i-1
      if( dz > 0 ) : theta[ii]=np.arctan((np.sqrt(dx**2+dy**2)/dz))
      if( dz == 0 ) : theta[ii]=np.pi/2
      if( dz < 0 ) : theta[ii]=np.arctan((np.sqrt(dx**2+dy**2)/dz)) + np.pi
      if( dx == 0 ) :
         if( dy > 0 ) : phi[ii]= np.pi/2
         if( dy < 0 ) : phi[ii]= np.pi*3/2
      if( dx>0 and dy>=0 ): phi[ii]=np.arctan(dy/dx)
      if( dx<0 ): phi[ii]=np.arctan(dy/dx) + np.pi
      if( dx>0 and dy<0 ): phi[ii]=np.arctan(dy/dx) + 2*np.pi
   ego=[]
#   for i in range(0,int(l/2)+1):
#      for j in range(i,l-i+1):
#         c=0-i-j
#         ego.append([i,j,c])
   for i in range(-l,l+1):
      for j in range(-l,l+1):
         c=0-i-j
         if(c>=-l and c<=l): ego.append([i,j,c])
   low=len(ego)
   wl=0
   for i in range(low):
    #  print(W3J(l,l,l,ego[i][0],ego[i][1],ego[i][2]),BOO_QLM(ego[i][0],l,phi,theta),
    #                    BOO_QLM(ego[i][1],l,phi,theta),BOO_QLM(ego[i][2],l,phi,theta))
      eff=W3J(l,l,l,ego[i][0],ego[i][1],ego[i][2])[0]
      wl += eff*BOO_QLM(ego[i][0],l,phi,theta)*\
            BOO_QLM(ego[i][1],l,phi,theta)*BOO_QLM(ego[i][2],l,phi,theta)
   return wl
#-------------------END:Bond orientation order: Wl------------------------------------

#----------------------Bond orientation order: Wl_bar---------------------------------
def BOO_WL_BAR(l,x,y,z):
   import numpy as np
   import scipy.special
   nb=len(x)-1  # neighbour number
   m=[i for i in range(-l,l+1)]
   lm=len(m)
   theta=[0 for i in range(nb)]  #polar angle
   phi=[0 for i in range(nb)]    #azimuthal angle
   for i in range(1,nb+1):
      dx=x[i]-x[0];dy=y[i]-y[0];dz=z[i]-z[0]
      ii=i-1
      if( dz > 0 ) : theta[ii]=np.arctan((np.sqrt(dx**2+dy**2)/dz))
      if( dz == 0 ) : theta[ii]=np.pi/2
      if( dz < 0 ) : theta[ii]=np.arctan((np.sqrt(dx**2+dy**2)/dz)) + np.pi
      if( dx == 0 ) :
         if( dy > 0 ) : phi[ii]= np.pi/2
         if( dy < 0 ) : phi[ii]= np.pi*3/2
      if( dx>0 and dy>=0 ): phi[ii]=np.arctan(dy/dx)
      if( dx<0 ): phi[ii]=np.arctan(dy/dx) + np.pi
      if( dx>0 and dy<0 ): phi[ii]=np.arctan(dy/dx) + 2*np.pi
   qlm=0
   for i in range(lm):
      qlm += abs(BOO_QLM(m[i],l,phi,theta))**2
   qlm=qlm**1.5
   w=BOO_WL(l,x,y,z)
   w_bar=w/qlm
   return w_bar
#----------------------END: Wl_bar---------------------------------------------------

#----------------------Frequency conuts----------------------------------------------
def freq(x,xbin,xstart,xend):
   if(xstart > min(x)):
       print("floor limit!", min(x))
   if(xend < max(x)):
       print("upper limit!", max(x))
   if(xstart==0 and xend == 0):
       xstart=(int(min(x)/xbin)-1)*xbin
       xend=(int(max(x)/xbin)+2)*xbin
   xlen=int(round((xend-xstart)/xbin,5))
   statcs=[0 for i in  range(xlen)]
   for a in x:
      statcs[int(round((a-xstart)/xbin,5))] += 1
   bincenter=[]
   for i in range(xlen):
      bincenter.append(round(xstart+xbin/2+xbin*i,5))
#      print(round(xstart+xbin/2+xbin*i,5),statcs[i])
   return(bincenter,statcs)
#--------------------END: freq-------------------------------------------------------

#--------------------Crystal Elongation---------------------------------------------
def cry_elong(file_in,elongate):
    fin=open(file_in,'r')
    fout=open('SPOSCAR.vasp','w+')
    print('Elongated: ',elongate[0],'*',elongate[1],'*',elongate[2],file=fout)
    fin.readline()
    print(fin.readline().strip('\n'),file=fout)
    v=[[0 for i in range(3)] for j in range(3)]
    for i in range(3):
        line=fin.readline().split()
        for j in range(3): 
            v[i][j]=float(line[j])
            print('%10.4f'%(float(line[j])*elongate[i]),end=' ',file=fout)
        print(' ',file=fout)
    print(fin.readline().strip('\n'),file=fout)
    line=fin.readline().split()
    natom=int(0)
    for i in range(len(line)):
        natom+=int(line[i])
        print('%6d'%(int(line[i])*elongate[0]*elongate[1]*elongate[2]),end=' ',file=fout)
    print(' ',file=fout)
    line=fin.readline().split()
    if(line[0] != 'Direct'):
        print('Use fractional coordinate!')
        return 0
    print('Direct',file=fout)
    atom=[]
    for i in range(natom):
        line=fin.readline().split()
        x=[]
        for j in range(3):
            x.append(float(line[j]))
        atom.append(x)
    for k0 in range(elongate[0]):
        for k1 in range(elongate[1]):
            for k2 in range(elongate[2]):
                for i in range(natom):
                    print('%10.4f%10.4f%10.4f'%((atom[i][0]+k0)/elongate[0],(atom[i][1]+k1)/elongate[1],(atom[i][2]+k2)/elongate[2]),file=fout)
    return 1
    fin.close()
    fout.close()
#--------------------END: Crystal Elongation---------------------------------------------


#--------------------dump file to poscar, average----------------------------------------- 
def lammps2posave(filein,atom_type):
    fin=open(filein,'r')
    for i in range(3): fin.readline()
    natom=int(fin.readline())
    line=fin.readline()
# vector
    steps=0
    v=[[0 for i in range(3)] for j in range(3)]
    while(line !=''):
        l1=fin.readline().split()
        l2=fin.readline().split()
        l3=fin.readline().split()
        if(len(l1)==2) :
            v[0][0] += float(l1[1])-float(l1[0])
            v[1][1] += float(l2[1])-float(l2[0])
            v[2][2] += float(l3[1])-float(l3[0])
            steps += 1
        else :
            v[0][0] += float(l1[1])-float(l1[0])
            v[1][1] += float(l2[1])-float(l2[0])
            v[2][2] += float(l3[1])-float(l3[0])
            v[1][0] += float(l1[2])
            v[2][0] += float(l2[2])
            v[2][1] += float(l3[2])
            steps += 1
        for i in range(natom+6): line=fin.readline()
    if(len(l1)==2):
        v[0][0] /= steps
        v[1][1] /= steps
        v[2][2] /= steps
    else : 
        v[0][0] /= steps
        v[1][1] /= steps
        v[2][2] /= steps
        v[1][0] /= steps
        v[2][0] /= steps
        v[2][1] /= steps
    fin.close()
    fin=open(filein,'r')
    fout=open('snapave.vasp','w')
    ntype=len(atom_type)
    print(filein,' -> averaged snapshot',steps, "steps.", file=fout)
    print('%5.3f'%(1.0),file=fout)
   # output vector
    for i in range(3):  print('%16.6f%16.6f%16.6f'%(v[i][0],v[i][1],v[i][2]),file=fout)
   # output atom species
    for i in range(ntype): print('%5s'%(atom_type[i]),end=' ',file=fout)
    print(' ',file=fout)

    for i in range(8): fin.readline()
    # xyz type ......
    arg_id=[];arg_type=[];arg_xs=[];arg_ys=[];arg_zs=[];arg_x=[];arg_y=[];arg_z=[]
    arg_fx=[];arg_fy=[];arg_fz=[];arg_cpa=[];arg_vx=[];arg_vy=[];arg_vz=[];arg_QConn=[]
 #  dic={arg_id:'id',arg_type:'type',arg_xs:'xs',arg_ys:'ys',arg_zs:'zs',arg_x:'x',arg_y:'y',arg_z:'z'}
    dic={'id':arg_id,'type':arg_type,'xs':arg_xs,'ys':arg_ys,'zs':arg_zs,'x':arg_x,'y':arg_y,'z':arg_z,
          'fx':arg_fx,'fy':arg_fy,'fz':arg_fz,'c_peratom':arg_cpa, 'c_pperatom':arg_cpa,'c_pea':arg_cpa,
            'vx':arg_vx, 'vy':arg_vy, 'vz':arg_vz, 'c_QConn':arg_QConn}
    l_arg=fin.readline().split()
    arg_num=len(l_arg)-2
    for i in range(natom):
       line=fin.readline().split()
       for j in range(arg_num):
          dic[l_arg[2+j]].append(line[j])
    x0=[0 for i in range(natom)];  xave=[0 for i in range(natom)]
    y0=[0 for i in range(natom)];  yave=[0 for i in range(natom)]
    z0=[0 for i in range(natom)];  zave=[0 for i in range(natom)]
    epave=[0 for i in range(natom)]
    for i in range(natom):
        x0[i]=float(arg_xs[i])
        y0[i]=float(arg_ys[i])
        z0[i]=float(arg_zs[i])
        xave[i]=float(arg_xs[i])
        yave[i]=float(arg_ys[i])
        zave[i]=float(arg_zs[i])
        if(len(arg_cpa) != 0):
            epave[i]=float(arg_cpa[i])
    n=[0 for i in range(ntype)]
    for i in range(ntype):  n[i]=arg_type.count(str(i+1))
   # output atom numbers for each species
    for i in range(ntype): print('%5d'%(n[i]),end=' ',file=fout)
    print(' ',file=fout)
    print('Direct',file=fout)
    if(len(arg_xs) != 0 ):
        for ii in range(steps-1):
            for jj in range(9):fin.readline()
            for i in range(natom):
                line=fin.readline().split()
                for j in range(arg_num):
                    dic[l_arg[2+j]][i]=line[j]
                if(float(arg_xs[i])-x0[i]>0.5):
                    xave[i] += float(arg_xs[i])-1
                elif(float(arg_xs[i])-x0[i]<-0.5):
                    xave[i] += float(arg_xs[i])+1
                else:
                    xave[i] += float(arg_xs[i])
                if(float(arg_ys[i])-y0[i]>0.5):
                    yave[i] += float(arg_ys[i])-1
                elif(float(arg_ys[i])-y0[i]<-0.5):
                    yave[i] += float(arg_ys[i])+1
                else:
                    yave[i] += float(arg_ys[i])
                if(float(arg_zs[i])-z0[i]>0.5):
                    zave[i] += float(arg_zs[i])-1
                elif(float(arg_zs[i])-z0[i]<-0.5):
                    zave[i] += float(arg_zs[i])+1
                else:
                    zave[i] += float(arg_zs[i])
                if(len(arg_cpa) != 0):
                    epave[i]+=float(arg_cpa[i])
        for k in range(ntype):
            for i in range(natom):
                if(arg_type[i]==str(k+1)):
                    if(len(arg_cpa) == 0):
                        print("%10.6f%10.6f%10.6f%10i"%(xave[i]/steps,yave[i]/steps,zave[i]/steps,int(arg_id[i])),file=fout)
                    else:
                        print("%10.6f%10.6f%10.6f%10i%10.6f"%(xave[i]/steps,yave[i]/steps,zave[i]/steps,int(arg_id[i]),epave[i]/steps),file=fout)
        return 1
    else:
        return ('use xs in lammps input')

#------crystal_takecluster---------------------------------------------
def cry_cluster(file_in, name_id, ntemp,center, neighbor,rrcut):
    import numpy as np
    import os
    if not os.path.exists("clusters"):
        os.makedirs("clusters")
    fin=open(file_in,'r')
    rcut=rrcut # angstrom
    fin.readline()
    factor=float(fin.readline().split()[0])
    v=[[0 for i in range(3)] for j in range(3)]
    elong=[0 for i in range(3)]
    for i in range(3):
        line=fin.readline().split()
        vl=0
        for j in range(3): 
            v[i][j]=float(line[j])*factor
#        vl += v[i][j]**2
#        vl=np.sqrt(vl)
#        #elong[i] = int(rcut/vl) + 1
#        elong[i] = int(rcut/v[i][i]) + 1
    vol=volume(v[0],v[1],v[2])
    b1=AcrossB(v[1],v[2])
    b2=AcrossB(v[2],v[1])
    b3=AcrossB(v[0],v[1])
    recv=[[0 for i in range(3)] for j in range(3)]
    for i in range(3):
        recv[0][i] = b1[i]/vol
        recv[1][i] = b2[i]/vol
        recv[2][i] = b3[i]/vol
    rvsq1= AdotB(recv[0],recv[0])
    rvsq2= AdotB(recv[1],recv[1])
    rvsq3= AdotB(recv[2],recv[2])

    elong[0] = int(rcut*np.sqrt(rvsq1))+1
    elong[1] = int(rcut*np.sqrt(rvsq2))+1
    elong[2] = int(rcut*np.sqrt(rvsq3))+1


    atype=fin.readline().split()
    ntype=len(atype)
    ll=fin.readline().split()
    na=[]
    for i in ll:
        na.append(int(i))
    tpe=[]
    for i in range(ntype):
        for k in range(na[i]):
            tpe.append( atype[i] )
    
    natom=sum(na)
    ntotal=natom*(2*elong[0]+1)*(2*elong[1]+1)*(2*elong[2]+1)
    if(fin.readline().split()[0] != "Direct"):
        print("Fractional coordinate!")
        exit()

    x=[]
    for i in range(ntype):
        for k in range(na[i]):
            ll=fin.readline().split()
            x.append([float(ll[0]),float(ll[1]), float(ll[2])])
    fin.close()

    # make elongations
    xall=[]
    all_type=[]
    all_id=[] # id in the primitive cell, starting from 1
    for k0 in range(int(-1)*elong[0],elong[0]+1):
        for k1 in range(int(-1)*elong[1],elong[1]+1):
            for k2 in range(int(-1)*elong[2],elong[2]+1):
                for ia in range(natom):
                    xall.append( [ x[ia][0] + k0, x[ia][1] + k1, x[ia][2] + k2 ] )
                    all_type.append(tpe[ia])
                    all_id.append(ia+1)
    if(len(xall) != ntotal): print("Warning!")
   
    print("elongated by" ,*elong,end=" ... " )
    # takecluster
    fout=open("clusters/cluster-total.xyz","w+")
    count=0
    for i in range(natom):
        if( tpe[i] != center ):
           continue
        count+=1
       
        print("%20d"%(ntemp), file=fout)
        print("%20s%10d"%(file_in,count), file=fout)
        ic=0
        nei_xyz=[]
        nei_dis=[]
        nei_tpe=[]
        nei_id=[]
        for j in range(ntotal):
            dd=[]
            for k in range(3):
                dd.append( (xall[j][0]-x[i][0])*v[0][k] + (xall[j][1]-x[i][1])*v[1][k] + (xall[j][2]-x[i][2])*v[2][k] )
            dis= np.sqrt( dd[0]**2 + dd[1]**2 + dd[2]**2)
            if( dis < rcut):
                ic += 1
                nei_xyz.append(dd)
                nei_dis.append(dis)
                nei_tpe.append(all_type[j])
                nei_id.append(all_id[j])
        # pai xu
        for ki in range(ic-1):
            for kj in range(ki+1,ic):
                if(nei_dis[ki] > nei_dis[kj]):
                    exchange=nei_xyz[ki]
                    nei_xyz[ki]=nei_xyz[kj]
                    nei_xyz[kj]=exchange

                    exchange=nei_dis[ki]
                    nei_dis[ki]=nei_dis[kj]
                    nei_dis[kj]=exchange
                    
                    exchange=nei_tpe[ki]
                    nei_tpe[ki]=nei_tpe[kj]
                    nei_tpe[kj]=exchange
        
                    exchange=nei_id[ki]
                    nei_id[ki]=nei_id[kj]
                    nei_id[kj]=exchange
        fsp=open("clusters/POSCAR_"+name_id+"-"+str(count)+"-"+nei_tpe[0]+".xyz","w+")
        print("%20d"%(ntemp), file=fsp)
        print("%20s%10d"%(file_in, count), file=fsp)
        j=0; nb_ct=0
        while nb_ct<ntemp:
            if(nei_tpe[j] in neighbor):
                print("%2s%12.6f%12.6f%12.6f"%(nei_tpe[j], nei_xyz[j][0], nei_xyz[j][1], nei_xyz[j][2]), nei_id[j],file=fout)
                print("%2s%12.6f%12.6f%12.6f"%(nei_tpe[j], nei_xyz[j][0], nei_xyz[j][1], nei_xyz[j][2]), nei_id[j],file=fsp)
                j+=1; nb_ct+=1
            else:
                j+=1
            if(j>=ntotal):
               print("Not enough neighbor atoms, enlarge cell!")
               return 0
        fsp.close()
    fout.close()
    return count
#--------------------END: Crystal Takecluster---------------------------------------------

#--------------------xdat2dump------------------------------------------------------------
def xdat2dump(file_in,nstep):
    print("Note: only for direct coordinate, only used boxsize from first step")
    fin=open(file_in,"r")
    fin.readline()
    scaler=float(fin.readline().split()[0])
    box=[]
    xy=0; yz=0; xz=0
    for k in range(3):
        ll=fin.readline().split()
        if(k==0):
           if( (float(ll[1]) != 0) or (float(ll[2]) != 0)):
              print("Error: lattice not triangle!")
        if(k==1):
           if( float(ll[0]) != 0 ):
              xy=float(ll[0])*scaler
           if( float(ll[2]) != 0 ):
              print("Error: lattice not triangle!")
        if(k==2):
           if( float(ll[0]) != 0 ):
              xz = float(ll[0])*scaler
           if( float(ll[1]) != 0 ):
              yz = float(ll[1])*scaler
        box.append(float(ll[k]))

    for i in range(len(box)):
        box[i]*=scaler
    fin.readline()
    ll=fin.readline().split()
    ntype=len(ll)
    na=[ int(ll[k]) for k in range(ntype)]
    fout=open("xdat2dump.atom","w+")
    for istep in range(nstep):
        fin.readline()
        print("ITEM: TIMESTEP",file=fout)
        print(istep,file=fout)
        print("ITEM: NUMBER OF ATOMS",file=fout)
        print(sum(na),file=fout)
        print("ITEM: BOX BOUNDS xy xz yz pp pp pp",file=fout)
        print(0,box[0]+xy+xz,xy,file=fout)
        print(0,box[1]+yz,xz,file=fout)
        print(0,box[2],yz,file=fout)
        print("ITEM: ATOMS id type xs ys zs",file=fout)
        count = 1
        flag=0
        for itype in range(ntype):
            for iatom in range(na[itype]):
                ll=fin.readline().split()
                if(len(ll) > 3):
                    print(ll[3], itype+1, ll[0], ll[1],ll[2],file=fout)
                    flag=1
                else:
                    print(count, itype+1, ll[0], ll[1], ll[2], file=fout)
                count += 1
    if(flag==1):
        print("Colume 4 is used for the atom ID.")
    else:
        print("Atom ID is automatically generated!")
    fin.close()
    fout.close()
    return 1
#---------------END: xdat2dump------------------------------------------------------------

#--------lammps2xml--------#
def lammps2xml(filein,step,atom_type):
   fin=open(filein,'r')
   for i in range(step-1):
      for k in range(3):
         fin.readline()
      natom=int(fin.readline())
      for k in range(5+natom):
         fin.readline()
   fin.readline()
   atep=int(fin.readline().strip('\n'))
   fin.readline()
   natom=int(fin.readline())
   fin.readline()
   print(atep)
# vector
   v=[[0 for i in range(3)] for j in range(3)]
   vlow=[0 for i in range(3)]
   vhigh=[ 0 for i in range(3)]
   l1=fin.readline().split()
   l2=fin.readline().split()
   l3=fin.readline().split()
   vlow[0]=float(l1[0]); vlow[1]=float(l2[0]); vlow[2]=float(l3[0])
   vhigh[0]=float(l1[1]); vhigh[1]=float(l2[1]); vhigh[2]=float(l3[1])
   if(len(l1)==2) :
      v[0][0]=float(l1[1])-float(l1[0])
      v[1][1]=float(l2[1])-float(l2[0])
      v[2][2]=float(l3[1])-float(l3[0])
   else :
      v[0][0]=float(l1[1])-float(l1[0])
      v[1][1]=float(l2[1])-float(l2[0])
      v[2][2]=float(l3[1])-float(l3[0])
      v[1][0]=float(l1[2])
      v[2][0]=float(l2[2])
      v[2][1]=float(l3[2])
# xyz type ......
   arg_id=[];arg_type=[];arg_xs=[];arg_ys=[];arg_zs=[];arg_x=[];arg_y=[];arg_z=[]
   arg_fx=[];arg_fy=[];arg_fz=[];arg_cpa=[]; arg_vx=[]; arg_vy=[]; arg_vz=[];arg_cpa2=[];
   arg_cka=[]; arg_diff=[]; arg_nb=[]; arg_V=[]; arg_sro=[]; arg_rmsd=[]
   dic={'id':arg_id,'type':arg_type,'xs':arg_xs,'ys':arg_ys,'zs':arg_zs,'x':arg_x,'y':arg_y,'z':arg_z,
           'fx':arg_fx,'fy':arg_fy,'fz':arg_fz,'c_peratom':arg_cpa, 'vx':arg_vx, 'vy':arg_vy, 'vz':arg_vz,
           'c_pea':arg_cpa2,'c_kea':arg_cka, 'diffusion':arg_diff, 'nb':arg_nb, 'V':arg_V, 'sro':arg_sro,
           'RMSD':arg_rmsd}
   l_arg=fin.readline().split()
   arg_num=len(l_arg)-2
   for i in range(natom):
      line=fin.readline().split()
      for j in range(arg_num):
         dic[l_arg[2+j]].append(line[j])
# sort
   if( len( arg_id) != 0):
      for k in range(natom):
         arg_id[k] =  int(arg_id[k])
      import numpy
      sort=numpy.argsort(arg_id)
   else:
      sort=[k for k in range(natom)]
# output
   if(len(arg_type)==0): 
      print('Error: Include "type" in your dump commond!')
      return 0
 #  ntype=int(arg_type[natom-1])
   ntype=len(atom_type)
   n=[0 for i in range(ntype)]
   for i in range(ntype):  n[i]=arg_type.count(str(i+1))
   fout=open(str(atep)+'.xml','w')
   print('<?xml version="1.0" encoding="UTF-8"?>',file=fout)
   print('<hoomd_xml version="1.7">',file=fout)
   print('<configuration time_step="'+str(atep)+'" dimensions="3" natoms="'+str(natom)+'" >', file=fout)
   print('<box lx="'+str(v[0][0])+'" ly="'+str(v[1][1])+'" lz="'+str(v[2][2])+'"',end=" ",file=fout)
   if(len(l1)==2):
      print('xy="0" xz="0" yz="0"/>',file=fout)
   else:
      print('xy="'+str(v[1][0])+' xz="'+str(v[2][0])+' yz="'+str(v[2][1])+'"/>',file=fout)
   # output vector
   print('<position num="'+str(natom)+'">',file=fout)
   if(len( arg_x) != 0):
      for k in range(natom):
         i=sort[k]
         print(arg_x[i],arg_y[i], arg_z[i], file=fout)
   else :
      for k in range(natom):
         i=sort[k]
         xx=vlow[0]+float(arg_xs[i])*v[0][0]+float(arg_ys[i])*v[1][0]+float(arg_zs[i])*v[2][0] - (vlow[0]+vhigh[0])/2
         yy=vlow[1]+float(arg_xs[i])*v[0][1]+float(arg_ys[i])*v[1][1]+float(arg_zs[i])*v[2][1] - (vlow[1]+vhigh[1])/2
         zz=vlow[2]+float(arg_xs[i])*v[0][2]+float(arg_ys[i])*v[1][2]+float(arg_zs[i])*v[2][2] - (vlow[2]+vhigh[2])/2
         print(xx,yy,zz,file=fout)
   print('</position>',file=fout)
   factor=0.01*1.01805
   if (len( arg_vx) !=0):
      print('<velocity num="'+str(natom)+'">',file=fout)
      for k in range(natom):
         i=sort[k]
         print(factor*float(arg_vx[i]),factor*float(arg_vy[i]), factor*float(arg_vz[i]), file=fout)
      print('</velocity>',file=fout)
   element_mass={'Co':58.933,'Zr':91.224,'Hf':178.490,'Al':26.983,'Sm':150.360,'Cu':63.546,'Ni':58.71}
   print('<mass num="'+str(natom)+'">',file=fout)
   for k in range(natom):
      i=sort[k]
      print(element_mass[atom_type[int(arg_type[i])-1]],file=fout)
   print('</mass>',file=fout)
   print('<type num="'+str(natom)+'">',file=fout)
   for k in range(natom):
      i=sort[k]
      print(atom_type[int(arg_type[i])-1],file=fout)
   print('</type>',file=fout)
   print('</configuration>\n</hoomd_xml>',file=fout)

#--------------------average dump file----------------------------------------- 
def dumpave(filein):
    fin=open(filein,'r')
    for i in range(3): fin.readline()
    natom=int(fin.readline())
    line=fin.readline()
# vector
    steps=0
    v=[[0 for i in range(3)] for j in range(3)]
    while(line !=''):
        l1=fin.readline().split()
        l2=fin.readline().split()
        l3=fin.readline().split()
        v[0][0] += float(l1[0])
        v[0][1] += float(l1[1])
        v[1][0] += float(l2[0])
        v[1][1] += float(l2[1])
        v[2][0] += float(l3[0])
        v[2][1] += float(l3[1])
        steps += 1
        for i in range(natom+6): line=fin.readline()
    v[0][0] /= steps
    v[0][1] /= steps
    v[1][0] /= steps
    v[1][1] /= steps
    v[2][0] /= steps
    v[2][1] /= steps
    fin.seek(0)
    fout=open('average.atom','w')
    print("ITEM: TIMESTEP\n0\nITEM: NUMBER OF ATOMS",file=fout)
    print(natom,file=fout)
    print("ITEM: BOX BOUNDS pp pp pp",file=fout)
    for i in range(3):  
       print('%16.6f%16.6f'%(v[i][0],v[i][1]),file=fout)
    print("ITEM: ATOMS id type xs ys zs",file=fout)
    for i in range(8): fin.readline()
    # xyz type ......
    arg_id=[];arg_type=[];arg_xs=[];arg_ys=[];arg_zs=[];arg_x=[];arg_y=[];arg_z=[]
    arg_fx=[];arg_fy=[];arg_fz=[];arg_cpa=[];arg_vx=[];arg_vy=[];arg_vz=[]
    dic={'id':arg_id,'type':arg_type,'xs':arg_xs,'ys':arg_ys,'zs':arg_zs,'x':arg_x,'y':arg_y,'z':arg_z,
            'fx':arg_fx,'fy':arg_fy,'fz':arg_fz,'c_peratom':arg_cpa, 'c_pperatom':arg_cpa,'vx':arg_vx, 'vy':arg_vy, 'vz':arg_vz}
    l_arg=fin.readline().split()
    arg_num=len(l_arg)-2
    for i in range(natom):
       line=fin.readline().split()
       for j in range(arg_num):
          dic[l_arg[2+j]].append(line[j])
    x0=[0 for i in range(natom)];  xave=[0 for i in range(natom)]
    y0=[0 for i in range(natom)];  yave=[0 for i in range(natom)]
    z0=[0 for i in range(natom)];  zave=[0 for i in range(natom)]
    epave=[0 for i in range(natom)]
    for i in range(natom):
        x0[i]=float(arg_xs[i])
        y0[i]=float(arg_ys[i])
        z0[i]=float(arg_zs[i])
        xave[i]=float(arg_xs[i])
        yave[i]=float(arg_ys[i])
        zave[i]=float(arg_zs[i])
        if(len(arg_cpa) != 0):
            epave[i]=float(arg_cpa[i])
    if(len(arg_xs) != 0 ):
        for ii in range(steps-1):
            for jj in range(9):fin.readline()
            for i in range(natom):
                line=fin.readline().split()
                for j in range(arg_num):
                    dic[l_arg[2+j]][i]=line[j]
                if(float(arg_xs[i])-x0[i]>0.5):
                    xave[i] += float(arg_xs[i])-1
                elif(float(arg_xs[i])-x0[i]<-0.5):
                    xave[i] += float(arg_xs[i])+1
                else:
                    xave[i] += float(arg_xs[i])
                if(float(arg_ys[i])-y0[i]>0.5):
                    yave[i] += float(arg_ys[i])-1
                elif(float(arg_ys[i])-y0[i]<-0.5):
                    yave[i] += float(arg_ys[i])+1
                else:
                    yave[i] += float(arg_ys[i])
                if(float(arg_zs[i])-z0[i]>0.5):
                    zave[i] += float(arg_zs[i])-1
                elif(float(arg_zs[i])-z0[i]<-0.5):
                    zave[i] += float(arg_zs[i])+1
                else:
                    zave[i] += float(arg_zs[i])
                if(len(arg_cpa) != 0):
                    epave[i]+=float(arg_cpa[i])
        for i in range(natom):
            xout=xave[i]/steps
            if(xout > 1): xout -= 1
            if(xout < 0): xout += 1
            yout=yave[i]/steps
            if(yout > 1): yout -= 1
            if(yout < 0): yout += 1
            zout=zave[i]/steps
            if(zout > 1): zout -= 1
            if(zout < 0): zout += 1
            print(arg_id[i],arg_type[i],xout, yout, zout, file=fout)
        return 1
    else:
        return ('use xs in lammps input')


##---------draw_ternary------------------------------------------------------##
def draw_ternary(pts,ele,string):
    #https://github.com/marcharper/python-ternary
    import matplotlib
    import ternary
    matplotlib.rcParams['figure.dpi'] = 200
    matplotlib.rcParams['figure.figsize'] = (4, 4)

    # scales
    figure, tax = ternary.figure(scale=1.0)
    # boundary
    tax.boundary(linewidth=1.0)
    tax.gridlines(color="grey", multiple=0.1)
    #tax.gridlines(color="blue", multiple=2, linewidth=0.5)

    # labels and title                                                                 57
    fontsize=12
    #tax.set_title("ternary phase diagram\n", fontsize=fontsize)
    tax.set_title(" ", fontsize=fontsize)
    tax.left_axis_label("$x_C$", fontsize=fontsize, offset=0.14)
    tax.right_axis_label("$x_B$", fontsize=fontsize, offset=0.14)
    tax.bottom_axis_label("$x_A$", fontsize=fontsize, offset=0.14)

    # ticks
    tax.ticks(axis='lbr', linewidth=1, multiple=0.2, tick_formats="%.1f", offset=0.03)

    # text
    figure.text(.89,0.12,ele[0],fontsize=12)
    figure.text(.48,0.93,ele[1],fontsize=12)
    figure.text(.02,0.15,ele[2],fontsize=12)

    # plot data
    pdata=[]
    for pt in pts:
        mm=pt[:3]
        pdata.append(1.0*mm/sum(mm))
    tax.scatter(pdata, color='red', s=8.0)

    # remove matplotlib axes
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')
    #figure.tight_layout()

#    ternary.plt.show()
    ternary.plt.savefig(string+".png")


##---------draw_ternary_convex-------------------------------------------------##
# input: pts as list with element [nA,nB,nC,Ef]
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
def draw_ternary_convex(pts,ele,string):
    import matplotlib
    import ternary
    import numpy as np
    from scipy.spatial import ConvexHull
    matplotlib.rcParams['figure.dpi'] = 200
    matplotlib.rcParams['figure.figsize'] = (4, 4)

    # scales
    figure, tax = ternary.figure(scale=1.0)
    # boundary
    tax.boundary(linewidth=0.5)
    tax.gridlines(color="grey", multiple=0.1)
    #tax.gridlines(color="blue", multiple=2, linewidth=0.5)

    # labels and title                                                                 57
    fontsize=12
    #tax.set_title("ternary phase diagram\n", fontsize=fontsize)
    #tax.set_title(" ", fontsize=fontsize)
    tax.left_axis_label(  "$x_{"+str(ele[2])+"}$", fontsize=fontsize, offset=0.14)
    tax.right_axis_label( "$x_{"+str(ele[1])+"}$", fontsize=fontsize, offset=0.14)
    tax.bottom_axis_label("$x_{"+str(ele[0])+"}$", fontsize=fontsize, offset=0.14)
    tax.right_corner_label(ele[0], fontsize=fontsize+2)
    tax.top_corner_label(ele[1],   fontsize=fontsize+2)
    tax.left_corner_label(ele[2],  fontsize=fontsize+2)

    # ticks
    tax.ticks(axis='lbr', linewidth=1, multiple=0.2, tick_formats="%.1f", offset=0.03)

    # text
    #figure.text(.91,0.12, ele[0], fontsize=14)
    #figure.text(.5,0.93 , ele[1], fontsize=14)
    #figure.text(.05,0.14, ele[2], fontsize=14)
    
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
    print("# of stable structures",len(hull.vertices),":")
    print(*ele,"Ef(eV/atom)")
    stables=[]
    stable_info=""
    for iv in hull.vertices:
        print("%2d%3d%3d%12.6f"%(int(pts[iv][0]),int(pts[iv][1]),int(pts[iv][2]),pts[iv][3]))
        comp_info=""
        zero_tag=0
        for i in range(3):
            if(int(pts[iv][i])==0):
                zero_tag+=1
                continue
            elif(int(pts[iv][i])==1):
                comp_info+=ele[i]
            else:
                comp_info+=ele[i]+"$_{"+str(int(pts[iv][i]))+"}$"
        if(zero_tag!=2):
            stable_info+=comp_info+"\n"  # ignore pure elements

        name=ele[0]+str(int(pts[iv][0]))+ele[1]+str(int(pts[iv][1]))+ele[2]+str(int(pts[iv][2]))
        stables.append([tpts[iv][0],tpts[iv][1],name])  # still not sure how to plot names on the figure 06/24
    #print(stables[0])
    #tax.text(stables[0][0], stables[0][1], stables[0][2], fontsize=14)
    #tax.set_title(stable_info, fontsize=int(fontsize*0.5), loc='left')

    # plot data
    pdata=[]
    for pt in pts:
        mm=pt[:3]
        pdata.append(1.0*mm/sum(mm))
    #1 plot stable and connect them
    label_tag=0
    for isimp in hull.simplices:
        if(label_tag==0):
            #tax.line(pdata[isimp[0]],pdata[isimp[1]],linewidth=1.,marker='.',markersize=6.,color='black', label=stable_info[:-1])
            tax.scatter([pdata[isimp[0]],pdata[isimp[1]]],color='k', marker='.', s=15.0, linewidth=1.5,label=stable_info[:-1])
        tax.line(pdata[isimp[0]],pdata[isimp[1]],linewidth=1.,marker='.',markersize=6.,color='black')
        tax.line(pdata[isimp[0]],pdata[isimp[2]],linewidth=1.,marker='.',markersize=6.,color='black')
        tax.line(pdata[isimp[1]],pdata[isimp[2]],linewidth=1.,marker='.',markersize=6.,color='black')
        label_tag=1
    
    
    #2 get meta-stable phases
    mstables=[]
    for i in range(len(pdata)):
        if(i not in hull.vertices):
            mstables.append(pdata[i])
    
    #3 find the distance to the convex hull
    print("# of metastable structures",len(mstables),":")
    print(*ele,"Ef(eV/atom) E_to_convex_hull(eV/atom)")
    mst_info=""
    for k in range(len(tpts)):
        if(k in hull.vertices):
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
                heights.append(h) # same data for the one at edges
        h=max(heights)
        print("%2d%3d%3d%12.6f%11.6f"%(int(pts[k][0]),int(pts[k][1]),int(pts[k][2]), pts[k][3],h))
        comp_info=""
        for i in range(3):
            if(int(pts[k][i])==0):
                continue
            elif(int(pts[k][i])==1):
                comp_info+=ele[i]
            else:
                comp_info+=ele[i]+"$_{"+str(int(pts[k][i]))+"}$"
        mst_info+=comp_info+",%.3f"%(h)+"eV\n"  # ignore pure elements

    #4 plot meta-stable phases
    mstables=np.array(mstables)
    if(len(mstables) != 0):
        tax.scatter(mstables, color='r', marker='x', s=20.0, linewidth=1.5,label=mst_info[:-1])

    # remove matplotlib axes
    tax.clear_matplotlib_ticks()
    tax.get_axes().axis('off')
    #tax.legend(fontsize=int(fontsize*0.5),bbox_to_anchor=(0.2,0.8), frameon=False)
    tax.legend(fontsize=int(fontsize*0.5),loc='upper left', frameon=False)

    ternary.plt.savefig(string+".png",bbox_inches = "tight")
    #ternary.plt.savefig(string+".png")
    #ternary.plt.show()

def vasp2spin(xml_path):
    # output, type:list, data: E, total[up, down], s[up,down], p[up, down], d[up, down]
    import xml.etree.ElementTree as ET
    import numpy as np
    tree =  ET.parse(xml_path)
    root = tree.getroot()
    # atom number
    natom=int(root.find('.//atominfo/atoms').text)
    # Efermi
    Efermi=float(root.find('.//calculation/dos/i[@name="{}"]'.format('efermi')).text)
    # total spin
    spin =  root.find('.//calculation/dos/total/array/set/set[@comment="{}"]'.format('spin 1'))
    E=[]; t_up=[]; t_dn=[]
    for ele in spin:
        ll=ele.text.split()
        E.append(float(ll[0])-Efermi)
        t_up.append(float(ll[1]))
    spin =  root.find('.//calculation/dos/total/array/set/set[@comment="{}"]'.format('spin 2'))
    for ele in spin:
        ll=ele.text.split()
        t_dn.append(float(ll[1]))
    # partial spin 
    su=np.array([0.0 for i in range(len(E))]); sd=np.array([0.0 for i in range(len(E))]); 
    pu=np.array([0.0 for i in range(len(E))]); pd=np.array([0.0 for i in range(len(E))]); 
    du=np.array([0.0 for i in range(len(E))]); dd=np.array([0.0 for i in range(len(E))]); 
    for i in range(natom):
        spin =  root.find('.//calculation/dos/partial/array/set/set[@comment="{}"]'.format('ion '+str(i+1)))
        # this finds "spin 1" and "spin 2" for "ion 1"
        # E, s, py, pz, px, dxy, dyz, dz2, dxz, dx2
        for j in range(len(spin[0])):
            l1=spin[0][j].text.split()
            l2=spin[1][j].text.split()
            # s 
            su[j]+=float(l1[1])
            sd[j]+=float(l2[1])
            # p
            pu[j]+=float(l1[2])+float(l1[3])+float(l1[4])
            pd[j]+=float(l2[2])+float(l2[3])+float(l2[4])
            # d
            du[j]+=float(l1[5])+float(l1[6])+float(l1[7])+float(l1[8])+float(l1[9])
            dd[j]+=float(l2[5])+float(l2[6])+float(l2[7])+float(l2[8])+float(l2[9])
    # unit to f.u.
    Tup=np.array(t_up)
    Tdn=np.array(t_dn)
    Tup = Tup/natom 
    Tdn = Tdn/natom*(-1) 
    su = su/natom
    sd = sd/natom*(-1)
    pu = pu/natom
    pd = pd/natom*(-1)
    du = du/natom
    dd = dd/natom*(-1)
    
    total=[]; s=[]; p=[]; d=[]
    total.append(Tup.tolist())
    total.append(Tdn.tolist())
    s.append(su.tolist())
    s.append(sd.tolist())
    p.append(pu.tolist())
    p.append(pd.tolist())
    d.append(du.tolist())
    d.append(dd.tolist())
    
    print("natom:", natom, "; Efermi:", Efermi)
    print("units: E(eV), DOS(/atom)")
    return E, total, s, p, d

##---------draw_binary_convex(pts,ele,string)-------------------------------##
def draw_binary_convex(pts,ele,string):
    from scipy.spatial import ConvexHull
    import numpy as np
    import pylab as plt
    from matplotlib.ticker import MultipleLocator
   
    apts=[]; apts_names=[]; 
    unstables=[]; unstables_names=[]
    idx=0
    for ipt in pts:
        if(ipt[2]<=0):
            apts.append([ipt[1]*1.0/(ipt[0]+ipt[1]) ,ipt[2]])
            apts_names.append(ele[0]+str(int(ipt[0]))+ele[1]+str(int(ipt[1])))
            idx+=1
        else:
            unstables.append([ipt[1]/(ipt[0]+ipt[1]) ,ipt[2]])
            unstables_names.append(ele[0]+str(int(ipt[0]))+ele[1]+str(int(ipt[1])))
    print(apts)
    # get stables
    hull=ConvexHull(apts)
    stables=[]
    stables_names=[]
    print("stable structures:",len(hull.vertices))
    for iv in hull.vertices:
#        stables.append(apts[iv])
#        stables_names.append(apts_names[iv])
        stables.append([apts[iv][0],apts[iv][1],apts_names[iv]])
        print(pts[iv])
    stables_sort=sorted(stables, key=lambda l:l[0],reverse=False)
    for vvv in stables_sort:
        print(vvv)
    # get metastables
    mstables=[]
    mstables_names=[]
    #print("metastable:")
    for i in range(len(apts)):
        if(i not in hull.vertices):
            mstables.append(apts[i])
            mstables_names.append(apts_names[iv])
        #    print(apts[i])
    # draw
    plt.figure(figsize=(8,6))
    mstables_plt=np.array(mstables)
    unstables_plt=np.array(unstables)
    # meta_stable
    if(len(mstables)!=0):
        plt.plot(mstables_plt[:,0],mstables_plt[:,1],"x", color="g", markersize=6)
        #for kkk in range(len(mstables_names)):
        #    plt.annotate(mstables_names[kkk], ( mstables_plt[kkk,0], mstables_plt[kkk,1] ), size=16)
    # unstable
    if(len(unstables)!=0):
        plt.plot(unstables_plt[:,0],unstables_plt[:,1],"x", color='r',markersize=6)
        #for kkk in range(len(unstables_names)):
        #    plt.annotate(unstables_names[kkk], ( unstables_plt[kkk,0], unstables_plt[kkk,1] ), size=12)
    # stables
    stables_plt=np.array(stables_sort,dtype=object)
    plt.plot(stables_plt[:,0],stables_plt[:,1],"o-", color='k', markersize=6)
    #print(len(stables))
    for kkk in range(len(stables)):
        #idx=int(stables_plt[kkk,2]+0.5)
        #print(idx)
        #plt.annotate(stables_names[idx], ( stables_plt[kkk,0], stables_plt[kkk,1] ), size=12)
        plt.annotate(stables_plt[kkk,2], ( stables_plt[kkk,0], stables_plt[kkk,1] ), size=12)
    # dash line 0 to 1
    dash=np.array([[0,0],[1,0]])
    plt.plot(dash[:,0], dash[:,1], '--', color='k')
    
    fc=16
    xr = [k*0.1 for k  in range(0,11,1)]
    ir=np.array(xr)
    xr = ["%.1f"%(k*0.1) for k  in range(0,11,1)]
    xr[0]=ele[0]; xr[-1]=ele[1]
    plt.xticks(ir,xr,fontsize=fc)
    plt.yticks(fontsize=fc)
    plt.axes().xaxis.set_minor_locator(MultipleLocator(0.05))
    plt.grid(axis="x",which="major")
    plt.ylabel(r'$E_f\ (eV/atom)$',fontsize=fc)
#    plt.xlabel(f'x{ele[1]}',fontsize=fc)
    plt.xlabel(r'$x_{'+ele[1]+'}$',fontsize=fc)
    plt.tight_layout()
    plt.savefig(string+".png")
    plt.show()
##---------------------------------------------------------------------------##

##---------count_n_plot(list)------------------------------------------------##
def count_n_plot(mylist):
    from collections import Counter
    import numpy as np
    import matplotlib.pyplot as plt
    labels, values = zip(*Counter(mylist).items())
    for kk in range(len(labels)):
        print(labels[kk], values[kk])
    index = np.arange(len(labels))
    width = .6
    fs=16
    plt.bar(index, values, width, color='b')
    plt.xticks(index, labels, fontsize=fs)
    plt.yticks(fontsize=fs)
    plt.ylabel("Count",fontsize=fs+2)
    plt.tight_layout()
    plt.savefig("count_n_plot.png")
    plt.show()
    return 0
##---------------------------------------------------------------------------##

##-------------cry_cluster_cut-----------------------------------------------##
# obtain clusters from crystal only based on a radius cutoff
def cry_cluster_cut(file_in, name_id, center, neighbor,rrcut):
    import numpy as np
    import os
    if not os.path.exists("clusters"):
        os.makedirs("clusters")
    fin=open(file_in,'r')
    rcut=rrcut # angstrom
    fin.readline()
    factor=float(fin.readline().split()[0])
    v=[[0 for i in range(3)] for j in range(3)]
    elong=[0 for i in range(3)]
    for i in range(3):
        line=fin.readline().split()
        vl=0
        for j in range(3): 
            v[i][j]=float(line[j])*factor
    vol=np.dot(v[0],np.cross(v[1],v[2]))
    b1=np.cross(v[1],v[2])
    b2=np.cross(v[2],v[1])
    b3=np.cross(v[0],v[1])
    recv=[[0 for i in range(3)] for j in range(3)]
    for i in range(3):
        recv[0][i] = b1[i]/vol
        recv[1][i] = b2[i]/vol
        recv[2][i] = b3[i]/vol
    rvsq1= np.dot(recv[0],recv[0])
    rvsq2= np.dot(recv[1],recv[1])
    rvsq3= np.dot(recv[2],recv[2])

    elong[0] = int(rcut*np.sqrt(rvsq1))+1
    elong[1] = int(rcut*np.sqrt(rvsq2))+1
    elong[2] = int(rcut*np.sqrt(rvsq3))+1


    atype=fin.readline().split()
    ntype=len(atype)
    ll=fin.readline().split()
    na=[]
    for i in ll:
        na.append(int(i))
    tpe=[]
    for i in range(ntype):
        for k in range(na[i]):
            tpe.append( atype[i] )
    
    natom=sum(na)
    ntotal=natom*(2*elong[0]+1)*(2*elong[1]+1)*(2*elong[2]+1)
    if(fin.readline().split()[0] != "Direct"):
        print("Fractional coordinate!")
        exit()

    x=[]
    for i in range(ntype):
        for k in range(na[i]):
            ll=fin.readline().split()
            x.append([float(ll[0]),float(ll[1]), float(ll[2])])
    fin.close()

    # make elongations
    xall=[]
    all_type=[]
    for k0 in range(int(-1)*elong[0],elong[0]+1):
        for k1 in range(int(-1)*elong[1],elong[1]+1):
            for k2 in range(int(-1)*elong[2],elong[2]+1):
                for ia in range(natom):
                    xall.append( [ x[ia][0] + k0, x[ia][1] + k1, x[ia][2] + k2 ] )
                    all_type.append(tpe[ia])
    if(len(xall) != ntotal): print("Warning!")
   
    print("elongated by" ,*elong,end=" ... " )
    # takecluster
    fout=open("clusters/cluster-total.xyz","w+")
    count=0
    for i in range(natom):
        if( tpe[i] != center ):
           continue
        count+=1
        
        # count # of atoms based only on the cutoff
        ic=0
        nei_xyz=[]
        nei_dis=[]
        nei_tpe=[]
        for j in range(ntotal):
            dd=[]
            for k in range(3):
                dd.append( (xall[j][0]-x[i][0])*v[0][k] + (xall[j][1]-x[i][1])*v[1][k] + (xall[j][2]-x[i][2])*v[2][k] )
            dis= np.sqrt( dd[0]**2 + dd[1]**2 + dd[2]**2)
            if( dis < rcut and all_type[j] in neighbor):
                ic += 1
                nei_xyz.append(dd)
                nei_dis.append(dis)
                nei_tpe.append(all_type[j])


        # pai xu
        for ki in range(ic-1):
            for kj in range(ki+1,ic):
                if(nei_dis[ki] > nei_dis[kj]):
                    exchange=nei_xyz[ki]
                    nei_xyz[ki]=nei_xyz[kj]
                    nei_xyz[kj]=exchange

                    exchange=nei_dis[ki]
                    nei_dis[ki]=nei_dis[kj]
                    nei_dis[kj]=exchange
                    
                    exchange=nei_tpe[ki]
                    nei_tpe[ki]=nei_tpe[kj]
                    nei_tpe[kj]=exchange
        
        if(ic==0):
            continue
        fsp=open("clusters/POSCAR_"+name_id+"-"+str(count)+"-"+nei_tpe[0]+".xyz","w+")
        print("%20d"%(ic+1), file=fout)
        print("%20s%10d"%(file_in,count), file=fout)
        print("%20d"%(ic+1), file=fsp)
        print("%20s%10d"%(file_in, count), file=fsp)
        print("%2s%12.6f%12.6f%12.6f"%(tpe[i], 0.0, 0.0, 0.0), file=fout)
        print("%2s%12.6f%12.6f%12.6f"%(tpe[i], 0.0, 0.0, 0.0), file=fsp)
        #j=0; nb_ct=0
        for j in range(ic):
            if(nei_tpe[j] in neighbor):
                print("%2s%12.6f%12.6f%12.6f%12.6f"%(nei_tpe[j], nei_xyz[j][0], nei_xyz[j][1], nei_xyz[j][2], nei_dis[j]), file=fout)
                print("%2s%12.6f%12.6f%12.6f%12.6f"%(nei_tpe[j], nei_xyz[j][0], nei_xyz[j][1], nei_xyz[j][2], nei_dis[j]), file=fsp)
        fsp.close()
    fout.close()
    return count
##---------------------------------------------------------------------------##

##---------------------RAS---------------------------------------------------##
def RAS(cluster): #cluster[0] must be the center atom
    import numpy as np
    import pylab as plt
    vc=np.array(cluster)
    nc=len(cluster)
    ang=[]
    for i in range(1,nc-1):
        v1=vc[i]-vc[0]
        abs_v1=np.sqrt(np.dot(v1,v1))
        for j in range(i+1,nc):
            v2=vc[j]-vc[0]
            abs_v2=np.sqrt(np.dot(v2,v2))
            c_theta=np.dot(v1,v2)/abs_v1/abs_v2
            if(c_theta>1.0): c_theta=1.0
            if(c_theta<-1.0): c_theta=-1.0
            theta = np.arccos(c_theta)/np.pi*180
            ang.append(theta)
    ave = np.average(ang)
    ang_max=[]; ang_min=[]
    for ia in ang:
        if(ia>ave):
            ang_max.append(ia)
        else:
            ang_min.append(ia)
    fout=open("ang-hist.dat","w+")
    bins = np.arange(0, 185, 2)
    hy,hx,patches=plt.hist(ang,bins=bins,density=True)
    for k in range(len(hy)):
        print(hx[k],hy[k],file=fout)
    fout.close()
    return np.average(ang_min), np.average(ang_max)
##---------------------------------------------------------------------------##
