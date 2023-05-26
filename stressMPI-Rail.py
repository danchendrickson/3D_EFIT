import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math
import time
import functools
from distBox import distBox
import sys
from mpi4py import MPI
from os import environ 
from typing import *                                                                       
#MPIComm = Union[MPI.Intracomm, MPI.Intercomm]
mpi_comm = MPI.COMM_WORLD
myid = mpi_comm.Get_rank()                                                         
mpi_size = mpi_comm.Get_size()        
nprocs=mpi_size

# for overlapping slabs:  
# # points per proc along z = npz = gh1/nproc (+2 to ghost boundaries)
# glob_index = loc_index-1 + npz*myid
# loc_index = glob_index - npz*myid + 1
# myid given glob_index = glob_index/npz = ghloc-2

# set Constants
#MATERIAL 1 ((steel))
pRatio1 = 0.29                                    #poission's ratio in 
yModulus1 = 200 * (10**9)                           #youngs modulus in pascals
rho1 = 7800                                        #density in kg/m^3

mu1 = yModulus1/(2*(1+pRatio1))                    #second Lame Parameter
lmbda1 = 2 * mu1 * pRatio1 / (1 - 2 * pRatio1)     #first Lame Parameter

#Calculate speed of longitudinal and transverse waves in material 1
cl1 = np.sqrt((lmbda1 + 2* mu1)/rho1)
ct1 = np.sqrt(mu1/rho1)

if myid == 0:
    print('MPI info: \n')
    print('MPI rank, size',myid,mpi_size)


if myid == 0: 
    print('material 1 wave speeds:' ,cl1,ct1)

#Choose ferquency to be used for excitment
frequency = 18000

#calculate wave lengths for material 1
omegaL1 = cl1 / frequency
omegaT1 = ct1 / frequency


#MATERIAL 2  (air)
pRatio2= 0.99
yModulus2= 1.13*(10**5)
rho2 = 1.15       
mu2 = yModulus2/(2*(1+pRatio2))                    
lmbda2 = abs(2 * mu2 * pRatio2 / (1 - 2 * pRatio2))

#Calculate speed of longitudinal and transverse waves in material 1
cl2= np.sqrt((lmbda2 + 2* mu2)/rho2)
ct2 = np.sqrt(mu2/rho2)

if myid == 0:
    print('material 2 wave speeds:' ,cl2,ct2)

#calculate wavelengths in material 2
omegaL2 = cl2 / frequency
omegaT2 = ct2 / frequency

#dimensions of materials in meters
#the dimensions of material 1 should be greater than material 2
#diensions of material 1 

length1 = 1.0
width1 = 0.1524
height1 = 0.1524

#is the rail supported by 0, 1 or 2 ties
Ties = 2

#Run for 4 Cycles:
runtime = 5 / frequency 

#Set time step and grid step to be 10 steps per frequency and ten steps per wavelength respectively
#ts = 1 / frequency / 10    #time step
gs = (min(omegaL1, omegaT1,omegaL2,omegaT2) /12)    #grid step
ts = gs/((max(cl1,ct1,cl2,ct2))*(np.sqrt(3)))*0.95 #time step

Tsteps = int(math.ceil(runtime / ts)) + 1       #total Time Steps

#number of grid points
gl1 = int(math.ceil(length1 / gs)) +1       #length 
gw1 = int(math.ceil(width1 / gs)) +1       #width
gh1 = int(math.ceil(height1 / gs)) +1       #height

frequency = 16355

# Keep these as the global values
xmax=gl1-1
ymax=gw1-1
zmax=gh1-1


#MPI EJW Section 1
#extend the length of the beam so that the number of nodes in the x dimmension 
#is the evenly divisible by the number of processors
if (gl1 % nprocs) != 0:
    gl1 += nprocs - (gl1 % nprocs)

#check you did it right
if (gl1 % nprocs) != 0:
    if myid == 0:
        print("Hey, gl1 not divisible by nproc",gl1,nprocs)
        sys.exit()
npx=int(gl1/nprocs)


if myid == 0:
    print("gl1,npx,nproc",gl1,npx,nprocs)

#print(runtime, ts, gs, Tsteps, gl, gh)

if myid == 0:
    print('runtime (s), time step size (s), total # of time steps:', runtime, ts, Tsteps)
    print('grid step size, # of length pts, # of height pts, # of width pts, gl1 loc pts:', gs,gl1,gw1,gh1,npx)

#tensor to store material properties for each point
#0 index is density
#1 index is first Lame Parmaeter
#2 index is second lame parameter

#MPI EJW Section 2 changes
matPropsglob=np.zeros((4,gl1,gw1,gh1))
matPropsglob[0,:,:,:]=rho1
matPropsglob[1,:,:,:]=lmbda1
matPropsglob[2,:,:,:]=mu1
matPropsglob[3,:,:,:]=0

#Make top surface work hardened:
whlayer = int(0.0001/gs)

matPropsglob[0,:,:,gh1-whlayer:gh1-1]=rho1*1.25
matPropsglob[1,:,:,gh1-whlayer:gh1-1]=lmbda1*1.5

if myid == 0:
    print('globs made, line 145')

AirCut = True
if AirCut:
    #zone 1 of air, left of web
    for yy in range(int(13/36*gw1)):
        y = yy + 0
        for zz in range(int(16/36*gh1)):
            z = zz + int(8/36*gh1)
            matPropsglob[0,:,y,z]=rho2
            matPropsglob[1,:,y,z]=lmbda2
            matPropsglob[2,:,y,z]=mu2
            matPropsglob[3,:,y,z]=99
            
    # zone 2 of air left of head
    for yy in range(int(4/36*gw1)):
        y = yy + 0
        for zz in range(int(12/36*gh1)):
            z = zz + int(24/36*gh1)
            matPropsglob[0,:,y,z]=rho2
            matPropsglob[1,:,y,z]=lmbda2
            matPropsglob[2,:,y,z]=mu2
            matPropsglob[3,:,y,z]=99
            
    # zone 3 of air, right of web
    for yy in range(int(13/36*gw1)):
        y = ymax - yy
        for zz in range(int(16/36*gh1)):
            z = zz + int(8/36*gh1)
            matPropsglob[0,:,y,z]=rho2
            matPropsglob[1,:,y,z]=lmbda2
            matPropsglob[2,:,y,z]=mu2
            matPropsglob[3,:,y,z]=99
            
    # zone 4 of air, right of head
    for yy in range(int(4/36*gw1)):
        y = ymax - yy
        for zz in range(int(12/36*gh1)):
            z = zz + int(24/36*gh1)
            matPropsglob[0,:,y,z]=rho2
            matPropsglob[1,:,y,z]=lmbda2
            matPropsglob[2,:,y,z]=mu2
            matPropsglob[3,:,y,z]=99

#set the boundary conditions in material props4
Ties = 2

# top
matPropsglob[3,:,:,zmax]=2

# bottom 
matPropsglob[3,:,:,0]=1

# left
matPropsglob[3,:,0,:]=3

# right
matPropsglob[3,:,ymax,:]=4

# bottom right
matPropsglob[3,:,0,zmax]=13

# top rigight
matPropsglob[3,:,zmax,zmax]=14

#top of footing
z = int(8/36*gh1)
for yy in range(int(13/36*gw1)):
    y = yy + 0
    matPropsglob[3,:,y,z]=2
for yy in range(int(13/36*gw1)):
    y = yy + int(23/36*gw1)
    matPropsglob[3,:,y,z]=2

#sides of web
yl=int(13/36*gw1)
yr=int(23/36*gw1)
for zz in range(int(16/36*gh1)):
    z = zz + int(8/36*gh1)
    matPropsglob[3,:,yl,z]=3
    matPropsglob[3,:,yr,z]=4


#bottom of head
z=int(24/36*gh1)
for yy in range(int(8/36*gw1)):
    y=yy+int(5/36*gw1)
    matPropsglob[3,:,y,z]=1
    y=yy +int(23/36*gw1)
    matPropsglob[3,:,y,z]=1

#sides of head
yl=int(4/36*gw1)
yr=int(32/36*gw1)
for zz in range(int(12/36*gh1)):
    z = zmax - zz
    matPropsglob[3,:,yl,z]=3
    matPropsglob[3,:,yr,z]=4
    
#sides of foot
for zz in range(int(7/36*gh1)):
    z = zz
    matPropsglob[3,:,0,z]=3
    matPropsglob[3,:,ymax,z]=4


#top edge of foot on left
matPropsglob[3,:,0,int(8/36*gh1)]=13
#top endge of foot on right
matPropsglob[3,:,ymax,int(8/36*gh1)]=14

#bottom of head on left
matPropsglob[3,:,int(4/36*gw1),int(24/36*gh1)]=9
#bottom of head on right
matPropsglob[3,:,int(32/36*gw1),int(24/36*gh1)]=10

#Top of head on left
matPropsglob[3,:,int(4/36*gw1),zmax]=13
#top of head on right
matPropsglob[3,:,int(32/36*gw1),zmax]=14


if Ties ==2:  #tie on both end, absorbing all vertical velocity in square at end of track
    #face
    matPropsglob[3,0:gw1,:,0]=35
    matPropsglob[3,xmax-gw1:xmax,:,0]=35

    #Edges
    matPropsglob[3,0,:,0]=27
    matPropsglob[3,xmax,:,0]=28

    matPropsglob[3,0:gw1,0,0]=29
    matPropsglob[3,xmax-gw1:xmax,0,0]=29
    matPropsglob[3,0:gw1,ymax,0]=30
    matPropsglob[3,xmax-gw1:xmax,:,0]=30

    #corners
    matPropsglob[3,0,0,0] = 31
    matPropsglob[3,0,ymax,0]= 32
    matPropsglob[3,xmax,0,0]=33
    matPropsglob[3,xmax,ymax,0]=34

elif Ties == 1:  #tie in the middle
    half=int(gl1/2)
    halfwidth = int(gw1/2)
    start =half - halfwidth
    end = half+halfwidth

    #face
    matPropsglob[3,start:end,:,0]=35

    #edge
    matPropsglob[3,start:end,0,0]=29
    matPropsglob[3,start:end,ymax,0]=30

elif Ties == 3: #time on both end and in middle
    #end ties
    #face
    matPropsglob[3,0:gw1,:,0]=35
    matPropsglob[3,xmax-gw1:xmax,:,0]=35

    #Edges
    matPropsglob[3,0,:,0]=27
    matPropsglob[3,xmax,:,0]=28

    matPropsglob[3,0:gw1,0,0]=29
    matPropsglob[3,xmax-gw1:xmax,0,0]=29
    matPropsglob[3,0:gw1,ymax,0]=30
    matPropsglob[3,xmax-gw1:xmax,:,0]=30

    #corners
    matPropsglob[3,0,0,0] = 31
    matPropsglob[3,0,ymax,0]= 32
    matPropsglob[3,xmax,0,0]=33
    matPropsglob[3,xmax,ymax,0]=34   

    #middle Tie
    half=int(gl1/2)
    halfwidth = int(gw1/2)
    start =half - halfwidth
    end = half+halfwidth

    #face
    matPropsglob[3,start:end,:,0]=35

    #edge
    matPropsglob[3,start:end,0,0]=29
    matPropsglob[3,start:end,ymax,0]=30

if myid == 0:
    print('air cuts made, line 310')

#define sine-exponential wave excitation

timeVec=np.linspace(0,runtime,Tsteps)

#MPI EJW Section #3 changes
#radius
r=3
inputx=2
inputy=int(gw1/2)
inputz=int(gh1/2)

# get loc by formula

inputid=int(inputx / npx)
inputlocx=int(inputx - inputid*npx+1)

if (myid == 0) :
    print("line 369: glb inputx, local inputx id, local inputx:  ",inputx,inputid,inputlocx)


szzConst=2*ts/(gs*rho1)

amp=10000
decayRate= 0
sinConst=ts*amp/rho1

sinInputSignal=sinConst*np.sin(2*np.pi*frequency*timeVec)*np.exp(-decayRate*timeVec)

# MPI EJW Section #4 changes 

#initialize fields
vx=np.zeros((npx+2,gw1,gh1))
vy=np.zeros((npx+2,gw1,gh1))
vz=np.zeros((npx+2,gw1,gh1))

sxx=np.zeros((npx+2,gw1,gh1))
syy=np.zeros((npx+2,gw1,gh1))
szz=np.zeros((npx+2,gw1,gh1))
sxy=np.zeros((npx+2,gw1,gh1))
sxz=np.zeros((npx+2,gw1,gh1))
syz=np.zeros((npx+2,gw1,gh1))

#record the signal at a specified location
### ADD map function for this
signalLocx=int(gl1/2)
signalLocy=int(gw1/2)
signalLocz=int(gh1/2)

#SAME AS INPUTZ?
signalLocxid=int(signalLocx / npx)
signalLocxlocx=int(signalLocx - myid*npx+1)

vxSignal=np.zeros(Tsteps)
vySignal=np.zeros(Tsteps)
vzSignal=np.zeros(Tsteps)

# Grab splits and offsets for scattering arrays
# Only thing to scatter is matPropsglob
# v's and s's are zero to start + source applied later 
# in single proc's array
if myid == 0:
    split=np.zeros(nprocs)
    split[:]=gw1*gh1*npx

    offset=np.zeros(nprocs)
    for i in range(nprocs):
        offset[i]=i*gw1*gh1*npx
else:
    split=None
    offset=None

split=mpi_comm.bcast(split)
offset=mpi_comm.bcast(offset)

matProps0 = np.zeros((npx,gw1,gh1))
matProps1 = np.zeros((npx,gw1,gh1))
matProps2 = np.zeros((npx,gw1,gh1))
matProps3 = np.zeros((npx,gw1,gh1))

mpi_comm.Scatterv([matPropsglob[0,:,:,:],split,offset,MPI.DOUBLE], matProps0)
mpi_comm.Scatterv([matPropsglob[1,:,:,:],split,offset,MPI.DOUBLE], matProps1)
mpi_comm.Scatterv([matPropsglob[2,:,:,:],split,offset,MPI.DOUBLE], matProps2)
mpi_comm.Scatterv([matPropsglob[3,:,:,:],split,offset,MPI.DOUBLE], matProps3)

matProps0=distBox(matProps0,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
matProps1=distBox(matProps1,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
matProps2=distBox(matProps2,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
matProps3=distBox(matProps3,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        

#Now slab has local versions with ghosts of matProps
if (myid == 0) :
    print('split matprops to globs, scratttered parameters to processors, line 449')


def updateStress(x,y,z):
        
    #Calculate constants for stress equations
    norm1=(1/gs)*(matProps1[x,y,z]+2*matProps2[x,y,z])
    norm2=(1/gs)*(matProps1[x,y,z])


    try:
        shearDenomxy=(1/matProps2[x,y,z])+(1/matProps2[x+1,y,z])+(1/matProps2[x,y+1,z])+(1/matProps2[x+1,y+1,z])
        shearxy=4*(1/gs)*(1/shearDenomxy)
    except:
        pass
    
    try:
        shearDenomxz=(1/matProps2[x,y,z])+(1/matProps2[x+1,y,z])+(1/matProps2[x,y,z+1])+(1/matProps2[x+1,y,z+1])
        shearxz=4*(1/gs)*(1/shearDenomxz)
    except:
        pass
    
    try:
        shearDenomyz=(1/matProps2[x,y,z])+(1/matProps2[x,y+1,z])+(1/matProps2[x,y,z+1])+(1/matProps2[x,y+1,z+1])
        shearyz=4*(1/gs)*(1/shearDenomyz)
    except:
        pass
    try:
        #FACES
        if matProps3[x,y,z] == 0:
            norm1=(1/gs)*(matProps1[x,y,z]+2*matProps2[x,y,z])
            norm2=(1/gs)*(matProps1[x,y,z])

            shearDenomxy=(1/matProps2[x,y,z])+(1/matProps2[x+1,y,z])+(1/matProps2[x,y+1,z])+(1/matProps2[x+1,y+1,z])
            shearxy=4*(1/gs)*(1/shearDenomxy)

            shearDenomxz=(1/matProps2[x,y,z])+(1/matProps2[x+1,y,z])+(1/matProps2[x,y,z+1])+(1/matProps2[x+1,y,z+1])
            shearxz=4*(1/gs)*(1/shearDenomxz)

            shearDenomyz=(1/matProps2[x,y,z])+(1/matProps2[x,y+1,z])+(1/matProps2[x,y,z+1])+(1/matProps2[x,y+1,z+1])
            shearyz=4*(1/gs)*(1/shearDenomyz)

            ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
            sxx[x,y,z]=sxx[x,y,z]+ds*ts

            ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
            syy[x,y,z]=syy[x,y,z]+ds*ts

            ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
            szz[x,y,z]=szz[x,y,z]+ds*ts

            ds=shearxy*(vx[x,y+1,z]-vx[x,y,z]+vy[x+1,y,z]-vy[x,y,z])
            sxy[x,y,z]=sxy[x,y,z]+ds*ts

            ds=shearxz*(vx[x,y,z+1]-vx[x,y,z]+vz[x+1,y,z]-vz[x,y,z])
            sxz[x,y,z]=sxz[x,y,z]+ds*ts   

            ds=shearyz*(vy[x,y,z+1]-vy[x,y,z]+vz[x,y+1,z]-vz[x,y,z])
            syz[x,y,z]=syz[x,y,z]+ds*ts

        elif matProps3[x,y,z] == 1 or matProps3[x,y,z] == 35:

            ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z])
            sxx[x,y,z]=sxx[x,y,z]+ds*ts

            ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z])
            syy[x,y,z]=syy[x,y,z]+ds*ts

            szz[x,y,z]=-szz[x,y,z+1]

            ds=shearxy*(vx[x,y+1,z]-vx[x,y,z]+vy[x+1,y,z]-vy[x,y,z])
            sxy[x,y,z]=sxy[x,y,z]+ds*ts

            sxz[x,y,z]=0
            syz[x,y,z]=0

        elif matProps3[x,y,z] == 5:


            sxx[x,y,x]=-sxx[x+1,y,z]

            ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vz[x,y,z]-vz[x,y,z-1])
            syy[x,y,z]=syy[x,y,z]+ds*ts

            ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vy[x,y,z]-vy[x,y-1,z])
            szz[x,y,z]=szz[x,y,z]+ds*ts

            sxy[x,y,z]=0
            sxz[x,y,z]=0

            ds=shearyz*(vy[x,y,z+1]-vy[x,y,z]+vz[x,y+1,z]-vz[x,y,z])
            syz[x,y,z]=syz[x,y,z]+ds*ts

        elif matProps3[x,y,z] == 6:
            sxx[x,y,z]=-sxx[x-1,y,z]

            ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
            syy[x,y,z]=syy[x,y,z]+ds*ts

            ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
            szz[x,y,z]=szz[x,y,z]+ds*ts

            sxy[x,y,z]=0
            sxz[x,y,z]=0

            ds=shearyz*(vy[x,y,z+1]-vy[x,y,z]+vz[x,y+1,z]-vz[x,y,z])
            syz[x,y,z]=syz[x,y,z]+ds*ts

        elif matProps3[x,y,z] == 3:

            ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vz[x,y,z]-vz[x,y,z-1])
            sxx[x,y,z]=sxx[x,y,z]+ds*ts

            syy[x,y,z]=-syy[x,y+1,z]

            ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z])
            szz[x,y,z]=szz[x,y,z]+ds*ts

            sxy[x,y,z]=0

            ds=shearxz*(vx[x,y,z+1]-vx[x,y,z]+vz[x+1,y,z]-vz[x,y,z])
            sxz[x,y,z]=sxz[x,y,z]+ds*ts

            syz[x,y,z]=0

        elif matProps3[x,y,z] == 4:

            ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
            sxx[x,y,z]=sxx[x,y,z]+ds*ts

            syy[x,y,z]=-syy[x,y-1,z]

            ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
            szz[x,y,z]=szz[x,y,z]+ds*ts

            sxy[x,y,z]=0

            ds=shearxz*(vx[x,y,z+1]-vx[x,y,z]+vz[x+1,y,z]-vz[x,y,z])
            sxz[x,y,z]=sxz[x,y,z]+ds*ts

            syz[x,y,z]=0

        elif matProps3[x,y,z] == 2:

            ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
            sxx[x,y,z]=sxx[x,y,z]+ds*ts

            ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
            syy[x,y,z]=syy[x,y,z]+ds*ts

            szz[x,y,z]=-szz[x,y,z-1]

            ds=shearxy*(vx[x,y+1,z]-vx[x,y,z]+vy[x+1,y,z]-vy[x,y,z])
            sxy[x,y,z]=sxy[x,y,z]+ds*ts

            sxz[x,y,z]=0
            syz[x,y,z]=0


        #EDGES
        #bottom edges
        elif matProps3[x,y,z] == 7 or matProps3[x,y,z] == 27:

            sxx[x,y,z]=-sxx[x+1,y,z]

            ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
            syy[x,y,z]=syy[x,y,z]+ds*ts

            szz[x,y,z]=-szz[x,y,z+1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        elif matProps3[x,y,z] == 8 or matProps3[x,y,z] == 28:

            sxx[x,y,z]=-sxx[x-1,y,z]

            ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
            syy[x,y,z]=syy[x,y,z]+ds*ts

            szz[x,y,z]-szz[x,y,z+1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0


        elif matProps3[x,y,z] == 9 or matProps3[x,y,z] == 29:

            ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
            sxx[x,y,z]=sxx[x,y,z]+ds*ts

            syy[x,y,z]=-syy[x,y+1,z]

            szz[x,y,z]=-szz[x,y,z+1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0


        elif matProps3[x,y,z] == 10 or matProps3[x,y,z] == 30:

            ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
            sxx[x,y,z]=sxx[x,y,z]+ds*ts

            syy[x,y,z]=-syy[x,y-1,z]

            szz[x,y,z]=-szz[x,y,z+1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        #side edges
        elif matProps3[x,y,z] == 15:

            sxx[x,y,z]=-sxx[x+1,y,z]

            syy[x,y,z]=-syy[x,y+1,z]

            ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
            szz[x,y,z]=szz[x,y,z]+ds*ts

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0


        elif matProps3[x,y,z] == 17:

            sxx[x,y,z]=-sxx[x-1,y,z]

            syy[x,y,z]=-syy[x,y+1,z]

            ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
            szz[x,y,z]=szz[x,y,z]+ds*ts

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        elif matProps3[x,y,z] == 16:

            sxx[x,y,z]=-sxx[x+1,y,z]

            syy[x,y,z]=-syy[x,y-1,z]

            ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
            szz[x,y,z]=szz[x,y,z]+ds*ts

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0


        elif matProps3[x,y,z] == 18:

            sxx[x,y,z]=-sxx[x-1,y,z]

            syy[x,y,z]=-syy[x,y-1,z]

            ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
            szz[x,y,z]=szz[x,y,z]+ds*ts

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0


        #top edges
        elif matProps3[x,y,z] == 11:

            sxx[x,y,z]=-sxx[x+1,y,z]

            ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
            syy[x,y,z]=syy[x,y,z]+ds*ts

            szz[x,y,z]=-szz[x,y,z-1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        elif matProps3[x,y,z] == 12:

            sxx[x,y,z]=-sxx[x-1,y,z]

            ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
            syy[x,y,z]=syy[x,y,z]+ds*ts

            szz[x,y,z]=-szz[x,y,z-1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        elif matProps3[x,y,z] == 13:
            ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
            sxx[x,y,z]=sxx[x,y,z]+ds*ts

            syy[x,y,z]=-syy[x,y+1,z]

            szz[x,y,z]=-szz[x,y,z-1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0


        elif matProps3[x,y,z] == 14:
            ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
            sxx[x,y,z]=sxx[x,y,z]+ds*ts

            syy[x,y,z]=-syy[x,y-1,z]

            szz[x,y,z]=-szz[x,y,z-1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0


        #CORNERS

        elif matProps3[x,y,z] == 19 or matProps3[x,y,z] == 31:

            sxx[x,y,z]=-sxx[x+1,y,z]

            syy[x,y,z]=-syy[x,y+1,z]

            szz[x,y,z]=-szz[x,y,z+1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        elif matProps3[x,y,z] ==  20:

            sxx[x,y,z]=-sxx[x+1,y,z]

            syy[x,y,z]=-syy[x,y+1,z]

            szz[x,y,z]=-szz[x,y,z-1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        elif matProps3[x,y,z] ==  21 or matProps3[x,y,z] == 32:

            sxx[x,y,z]=-sxx[x+1,y,z]

            syy[x,y,z]=-syy[x,y-1,z]

            szz[x,y,z]=-szz[x,y,z+1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        elif matProps3[x,y,z] == 22:

            sxx[x,y,z]=-sxx[x+1,y,z]

            syy[x,y,z]=-syy[x,y-1,z]

            szz[x,y,z]=-szz[x,y,z-1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        elif matProps3[x,y,z] == 23 or matProps3[x,y,z] == 33:

            sxx[x,y,z]=-sxx[x-1,y,z]

            syy[x,y,z]=-syy[x,y+1,z]

            szz[x,y,z]=-szz[x,y,z+1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        elif matProps3[x,y,z] == 24:

            sxx[x,y,z]=-sxx[x-1,y,z]

            syy[x,y,z]=-syy[x,y+1,z]

            szz[x,y,z]=-szz[x,y,z-1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        elif matProps3[x,y,z] == 25 or matProps3[x,y,z] == 34:

            sxx[x,y,z]=-sxx[x-1,y,z]

            syy[x,y,z]=-syy[x,y-1,z]

            szz[x,y,z]=-szz[x,y,z+1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        elif matProps3[x,y,z] ==  26:

            sxx[x,y,z]=-sxx[x-1,y,z]

            syy[x,y,z]=-syy[x,y-1,z]

            szz[x,y,z]=-szz[x,y,z-1]

            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0


        elif matProps3[x,y,z] == 99:
            sxx[x,y,z]=0
            syy[x,y,z]=0
            szz[x,y,z]=0
            sxy[x,y,z]=0
            sxz[x,y,z]=0
            syz[x,y,z]=0

        else: print('error:', str(x), str(y), str(z))
    except:
        print('Boundary Conditon isssue Stress: ', str(x), str(y), str(z), str(matProps3[x,y,z]))

# %%
def updateVelocity(x,y,z):
    try:
        if x!=xmax:
            dvxConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x+1,y,z]))
            vminx=(2*ts)/(matProps0[x+1,y,z]*gs)

        if y!=ymax:
            dvyConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y+1,z]))
            vminy=(2*ts)/(matProps0[x,y+1,z]*gs)

        if z!=zmax:
            dvzConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y,z+1]))
            vminz=(2*ts)/(matProps0[x,y,z+1]*gs)

        vmax=(2*ts)/(matProps0[x,y,z]*gs)

        #FACES
        if matProps3[x,y,z] == 0:
            dvxConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x+1,y,z]))
            dvyConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y+1,z]))
            dvzConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y,z+1]))

            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts

            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts

            dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
            vz[x,y,z]=vz[x,y,z]+dv*ts

        elif matProps3[x,y,z] == 1:
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts

            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]

        elif matProps3[x,y,z] == 5:

            vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]

            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts

            dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
            vz[x,y,z]=vz[x,y,z]+dv*ts

        elif matProps3[x,y,z] == 6:

            vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]

            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts

            dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
            vz[x,y,z]=vz[x,y,z]+dv*ts

        elif matProps3[x,y,z] == 3:

            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts

            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]

            dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
            vz[x,y,z]=vz[x,y,z]+dv*ts


        elif matProps3[x,y,z] == 4:
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts

            vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]

            dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
            vz[x,y,z]=vz[x,y,z]+dv*ts                     


        elif matProps3[x,y,z] == 2:
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts

            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]

        #EDGES
        #bottom edges
        elif matProps3[x,y,z] == 7:

            vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]

            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]

        elif matProps3[x,y,z] == 8:

            vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]

            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]

        elif matProps3[x,y,z] == 9:

            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts

            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]

            vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]

        elif matProps3[x,y,z] == 10:

            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts

            vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]

            vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]

        #side edges
        elif matProps3[x,y,z] == 15:

            vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]

            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]

            dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
            vz[x,y,z]=vz[x,y,z]+dv*ts

        elif matProps3[x,y,z] == 17:

            vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]

            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]

            dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
            vz[x,y,z]=vz[x,y,z]+dv*ts


        elif matProps3[x,y,z] == 16:

            vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]

            vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]

            dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
            vz[x,y,z]=vz[x,y,z]+dv*ts



        elif matProps3[x,y,z] == 18:

            vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]

            vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]

            dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
            vz[x,y,z]=vz[x,y,z]+dv*ts


        #top edges
        elif matProps3[x,y,z] == 11:

            vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]

            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]


        elif matProps3[x,y,z] == 12:

            vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]

            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]

        elif matProps3[x,y,z] == 13:

            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts

            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]

            vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]


        elif matProps3[x,y,z] == 14:

            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts

            vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]

            vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]


        #CORNERS
        elif matProps3[x,y,z] == 19:

            vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
            vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]

        elif matProps3[x,y,z] == 20:

            vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
            vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]

        elif matProps3[x,y,z] == 21:

            vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
            vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
            vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]

        elif matProps3[x,y,z] == 22:

            vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
            vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
            vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]

        elif matProps3[x,y,z] == 23:

            vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
            vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]

        elif matProps3[x,y,z] == 24:

            vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
            vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]

        elif matProps3[x,y,z] == 25:

            vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
            vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
            vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]

        elif matProps3[x,y,z] == 26:

            vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
            vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
            vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]

        #veleocity blocking Boundaries
        #face
        elif matProps3[x,y,z] == 35:
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts

            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z]=0
        #edges
        elif matProps3[x,y,z] == 27: #7
            vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts
            vz[x,y,z]=0

        elif matProps3[x,y,z] == 28: #8
            vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]

            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts
            vz[x,y,z]=0

        elif matProps3[x,y,z] == 29: #9
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts
            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
            vz[x,y,z]=0

        elif matProps3[x,y,z] == 30: #10
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts
            vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
            vz[x,y,z]=0

        #corners
        elif matProps3[x,y,z] == 31: #19
            vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
            vz[x,y,z]=0
        elif matProps3[x,y,z] == 32: #21
            vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
            vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
            vz[x,y,z]=0
        elif matProps3[x,y,z] == 33: #23
            vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
            vz[x,y,z]=0
        elif matProps3[x,y,z] == 34: #25
            vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
            vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
            vz[x,y,z]=0


        #incase non boundaries or unknown areas snuck through
        elif matProps3[x,y,z] == 99:
            vx[x,y,z]=0
            vy[x,y,z]=0
            vz[x,y,z]=0

        else: print('error: ',x,y,z, matProps3[x,y,z])
    except:
        print('Boundary Conditon isssue Velocity: ', str(x), str(y), str(z), str(matProps3[x,y,z]))
stime = time.time()

if (myid == 0 ):
    print('subs setup, line 1213.  About to start at ' + str(stime))
    

for t in range(0,Tsteps):
#for t in range(0,40):
    
    #sin-exponential input
    if (myid==inputid) :
        vz[inputlocx,inputy,inputz]=vz[inputlocx,inputy,inputz]-szzConst * szz[inputlocx,inputy,inputz] + sinInputSignal[t]

#    for x in range(gl1):
    for x in range(1,npx+1):
        for y in range(gw1):
            for z in range(gh1):
                updateStress(x,y,z)

    # cut boundaries off of arrays
    sxxt=sxx[1:npx+1,:,:]
    syyt=syy[1:npx+1,:,:]
    szzt=szz[1:npx+1,:,:]
    sxyt=sxy[1:npx+1,:,:]
    sxzt=sxz[1:npx+1,:,:]
    syzt=syz[1:npx+1,:,:]

    # redistrubute ghost/boundary values
    sxx=distBox(sxxt,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
    syy=distBox(syyt,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
    szz=distBox(szzt,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
    sxy=distBox(sxyt,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
    sxz=distBox(sxzt,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
    syz=distBox(syzt,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        

    for x in range(1,npx+1):
        for y in range(gw1):
            for z in range(gh1):
                updateVelocity(x,y,z)

    # cut boundaries off of arrays
    vxt=vx[1:npx+1,:,:]
    vyt=vy[1:npx+1,:,:]
    vzt=vz[1:npx+1,:,:]

    # redistrubute ghost/boundary values
    vx=distBox(vxt,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
    vy=distBox(vyt,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
    vz=distBox(vzt,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        

    #record signals
    if (myid==signalLocxid) :
        vxSignal[t]=vx[signalLocxlocx,signalLocy,signalLocz]
        vySignal[t]=vy[signalLocxlocx,signalLocy,signalLocz]
        vzSignal[t]=vz[signalLocxlocx,signalLocy,signalLocz]

    # save vx cut figure
    # ADD GATHER for plotting

    vzg = np.zeros((gl1,gw1,gh1))
    if t%10==0:
        vzt=vz[1:npx+1,:,:]        
        mpi_comm.Gatherv(vzt,[vzg,split,offset,MPI.DOUBLE])
        if (myid == 0 ) :
            fig=plt.figure()
            plt.contourf(np.transpose(vzg[:,:,int(gh1/2)]), cmap='seismic')
            plt.savefig('vzmidCut'+str(t)+'.png')
            plt.close(fig)        
    

    # Collect vx, sxx checksum contributions for printing
    vxt=vx[1:npx+1,:,:]
    sxxt=sxx[1:npx+1,:,:]

    ckvs=np.array(0.0,'d')
    ckss=np.array(0.0,'d')
    
    ckv=np.sum(np.absolute(vxt))
    cks=np.sum(np.absolute(sxxt))
    mpi_comm.Reduce(ckv,ckvs,op=MPI.SUM,root=0)
    mpi_comm.Reduce(cks,ckss,op=MPI.SUM,root=0)

    if (myid == 0 ):
        print(t,'/',Tsteps-1,'checksums vx, sxx:',ckvs,ckss, time.time()-stime)
    sys.stdout.flush()

if (myid == signalLocxid) :
    plt.clf()
    plt.plot(vxSignal)
    plt.savefig('vxsignal.png')

if (myid == signalLocxid) :
    plt.clf()
    plt.plot(vySignal)
    plt.savefig('vysignal.png')

if (myid == signalLocxid) :
    plt.clf()
    plt.plot(vzSignal)
    plt.savefig('vzsignal.png')


    
