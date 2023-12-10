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
AirCut = False
RailShape = False

#Dimmesnsion of simulation space in meters
length1 = 1
width1 = 0.2
height1 = 0.2

#Image Folder
imFolder = '/sciclone/scr10/dchendrickson01/EFIT/'
runName = 'JBStyleRodSupDenseAcu'

#is the rail supported by 0, 1 or 2 ties
Ties = 0

#Choose ferquency to be used for excitment
frequency = 640000
#frequency = 8000

figDPI = 600



#Forcing Function Location and type
# 1 for dropped wheel on top
# 2 for rubbing flange on side
# 3 for plane wave 

FFunction = 3

WheelLoad = 173000 #crane force in Neutons

#MATERIAL 1 ((steel))
pRatio1 = 0.29                                    #poission's ratio in 
yModulus1 = 200 * (10**9)                           #youngs modulus in pascals
rho1 = 7800                                        #density in kg/m^3


#CALCULATED PARAMETERS FROM INPUTS

mu1 = yModulus1/(2*(1+pRatio1))                    #second Lame Parameter
lmbda1 = 2 * mu1 * pRatio1 / (1 - 2 * pRatio1)     #first Lame Parameter
#Calculate speed of longitudinal and transverse waves in material 1
cl1 = np.sqrt((lmbda1 + 2* mu1)/rho1)
ct1 = np.sqrt(mu1/rho1)

#calculate wave lengths for material 1
omegaL1 = cl1 / frequency
omegaT1 = ct1 / frequency

#Image Folder
if FFunction == 1:
    imFolder += 'TopHit/'
elif FFunction == 2:
    imFolder += 'SideRub/'
elif FFunction == 3:
    imFolder += 'BiggerAcTii/'
elif FFunction == 4:
    imFolder += 'BiggerAcTii/'

'''
#MATERIAL 2  (made up)
pRatio2= 0.3
yModulus2= 100*(10**8)
rho2 = 3000       
mu2 = yModulus2/(2*(1+pRatio2))                    
lmbda2 = abs(2 * mu2 * pRatio2 / (1 - 2 * pRatio2))

#Calculate speed of longitudinal and transverse waves in material 1
cl2= np.sqrt((lmbda2 + 2* mu2)/rho2)
ct2 = np.sqrt(mu2/rho2)

#calculate wavelengths in material 2
omegaL2 = cl2 / frequency
omegaT2 = ct2 / frequency

if myid == 0:
    print('material 2 wave speeds:' ,cl2,ct2)
'''

#Set time step and grid step to be 10 steps per frequency and ten steps per wavelength respectively
#ts = 1 / frequency / 10    #time step
gs = (min(omegaL1, omegaT1) /12)    #grid step, omegaL2,omegaT2
ts = gs/((max(cl1,ct1))*(np.sqrt(3)))*0.95 #time step, cl2,ct2


#change to lower frequency but with dense grid
frequency = 16300

#Run for 4 Cycles:
runtime = 20 / frequency 

Tsteps = int(math.ceil(runtime / ts)) + 1       #total Time Steps

#number of grid points
gl1 = int(math.ceil(length1 / gs)) +1       #length 
gw1 = int(math.ceil(width1 / gs)) +1       #width
gh1 = int(math.ceil(height1 / gs)) +1       #height

print('runtime, gs, ts, gl, gw, gh, Tsteps, imFolder : ', runtime, gs, ts, gl1, gw1, gh1, Tsteps, imFolder)

# Keep these as the global values
xmax=gl1-1
ymax=gw1-1
zmax=gh1-1

#####




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
signalLocation=np.zeros((gl1,gw1,gh1))

matPropsglob[0,:,:,:]=rho1
matPropsglob[1,:,:,:]=lmbda1
matPropsglob[2,:,:,:]=mu1
matPropsglob[3,:,:,:]=0

#un comment for wedge test
'''
for x in range(gl1):
    for y in range(gw1):
        matPropsglob[0,x,y,:]=rho2
        matPropsglob[1,x,y,:]=lmbda2
        matPropsglob[2,x,y,:]=mu2
'''


#Make the Signal Location grid
if FFunction == 1:
    pnodes = max(int(whlayer / 2),3)
    contactLength = max(int(0.001 / gs),3)  #1 cm contact patch or 3 nodes, whichever is larger
    
    #starting at .25 down, to be between the first 2 ties
    WheelStartPoint = int(0.25 * gl1)
    
    signalLocation[WheelStartPoint:WheelStartPoint+contactLength,gridStartHeadWidth:gridEndHeadWidth, -3:] = 1
    
elif FFunction == 2:
    pnodes = int(whlayer / 4)
    contactLength = max(int(0.004 / gs),3)  #4 cm contact patch

    signalLocation[0:contactLength,gridStartHeadWidth:gridStartHeadWidth+3, gridStartHead:] = 1
    if myid == 0:
        print(FFunction, contactLength, np.sum(signalLocation))

    ## Find the share of the force per node for FF1

elif FFunction == 3:
    signalLocation[2:4,:,:] = 1


elif FFunction == 4:
    start = 2*int(gh1/5)
    end = 3 * int(gh1/5)
    signalLocation[2:4,start:end,start:end] = 1


if myid == 0:
    print('globs made, line 145')

#########
# FUnctions
def JBSU(x,y,z):
    if matProps3[x,y,z] in [0] : #,2,4,6]:
        norm1=(1/gs)*(matProps1[x,y,z]+2*matProps2[x,y,z])
        norm2=(1/gs)*(matProps1[x,y,z])

        ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
        sxx[x,y,z]=sxx[x,y,z]+ds*ts

        ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
        syy[x,y,z]=syy[x,y,z]+ds*ts

        ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
        szz[x,y,z]=szz[x,y,z]+ds*ts
    
    if matProps3[x,y,z] in [0,5]:
        shearDenomxy=(1/matProps2[x,y,z])+(1/matProps2[x+1,y,z])+(1/matProps2[x,y+1,z])+(1/matProps2[x+1,y+1,z])
        shearxy=4*(1/gs)*(1/shearDenomxy)
        ds=shearxy*(vx[x,y+1,z]-vx[x,y,z]+vy[x+1,y,z]-vy[x,y,z])
        sxy[x,y,z]=sxy[x,y,z]+ds*ts

    if matProps3[x,y,z] in [0,3]:
        shearDenomxz=(1/matProps2[x,y,z])+(1/matProps2[x+1,y,z])+(1/matProps2[x,y,z+1])+(1/matProps2[x+1,y,z+1])
        shearxz=4*(1/gs)*(1/shearDenomxz)
        ds=shearxz*(vx[x,y,z+1]-vx[x,y,z]+vz[x+1,y,z]-vz[x,y,z])
        sxz[x,y,z]=sxz[x,y,z]+ds*ts   

    if matProps3[x,y,z] in [0,1]:
        shearDenomyz=(1/matProps2[x,y,z])+(1/matProps2[x,y+1,z])+(1/matProps2[x,y,z+1])+(1/matProps2[x,y+1,z+1])
        shearyz=4*(1/gs)*(1/shearDenomyz)
        ds=shearyz*(vy[x,y,z+1]-vy[x,y,z]+vz[x,y+1,z]-vz[x,y,z])
        syz[x,y,z]=syz[x,y,z]+ds*ts
        
# %%
def JBUV(x,y,z):
    
    #Vx Cases
    #x in middle
    try:
        if matProps3[x,y,z] in [0,1,2,3,4,9,10,13,14]:
            dvxConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x+1,y,z]))
            dv=dvxConst*( sxx[x+1,y,z]-sxx[x,y,z]
                         +sxy[x,y,z]-sxy[x,y-1,z]
                         +sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts
        #x at 0
        elif matProps3[x,y,z] in [5,7,11,15,16,19,20,21,22]:
            vx[x,y,z] = 2* ((sxx[x+1,y,z])/(matProps0[x,y,z] * gs)) * ts
        # x at xmax
        elif matProps3[x,y,z] in [6,8,12,17,18,23,24,25,26]:
            vx[x,y,z] -= 2 * ((sxx[x,y,z])/(matProps0[x,y,z] * gs))*ts
        #outside sim space
        elif matProps3[x,y,z]  in [99]:
            vx[x,y,z] = 0
        #error handling
        else:
            print('Unrecognized BC x', matProps3[x,y,z],x,y,z)
    except:
        vx[x,y,z] = 0
    
    #Vy cases
    try:
        #y in middle
        if matProps3[x,y,z] in [0,1,2,5,6,7,8,11,12]:
            dvyConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y+1,z]))
            dv=dvyConst* ( sxy[x,y,z]-sxy[x-1,y,z]
                          +syy[x,y+1,z]-syy[x,y,z]
                          +syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts
        #y = 0
        elif matProps3[x,y,z] in [3,9,13,15,17,19,20,23,24]:
            vy[x,y,z] += 2* ((syy[x,y+1,z])/(matProps0[x,y,z] * gs)) * ts
        #y = ymax
        elif matProps3[x,y,z] in [4,10,14,16,18,21,22,25,26]:
            vy[x,y,z] -= 2 * ((syy[x,y,z])/(matProps0[x,y,z] * gs))*ts
        #outside of sim space
        elif matProps3[x,y,z]  in [99]:
            vy[x,y,z] = 0
        #error handling
        else:
            print('Unrecognized BC y', matProps3[x,y,z],x,y,z)
    except:
        vy[x,y,z] = 0

    #Vz cases
    try:
        #z in the middle
        if matProps3[x,y,z] in [0,3,4,5,6,15,16,17,18]:
            dvzConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y,z+1]))
            dv=dvzConst*( sxz[x,y,z]-sxz[x-1,y,z]
                         +syz[x,y,z]-syz[x,y-1,z]
                         +szz[x,y,z+1]-szz[x,y,z])
            vz[x,y,z]=vz[x,y,z]+dv*ts
        #z at 0
        elif matProps3[x,y,z] in [1,7,8,9,10,19,21,23,25]:
            vz[x,y,z] += 2* ((szz[x,y,z+1])/(matProps0[x,y,z] * gs)) * ts
        #z at zmax
        elif matProps3[x,y,z] in [2,11,12,13,14,20,22,24,26]:
            vz[x,y,z] -= 2 * ((szz[x,y,z])/(matProps0[x,y,z] * gs))*ts
        #outside sim space
        elif matProps3[x,y,z] in [99]:
            vz[x,y,z] = 0
        #error handling
        else:
            print('Unrecognized BC z', matProps3[x,y,z],x,y,z)
    except:
        vz[x,y,z] = 0

def setSimSpaceBC99(matPropsglob):
    #matPropsglob[3,:,:,zmax]=99
    #matPropsglob[3,:,:,0]=99
    #matPropsglob[3,:,0,:]=99
    #matPropsglob[3,:,ymax,:]=99
    #matPropsglob[3,0,:,:] = 99
    #matPropsglob[3,xmax,:,:] = 99
    
    matPropsglob[3,1:xmax,1:ymax,0] = 1
    matPropsglob[3,1:xmax,0,1:zmax] = 3
    matPropsglob[3,0,1:ymax,1:zmax] = 5
    
    matPropsglob[3,1:xmax,1:ymax,1] = 1
    matPropsglob[3,1:xmax,1,1:zmax] = 3
    matPropsglob[3,1,1:ymax,1:zmax] = 5
    
    matPropsglob[3,1:xmax,1:ymax,zmax] = 2
    matPropsglob[3,1:xmax,ymax,1:zmax] = 4
    matPropsglob[3,xmax,1:ymax,1:zmax] = 6

    # edges
    # front bottom 
    matPropsglob[3,0,1:ymax,0]=7
    # back bottom 
    matPropsglob[3,xmax,:,0]=8
    # bottom left 
    matPropsglob[3,1:xmax,0,0]=9 
    # bottom right 
    matPropsglob[3,:,ymax,0]=10
    # front top 
    matPropsglob[3,0,:,zmax]=11
    # back top 
    matPropsglob[3,xmax,:,zmax]=12
    # top left 
    matPropsglob[3,:,0,zmax]=13 
    # top rigight 
    matPropsglob[3,:,ymax,zmax]=14
    # front left 
    matPropsglob[3,0,0,1:zmax]=15
    # front right 
    matPropsglob[3,0,ymax,:]=16
    # back left 
    matPropsglob[3,xmax,0,:]=17
    # back right 
    matPropsglob[3,xmax,ymax,:]=18
    ## Corners
    # bottom left front 
    matPropsglob[3,0,0,0]=19
    # top left front 
    matPropsglob[3,0,0,zmax]=20
    # bottom right front 
    matPropsglob[3,0,ymax,0]=21
    # top right front 
    matPropsglob[3,0,ymax,zmax]=22
    # bottom left back 
    matPropsglob[3,xmax,0,0]=23
    # top left back 
    matPropsglob[3,xmax,0,zmax]=24
    # bottom right back 
    matPropsglob[3,xmax,zmax,0]=25
    # top right back 
    matPropsglob[3,xmax,ymax,zmax]=26
    
    return matPropsglob    

matPropsglob = setSimSpaceBC99(matPropsglob)    
if RailShape:
    matPropsglob = setAirCut(matPropsglob)
    matPropsglob = setRailBCs(matPropsglob)
    matPropsglob = addTies(matPropsglob,2)


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

amp=100
decayRate= 0
sinConst=ts*amp/rho1

sinInputSignal=sinConst*np.sin(2*np.pi*frequency*timeVec)*np.exp(-decayRate*timeVec)
sinInputSignal[int(.1*Tsteps+1):] = 0

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
signalloc = np.zeros((npx,gw1,gh1))

mpi_comm.Scatterv([matPropsglob[0,:,:,:],split,offset,MPI.DOUBLE], matProps0)
mpi_comm.Scatterv([matPropsglob[1,:,:,:],split,offset,MPI.DOUBLE], matProps1)
mpi_comm.Scatterv([matPropsglob[2,:,:,:],split,offset,MPI.DOUBLE], matProps2)
mpi_comm.Scatterv([matPropsglob[3,:,:,:],split,offset,MPI.DOUBLE], matProps3)
mpi_comm.Scatterv([signalLocation[:,:,:],split,offset,MPI.DOUBLE], signalloc)


matProps0=distBox(matProps0,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
matProps1=distBox(matProps1,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
matProps2=distBox(matProps2,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
matProps3=distBox(matProps3,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
signalloc=distBox(signalloc,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        

#Now slab has local versions with ghosts of matProps
if (myid == 0) :
    print('split matprops to globs, scratttered parameters to processors, line 449')



stime = time.time()

if (myid == 0 ):
    print('subs setup, line 1213.  About to start at ' + str(stime))

#MidMatrix = np.zeros((gl1,Tsteps))
MidMatrixX = np.zeros((gl1,Tsteps))
MidMatrixY = np.zeros((gl1,Tsteps))
MidMatrixZ = np.zeros((gl1,Tsteps))

Movements = np.zeros((gl1,gw1,int(Tsteps/5)))
DisX = np.zeros((gl1,gw1,gh1))
DisY = np.zeros((gl1,gw1,gh1))
DisZ = np.zeros((gl1,gw1,gh1))

inner = []
outer=[]
for x in range(1,npx+1):
    for y in range(gw1):
        for z in range(gh1):
            if matProps3[x,y,z] == 0:
                inner.append([x,y,z])
            else:
                outer.append([x,y,z])
                
for t in range(0,Tsteps):
    if FFunction ==2:
        vz += signalloc * sinInputSignal[t]
    elif FFunction ==3:
        vx += signalloc * sinInputSignal[t]
    elif FFunction ==4:
        vx += signalloc * sinInputSignal[t]


    for pt in inner:
        #updateStress(pt[0],pt[1],pt[2])
        JBSU(pt[0],pt[1],pt[2])
    for pt in outer:
        #updateStress(pt[0],pt[1],pt[2])
        JBSU(pt[0],pt[1],pt[2])

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

    #if the forcing function is a stress
    if FFunction == 1:
        szz -= signalloc * specificWheelLoad


    for pt in inner:
        #updateVelocity(pt[0],pt[1],pt[2])
        JBUV(pt[0],pt[1],pt[2])
    for pt in outer:
        #updateVelocity(pt[0],pt[1],pt[2])
        JBUV(pt[0],pt[1],pt[2])
        
        
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

    vxg = np.zeros((gl1,gw1,gh1))
    vxt=vx[1:npx+1,:,:]        
    mpi_comm.Gatherv(vxt,[vxg,split,offset,MPI.DOUBLE])

    vzg = np.zeros((gl1,gw1,gh1))
    vzt=vz[1:npx+1,:,:]        
    mpi_comm.Gatherv(vzt,[vzg,split,offset,MPI.DOUBLE])

    vyg = np.zeros((gl1,gw1,gh1))
    vyt=vy[1:npx+1,:,:]        
    mpi_comm.Gatherv(vyt,[vyg,split,offset,MPI.DOUBLE])

    
    if myid==0:
        MidMatrixX[:,t] = vxg[:,inputy,inputz]
        MidMatrixY[:,t] = vyg[:,inputy,inputz]
        MidMatrixZ[:,t] = vzg[:,inputy,inputz]

        #MidMatrixX.append(vxg[:,inputy,inputz])
        #MidMatrixY.append(vxg[:,inputy,inputz])
        #MidMatrixZ.append(vxg[:,inputy,inputz])
        
        DisX += vxg[:,:,:] * ts
        DisY += vyg[:,:,:] * ts
        DisZ += vzg[:,:,:] * ts

        if t%5==0:
            Movements[:,:,int(t/5)-1] = np.sqrt(DisX[:,:,2]**2 + DisY[:,:,2]**2 + DisZ[:,:,2]**2)
            
            fig=plt.figure()
            plt.contourf(np.transpose(vzg[3:,:,int(gh1/2)]), cmap='seismic')
            plt.savefig(imFolder+'Mid/vzWeb'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)
            
            fig=plt.figure()
            plt.contourf(np.transpose(vzg[3:,int(gw1/2),:]), cmap='seismic')
            plt.savefig(imFolder + 'Vert/vzVertCut'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)    
            
            fig=plt.figure()
            plt.contourf(np.transpose(vzg[int(gl1/2),:,:]), cmap='seismic')
            plt.savefig(imFolder + 'Head/vzHead'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)  
               
            fig=plt.figure()
            plt.contourf(np.transpose(vxg[3:,:,int(gh1/2)]), cmap='seismic')
            plt.savefig(imFolder+'Mid/vxWeb'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)
            
            fig=plt.figure()
            plt.contourf(np.transpose(vxg[3:,int(gw1/2),:]), cmap='seismic')
            plt.savefig(imFolder + 'Vert/vxVertCut'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)    
            
            
            fig=plt.figure()
            plt.contourf(np.transpose(vxg[3:,:,int(gh1/2)]), cmap='seismic')
            plt.savefig(imFolder + 'Mid/MidZ'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)  

            fig=plt.figure()
            plt.contourf(np.transpose(np.sqrt(DisX[3:,:,zmax-1]**2 + DisY[3:,:,zmax-1]**2 + DisZ[3:,:,zmax-1]**2)), cmap='seismic')
            plt.savefig(imFolder+'TopSurface/TopSurface'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)
            
            fig=plt.figure()
            plt.contourf(np.transpose(np.sqrt(DisX[3:,ymax-1,:]**2 + DisY[3:,ymax-1,:]**2 + DisZ[3:,ymax-1,:]**2)), cmap='seismic')
            plt.savefig(imFolder+'RightSurface/RightSurface'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)
            
            fig=plt.figure()
            plt.contourf(np.transpose(np.sqrt(DisX[3:,2,:]**2 + DisY[3:,2,:]**2 + DisZ[3:,2,:]**2)), cmap='seismic')
            plt.savefig(imFolder+'LeftSurface/LeftSurface'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
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
        print(t,'/',Tsteps-1,'checksums vx, sxx:',ckvs,ckss, (time.time()-stime)/60.0)
    sys.stdout.flush()
'''
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
'''

if myid ==0:
    #print(MidMatrix)
    
    #MidMatrixX = np.matrix(MidMatrixX)
    #MidMatrixY = np.matrix(MidMatrixY)
    #MidMatrixZ = np.matrix(MidMatrixZ)
    
    if np.shape(MidMatrixX)[0] == Tsteps:
        MidMatrixX = MidMatrixX.T
    
    MidDisplace = np.zeros(np.shape(MidMatrixX))
    
    for i in range(np.shape(MidMatrixX)[0]):
        for j in range(np.shape(MidMatrixX)[1]):
            if j == 0:
                MidDisplace[i,j]=MidMatrixX[i,j]*ts
            else:
                MidDisplace[i,j]=MidDisplace[i,j-1]+MidMatrixX[i,j]*ts
      
    pts = 8
    rng = int(gl1/pts)-1
    
    print(pts, rng)
    
    fig = plt.figure(dpi=600, figsize=(6,4))
    #for i in range(pts):
    plt.plot(MidMatrixX[0,:],label=str('x0'))
    plt.plot(MidMatrixY[0,:],label=str('y0'))
    plt.plot(MidMatrixZ[0,:],label=str('z0'))
    plt.plot(MidMatrixX[50,:],label=str('x50'))
    plt.plot(MidMatrixY[50,:],label=str('y50'))
    plt.plot(MidMatrixZ[50,:],label=str('z50'))
    #    print(str(i*rng)
    plt.title('Velocity')
    plt.legend()
    plt.savefig(imFolder+runName+'MidVelocities.png')
    
    plt.close(fig)
    fig = plt.figure(dpi=600, figsize=(6,4))
    for i in range(pts):
        plt.plot(MidDisplace[i*rng,:],label=str(i*rng))
    plt.legend()
    plt.title('Displacement')
    plt.savefig(imFolder+runName+'MidDisplacements.png')

    plt.close(fig)
    fig = plt.figure(dpi=600, figsize=(6,4))
    for i in range(pts):
        plt.plot(MidMatrixX[i*rng,:],label=str(i*rng))
    plt.legend()
    plt.title('Displacement')
    plt.savefig(imFolder+runName+'MidVel.png')

    print(np.shape(MidMatrixX), np.shape(MidDisplace))
    
    EMin = np.min(Movements[2:,:,8,:134])
    EMax = np.max(Movements[2:,:,8,:134])

    v = np.linspace(EMin, EMax, 15, endpoint=True)

    for i in range(np.shape(Movements)[2]):
        #plt.contour(xi, yi, topSurface[:,:,t].T, v, linewidths=0.5, colors='k')
        plt.contourf(Movements[3:,:,t].T, v, cmap=plt.cm.jet)
        plt.savefig(imFolder+'Energy/Energy'+str(t).zfill(5)+'.png')