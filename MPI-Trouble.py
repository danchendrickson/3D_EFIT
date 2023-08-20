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
length1 = 5
width1 = 0.2
height1 = 0.2

#Image Folder
imFolder = '/sciclone/scr10/dchendrickson01/EFIT/'
runName = 'FixVelGraph'

#is the rail supported by 0, 1 or 2 ties
Ties = 0

#Choose ferquency to be used for excitment
#frequency = 64000
#frequency = 16300
frequency = 8000

#Run for 4 Cycles:
runtime = 12 / frequency 

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
    imFolder += 'CubeS/'

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

Tsteps = int(math.ceil(runtime / ts)) + 1       #total Time Steps

#number of grid points
gl1 = int(math.ceil(length1 / gs)) +1       #length 
gw1 = int(math.ceil(width1 / gs)) +1       #width
gh1 = int(math.ceil(height1 / gs)) +1       #height

print('runtime, gs, ts, gl, gw, gh, Tsteps : ', runtime, gs, ts, gl1, gw1, gh1, Tsteps)

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
    signalLocation[:3,:,:] = 1


if myid == 0:
    print('globs made, line 145')

#########
# FUnctions

def updateStress(x,y,z):
        
    #Calculate constants for stress equations
    norm1=(1/gs)*(matProps1[x,y,z]+2*matProps2[x,y,z])
    norm2=(1/gs)*(matProps1[x,y,z])


    try:
        shearDenomxy=(1/matProps2[x,y,z])+(1/matProps2[x+1,y,z])+(1/matProps2[x,y+1,z])+(1/matProps2[x+1,y+1,z])
        shearxy=(1/gs)*(4/shearDenomxy)
    except:
        pass
    
    try:
        shearDenomxz=(1/matProps2[x,y,z])+(1/matProps2[x+1,y,z])+(1/matProps2[x,y,z+1])+(1/matProps2[x+1,y,z+1])
        shearxz=(1/gs)*(4/shearDenomxz)
    except:
        pass
    
    try:
        shearDenomyz=(1/matProps2[x,y,z])+(1/matProps2[x,y+1,z])+(1/matProps2[x,y,z+1])+(1/matProps2[x,y+1,z+1])
        shearyz=(1/gs)*(4/shearDenomyz)
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
            vx[x,y,z] += 2* ((sxx[x+1,y,z])/(matProps0[x,y,z] * gs)) * ts
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
        
# %%
def updateVelocity(x,y,z):
    try:
        #body
        if matProps3[x,y,z] == 0:
            dvxConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x+1,y,z]))
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
            vx[x,y,z]=vx[x,y,z]+dv*ts

            dvyConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y+1,z]))
            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts

            dvzConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y,z+1]))
            dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
            vz[x,y,z]=vz[x,y,z]+dv*ts


        #incase non boundaries or unknown areas snuck through
        elif (matProps3[x,y,z] == 99): # or matProps3[x,y,z] == 2 or matProps3[x,y,z] == 4 or matProps3[x,y,z] == 6):
            vx[x,y,z]=0
            vy[x,y,z]=0
            vz[x,y,z]=0

        #faces    
        elif matProps3[x,y,z] == 1:
            dvxConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x+1,y,z]))
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]) # z = 0, and sxz is set to 0, so remove this bit+sxz[x,y,z+1]-sxz[x,y,z])
            vx[x,y,z]=vx[x,y,z]+dv*ts

            dvyConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y+1,z]))
            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]) # not at z=0+syz[x,y,z]-syz[x,y,z-1])
            vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z] += 2* ((szz[x,y,z+1])/(matProps0[x,y,z] * gs)) * ts

        elif matProps3[x,y,z] == 2:
            dvxConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x+1,y,z]))
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]) #took out at 0, try at max +sxz[x,y,z]-sxz[x,y,z-1]
            #vx[x,y,z]=vx[x,y,z]+dv*ts

            dvyConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y+1,z]))
            dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]) #same +syz[x,y,z]-syz[x,y,z-1]
            #vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z] -= 2 * ((szz[x,y,z])/(matProps0[x,y,z] * gs))*ts

        elif matProps3[x,y,z] == 3:
            dvxConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x+1,y,z]))
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxz[x,y,z]-sxz[x,y,z-1]) #not at y=0  +sxy[x,y,z]-sxy[x,y-1,z]
            vx[x,y,z]=vx[x,y,z]+dv*ts

            vy[x,y,z] += 2* ((syy[x,y+1,z])/(matProps0[x,y,z] * gs)) * ts

            dvzConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y,z+1]))
            dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+szz[x,y,z+1]-szz[x,y,z]) # not at y=0+syz[x,y,z]-syz[x,y-1,z]
            vz[x,y,z]=vz[x,y,z]+dv*ts


        elif matProps3[x,y,z] == 4:
            dvxConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x+1,y,z]))
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxz[x,y,z]-sxz[x,y,z-1]) # took out at 0 try at max +sxy[x,y,z]-sxy[x,y-1,z]
            #vx[x,y,z]=vx[x,y,z]+dv*ts

            vy[x,y,z] -= 2 * ((syy[x,y,z])/(matProps0[x,y,z] * gs))*ts

            dvzConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y,z+1]))
            dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+szz[x,y,z+1]-szz[x,y,z]) # same +syz[x,y,z]-syz[x,y-1,z]
            #vz[x,y,z]=vz[x,y,z]+dv*ts     
            
        elif matProps3[x,y,z] == 5:
            vx[x,y,z] += 2* ((sxx[x+1,y,z])/(matProps0[x,y,z] * gs)) * ts

            dvyConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y+1,z]))
            dv=dvyConst*(syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1]) #not at x=0sxy[x,y,z]-sxy[x-1,y,z]+
            vy[x,y,z]=vy[x,y,z]+dv*ts

            dvzConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y,z+1]))
            dv=dvzConst*(syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z]) # not at x=0 sxz[x,y,z]-sxz[x-1,y,z]+
            vz[x,y,z]=vz[x,y,z]+dv*ts

        elif matProps3[x,y,z] == 6:
            vx[x,y,z] -= 2 * ((sxx[x,y,z])/(matProps0[x,y,z] * gs))*ts
            
            dvyConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y+1,z]))
            dv=dvyConst*(syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1]) #took out at 0, try at max sxy[x,y,z]-sxy[x-1,y,z]+
            #vy[x,y,z]=vy[x,y,z]+dv*ts

            dvzConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y,z+1]))
            dv=dvzConst*(syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z]) # sames sxz[x,y,z]-sxz[x-1,y,z]+
            #vz[x,y,z]=vz[x,y,z]+dv*ts
               


        #EDGES
        #bottom edges
        elif matProps3[x,y,z] == 7:
            vx[x,y,z] += 2* ((sxx[x+1,y,z])/(matProps0[x,y,z] * gs)) * ts

            dvyConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y+1,z]))
            dv=dvyConst*(syy[x,y+1,z]-syy[x,y,z]) #x and y = 0 sxy[x,y,z]-sxy[x-1,y,z]++syz[x,y,z]-syz[x,y,z-1]
            vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z] += 2* ((szz[x,y,z+1])/(matProps0[x,y,z] * gs)) * ts
            
        elif matProps3[x,y,z] == 8:
            vx[x,y,z] -= 2 * ((sxx[x,y,z])/(matProps0[x,y,z] * gs))*ts
            
            dvyConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y+1,z]))
            dv=dvyConst*(syy[x,y+1,z]-syy[x,y,z]) #sxy[x,y,z]-sxy[x-1,y,z]+ +syz[x,y,z]-syz[x,y,z-1]
            vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z] += 2* ((szz[x,y,z+1])/(matProps0[x,y,z] * gs)) * ts
            
        elif matProps3[x,y,z] == 9:
            dvxConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x+1,y,z]))
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]) # +sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1]
            vx[x,y,z]=vx[x,y,z]+dv*ts

            vy[x,y,z] += 2* ((syy[x,y+1,z])/(matProps0[x,y,z] * gs)) * ts
            
            vz[x,y,z] += 2* ((szz[x,y,z+1])/(matProps0[x,y,z] * gs)) * ts
            
        elif matProps3[x,y,z] == 10:
            dvxConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x+1,y,z]))
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]) # +sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1]
            vx[x,y,z]=vx[x,y,z]+dv*ts

            vy[x,y,z] -= 2 * ((syy[x,y,z])/(matProps0[x,y,z] * gs))*ts
            
            vz[x,y,z] += 2* ((szz[x,y,z+1])/(matProps0[x,y,z] * gs)) * ts
            
        #top edges
        elif matProps3[x,y,z] == 11:
            vx[x,y,z] += 2* ((sxx[x+1,y,z])/(matProps0[x,y,z] * gs)) * ts
            
            dvyConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y+1,z]))
            dv=dvyConst*(syy[x,y+1,z]-syy[x,y,z]) #sxy[x,y,z]-sxy[x-1,y,z]+ +syz[x,y,z]-syz[x,y,z-1]
            vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z] -= 2 * ((szz[x,y,z])/(matProps0[x,y,z] * gs))*ts

        elif matProps3[x,y,z] == 12:
            vx[x,y,z] -= 2 * ((sxx[x,y,z])/(matProps0[x,y,z] * gs))*ts
            
            dvyConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y+1,z]))
            dv=dvyConst*(syy[x,y+1,z]-syy[x,y,z]) #sxy[x,y,z]-sxy[x-1,y,z]+ +syz[x,y,z]-syz[x,y,z-1]
            vy[x,y,z]=vy[x,y,z]+dv*ts

            vz[x,y,z] -= 2 * ((szz[x,y,z])/(matProps0[x,y,z] * gs))*ts
            
            
        elif matProps3[x,y,z] == 13:
            dvxConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x+1,y,z]))
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]) # +sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1]
            vx[x,y,z]=vx[x,y,z]+dv*ts

            vy[x,y,z] += 2* ((syy[x,y+1,z])/(matProps0[x,y,z] * gs)) * ts
            
            vz[x,y,z] -= 2 * ((szz[x,y,z])/(matProps0[x,y,z] * gs))*ts
            

        elif matProps3[x,y,z] == 14:
            dvxConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x+1,y,z]))
            dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]) # +sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1]
            vx[x,y,z]=vx[x,y,z]+dv*ts

            vy[x,y,z] -= 2 * ((syy[x,y,z])/(matProps0[x,y,z] * gs))*ts
            
            vz[x,y,z] -= 2 * ((szz[x,y,z])/(matProps0[x,y,z] * gs))*ts

        #side edges
        elif matProps3[x,y,z] == 15:
            vx[x,y,z] += 2* ((sxx[x+1,y,z])/(matProps0[x,y,z] * gs)) * ts
            
            vy[x,y,z] += 2* ((syy[x,y+1,z])/(matProps0[x,y,z] * gs)) * ts
            
            dvzConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y,z+1]))
            dv=dvzConst*(szz[x,y,z+1]-szz[x,y,z]) # sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+
            vz[x,y,z]=vz[x,y,z]+dv*ts


        elif matProps3[x,y,z] == 16:
            vx[x,y,z] += 2* ((sxx[x+1,y,z])/(matProps0[x,y,z] * gs)) * ts

            vy[x,y,z] -= 2 * ((syy[x,y,z])/(matProps0[x,y,z] * gs))*ts
            
            dvzConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y,z+1]))
            dv=dvzConst*(szz[x,y,z+1]-szz[x,y,z]) # sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+
            vz[x,y,z]=vz[x,y,z]+dv*ts

        elif matProps3[x,y,z] == 17:
            vx[x,y,z] -= 2 * ((sxx[x,y,z])/(matProps0[x,y,z] * gs))*ts
            
            vy[x,y,z] += 2* ((syy[x,y+1,z])/(matProps0[x,y,z] * gs)) * ts
            
            dvzConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y,z+1]))
            dv=dvzConst*(szz[x,y,z+1]-szz[x,y,z]) # sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+
            vz[x,y,z]=vz[x,y,z]+dv*ts 



        elif matProps3[x,y,z] == 18:
            vx[x,y,z] -= 2 * ((sxx[x,y,z])/(matProps0[x,y,z] * gs))*ts
            
            vy[x,y,z] -= 2 * ((syy[x,y,z])/(matProps0[x,y,z] * gs))*ts
            
            dvzConst=2*(1/gs)*(1/(matProps0[x,y,z]+matProps0[x,y,z+1]))
            dv=dvzConst*(szz[x,y,z+1]-szz[x,y,z]) # sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+
            vz[x,y,z]=vz[x,y,z]+dv*ts



        #CORNERS
        elif matProps3[x,y,z] == 19:
            vx[x,y,z] += 2* ((sxx[x+1,y,z])/(matProps0[x,y,z] * gs)) * ts
            vy[x,y,z] += 2* ((syy[x,y+1,z])/(matProps0[x,y,z] * gs)) * ts
            vz[x,y,z] += 2* ((szz[x,y,z+1])/(matProps0[x,y,z] * gs)) * ts
            
        elif matProps3[x,y,z] == 20:
            vx[x,y,z] += 2* ((sxx[x+1,y,z])/(matProps0[x,y,z] * gs)) * ts
            vy[x,y,z] += 2* ((syy[x,y+1,z])/(matProps0[x,y,z] * gs)) * ts
            vz[x,y,z] -= 2 * ((szz[x,y,z])/(matProps0[x,y,z] * gs))*ts

        elif matProps3[x,y,z] == 21:
            vx[x,y,z] += 2* ((sxx[x+1,y,z])/(matProps0[x,y,z] * gs)) * ts
            vy[x,y,z] -= 2 * ((syy[x,y,z])/(matProps0[x,y,z] * gs))*ts
            vz[x,y,z] += 2* ((szz[x,y,z+1])/(matProps0[x,y,z] * gs)) * ts

        elif matProps3[x,y,z] == 22:
            vx[x,y,z] += 2* ((sxx[x+1,y,z])/(matProps0[x,y,z] * gs)) * ts
            vy[x,y,z] -= 2 * ((syy[x,y,z])/(matProps0[x,y,z] * gs))*ts
            vz[x,y,z] -= 2 * ((szz[x,y,z])/(matProps0[x,y,z] * gs))*ts

        elif matProps3[x,y,z] == 23:
            vx[x,y,z] -= 2 * ((sxx[x,y,z])/(matProps0[x,y,z] * gs))*ts
            vy[x,y,z] += 2* ((syy[x,y+1,z])/(matProps0[x,y,z] * gs)) * ts
            vz[x,y,z] += 2* ((szz[x,y,z+1])/(matProps0[x,y,z] * gs)) * ts

        elif matProps3[x,y,z] == 24:
            vx[x,y,z] -= 2 * ((sxx[x,y,z])/(matProps0[x,y,z] * gs))*ts
            vy[x,y,z] += 2* ((syy[x,y+1,z])/(matProps0[x,y,z] * gs)) * ts
            vz[x,y,z] -= 2 * ((szz[x,y,z])/(matProps0[x,y,z] * gs))*ts

        elif matProps3[x,y,z] == 25:
            vx[x,y,z] -= 2 * ((sxx[x,y,z])/(matProps0[x,y,z] * gs))*ts
            vy[x,y,z] -= 2 * ((syy[x,y,z])/(matProps0[x,y,z] * gs))*ts
            vz[x,y,z] += 2* ((szz[x,y,z+1])/(matProps0[x,y,z] * gs)) * ts

        elif matProps3[x,y,z] == 26:
            vx[x,y,z] -= 2 * ((sxx[x,y,z])/(matProps0[x,y,z] * gs))*ts
            vy[x,y,z] -= 2 * ((syy[x,y,z])/(matProps0[x,y,z] * gs))*ts
            vz[x,y,z] -= 2 * ((szz[x,y,z])/(matProps0[x,y,z] * gs))*ts

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
            vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
            vz[x,y,z]=0


        else: print('error: ',x,y,z, matProps3[x,y,z])
    except:
        print('Boundary Conditon isssue Velocity: ', str(x), str(y), str(z), str(matProps3[x,y,z]))
              
        
def setAirCut(matPropsglob):
        
    #zone 1 of air, left of web
    for yy in range(gridStartWeb):
        y = yy + 0
        for zz in range(gridStartHead - gridEndFoot):
            z = zz + gridEndFoot + 1
            matPropsglob[0,:,y,z]=rho2
            matPropsglob[1,:,y,z]=lmbda2
            matPropsglob[2,:,y,z]=mu2
            matPropsglob[3,:,y,z]=99
            
    # zone 2 of air left of head
    for yy in range(gridStartHeadWidth):
        y = yy + 0
        for zz in range(gridStartHead):
            z = zmax - zz
            matPropsglob[0,:,y,z]=rho2
            matPropsglob[1,:,y,z]=lmbda2
            matPropsglob[2,:,y,z]=mu2
            matPropsglob[3,:,y,z]=99
            
    # zone 3 of air, right of web
    for yy in range(gridStartWeb):
        y = ymax - yy
        for zz in range(gridStartHead - gridEndFoot):
            z = zz + gridEndFoot + 1
            matPropsglob[0,:,y,z]=rho2
            matPropsglob[1,:,y,z]=lmbda2
            matPropsglob[2,:,y,z]=mu2
            matPropsglob[3,:,y,z]=99
            
    # zone 4 of air, right of head
    for yy in range(gridStartHeadWidth):
        y = ymax - yy
        for zz in range(gridStartHead):
            z = zmax - zz
            matPropsglob[0,:,y,z]=rho2
            matPropsglob[1,:,y,z]=lmbda2
            matPropsglob[2,:,y,z]=mu2
            matPropsglob[3,:,y,z]=99
            
    return matPropsglob

def setSimSpaceBCs(matPropsglob):
    # Set the simulations pace boundarys
    
    #faces
    # top
    matPropsglob[3,:,:,zmax]=2

    # bottom 
    matPropsglob[3,:,:,0]=1

    # left
    matPropsglob[3,:,0,:]=3

    # right
    matPropsglob[3,:,ymax,:]=4
    
    # Front
    matPropsglob[3,0,:,:] = 5
    
    # back
    matPropsglob[3,xmax,:,:] = 6

    # edges
    # top left
    matPropsglob[3,:,0,zmax]=13

    # top rigight
    matPropsglob[3,:,ymax,zmax]=14
    
    # bottom left
    matPropsglob[3,:,0,0]=9
    
    # bottom right
    matPropsglob[3,:,ymax,0]=10
    
    # front left
    matPropsglob[3,0,0,:]=15
    
    # front right
    matPropsglob[3,0,ymax,:]=16
    
    # back left
    matPropsglob[3,xmax,0,:]=17
    
    # back right
    matPropsglob[3,xmax,ymax,:]=18
    
    # front top
    matPropsglob[3,0,:,zmax]=11
    
    # front bottom
    matPropsglob[3,0,:,0]=7
    
    # back top
    matPropsglob[3,xmax,:,zmax]=12
    
    # back bottom
    matPropsglob[3,xmax,:,0]=8
    
    ## Corners
    # top left front
    matPropsglob[3,0,0,zmax]=20
    
    # top right front
    matPropsglob[3,0,ymax,zmax]=22
    
    # top left back
    matPropsglob[3,xmax,0,zmax]=24
    
    # top right back
    matPropsglob[3,xmax,ymax,zmax]=26
    
    # bottom left front
    matPropsglob[3,0,0,0]=19
    
    # bottom right front
    matPropsglob[3,0,ymax,0]=21
    
    # bottom left back
    matPropsglob[3,xmax,0,0]=23
    
    # bottom right back
    matPropsglob[3,xmax,zmax,0]=25

    return matPropsglob
    
def setRailBCs(matPropsglob):
    #set the boundary conditions in material props4
    # Set the simulations pace boundarys
    # top
    
    #top of footing
    z = gridEndFoot
    for yy in range(gridStartWeb):
        y = yy + 0
        matPropsglob[3,:,y,z]=2
        y = ymax - yy
        matPropsglob[3,:,y,z]=2

    #sides of web
    yl=gridStartWeb
    yr=gridEndWeb
    for zz in range(gridStartHead - gridEndFoot):
        z = zz + gridEndFoot
        matPropsglob[3,:,yl,z]=3
        matPropsglob[3,:,yr,z]=4

    #bottom of head
    z=gridStartHead
    for yy in range(gridStartWeb - gridStartHeadWidth+1):
        y=yy+gridStartHeadWidth
        matPropsglob[3,:,y,z]=1
        y=gridEndHeadWidth-yy
        matPropsglob[3,:,y,z]=1

    #sides of head
    yl=gridStartHeadWidth
    yr=gridEndHeadWidth
    for zz in range(zmax - gridStartHead):
        z = zmax - zz
        matPropsglob[3,:,yl,z]=3
        matPropsglob[3,:,yr,z]=4

    #top edge of foot on left
    matPropsglob[3,:,0,gridEndFoot]=13
    #top endge of foot on right
    matPropsglob[3,:,ymax,gridEndFoot]=14

    #bottom of head on left
    matPropsglob[3,:,gridStartHeadWidth,gridStartHead]=9
    #bottom of head on right
    matPropsglob[3,:,gridEndHeadWidth,gridStartHead]=10

    #Top of head on left
    matPropsglob[3,:,gridStartHeadWidth,zmax]=13
    #top of head on right
    matPropsglob[3,:,gridEndHeadWidth,zmax]=14
    
    ## Special cases for front and back face
    #bottom of head
    z=gridStartHead
    for yy in range(gridStartWeb - gridStartHeadWidth+1):
        y=yy+gridStartHeadWidth
        matPropsglob[3,0,y,z]=7
        matPropsglob[3,xmax,y,z]=8
        y=gridEndHeadWidth-yy
        matPropsglob[3,0,y,z]=7
        matPropsglob[3,xmax,y,z]=8
    #top of foot
    z = gridEndFoot
    for yy in range(gridStartWeb):
        y = yy + 0
        matPropsglob[3,0,y,z]=11
        matPropsglob[3,xmax,y,z]=12
        y = ymax - yy
        matPropsglob[3,0,y,z]=11
        matPropsglob[3,xmax,y,z]=12    
    #sides of web
    yl=gridStartWeb
    yr=gridEndWeb
    for zz in range(gridStartHead - gridEndFoot):
        z = zz + gridEndFoot
        matPropsglob[3,0,yl,z]=15
        matPropsglob[3,0,yr,z]=16
        matPropsglob[3,xmax,yl,z]=17
        matPropsglob[3,xmax,yr,z]=18
    #sides of head
    yl=gridStartHeadWidth
    yr=gridEndHeadWidth
    for zz in range(zmax - gridStartHead):
        z = zmax - zz
        matPropsglob[3,0,yl,z]=15
        matPropsglob[3,0,yr,z]=16    
        matPropsglob[3,xmax,yl,z]=17
        matPropsglob[3,xmax,yr,z]=18   
    
    #front bottom left head corner
    matPropsglob[3,0,gridStartHeadWidth,gridStartHead]=19
    #front bottom right head corner
    matPropsglob[3,0,gridEndHeadWidth,gridStartHead]=21
    #front top left head corner
    matPropsglob[3,0,gridStartHeadWidth,zmax]=20
    #front top right head corner
    matPropsglob[3,0,gridEndHeadWidth,zmax]=22
    #front left top foot corner
    matPropsglob[3,0,0,gridEndFoot]=20
    #front right top foot corner
    matPropsglob[3,0,ymax,gridEndFoot]=22
    #back bottom left head corner
    matPropsglob[3,xmax,gridStartHeadWidth,gridStartHead]=23
    #back bottom right head corner
    matPropsglob[3,xmax,gridEndHeadWidth,gridStartHead]=25
    #back top left head corner
    matPropsglob[3,xmax,gridStartHeadWidth,zmax]=24
    #back top right head corner
    matPropsglob[3,xmax,gridEndHeadWidth,zmax]=26
    #back left top foot corner
    matPropsglob[3,xmax,0,gridEndFoot]=24
    #back right top foot corner
    matPropsglob[3,xmax,ymax,gridEndFoot]=26
    
        
    return matPropsglob
       


def addTies(matPropsglob, Ties):

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
    
    return matPropsglob
    

matPropsglob = setSimSpaceBCs(matPropsglob)
    
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
MidMatrixX=[]
MidMatrixY=[]
MidMatrixZ=[]

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

    if FFunction == 2:
        vz += signalloc * sinInputSignal[t]
    elif FFunction ==3:
        vz += signalloc * sinInputSignal[t]


    for pt in inner:
        updateStress(pt[0],pt[1],pt[2])
    for pt in outer:
        updateStress(pt[0],pt[1],pt[2])
      

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
        updateVelocity(pt[0],pt[1],pt[2])
    for pt in outer:
        updateVelocity(pt[0],pt[1],pt[2])
        
        
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
        MidMatrixX.append(vxg[:,inputy,inputz])
        MidMatrixY.append(vxg[:,inputy,inputz])
        MidMatrixZ.append(vxg[:,inputy,inputz])
    
        '''if t%10==0:
            fig=plt.figure()
            plt.contourf(np.transpose(vzg[:,:,int(gh1/2)]), cmap='seismic')
            plt.savefig(imFolder+'Mid/vzWeb'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)
            
            fig=plt.figure()
            plt.contourf(np.transpose(vzg[:,int(gw1/2),:]), cmap='seismic')
            plt.savefig(imFolder + 'Vert/vzVertCut'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)    
            
            fig=plt.figure()
            plt.contourf(np.transpose(vzg[int(gl1/2),:,:]), cmap='seismic')
            plt.savefig(imFolder + 'Head/vzHead'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)  
               
            fig=plt.figure()
            plt.contourf(np.transpose(vxg[:,:,int(gh1/2)]), cmap='seismic')
            plt.savefig(imFolder+'Mid/vxWeb'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)
            
            fig=plt.figure()
            plt.contourf(np.transpose(vxg[:,int(gw1/2),:]), cmap='seismic')
            plt.savefig(imFolder + 'Vert/vxVertCut'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)    
            
            fig=plt.figure()
            plt.contourf(np.transpose(vxg[int(gl1/2),:,:]), cmap='seismic')
            plt.savefig(imFolder + 'Head/vxHead'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)  
        '''
    

            

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
    
    MidMatrixX = np.matrix(MidMatrixX)
    MidMatrixY = np.matrix(MidMatrixY)
    MidMatrixZ = np.matrix(MidMatrixZ)
    
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