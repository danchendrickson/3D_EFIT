import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math
import time
import functools
import pickle
    

from distBox import distBox

import sys
from mpi4py import MPI
from os import environ 
import os
from typing import *                                                                       
#MPIComm = Union[MPI.Intracomm, MPI.Intercomm]
mpi_comm = MPI.COMM_WORLD
myid = mpi_comm.Get_rank()                                                         
mpi_size = mpi_comm.Get_size()        
nprocs=mpi_size

if myid==0:
    print("started")
# for overlapping slabs:  
# # points per proc along z = npz = gh1/nproc (+2 to ghost boundaries)
# glob_index = loc_index-1 + npz*myid
# loc_index = glob_index - npz*myid + 1
# myid given glob_index = glob_index/npz = ghloc-2

# set Constants
AirCut = False
RailShape = True
Flaw = False
Absorbing = False
CSVs = False
Pickles = True
Images = False

#Dimmesnsion of simulation space in meters
length1 = 10
width1 = 0.1524 # 0.1524
height1 = 0.1524
Wheel1Distance = 1.5 # wheel starts 1 meter down track

#Image Folder
imFolder = '/sciclone/scr10/dchendrickson01/EFIT/'
runName = '12mLongRail'

#is the rail supported by 0, 1 or 2 ties
Ties = 0

#Choose ferquency to be used for excitment
frequency = 74574.0  #brute forced this number to be where simulation frequency 
#            49720  is 2,000,000, alowing for %10 to equal laser 200k same as actual
#            74574  is 3,000,000 hz running, and sample rate %15 is 200k same as actual, if we need more dense
Signalfrequency = 16300

SaveSize = 250  #experimentally found for where we don't run out of memory.

cycles = 60

figDPI = 600

#Forcing Function Location and type
# 1 for dropped wheel on top toward start
# 2 for rubbing flange on side
# 3 for plane wave from end
# 4 is square patch rail 20% down
# 5 is a dropped wheel on the mid point
# 6 for two rubbing flanges, one on each side

FFunction = 6

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
    imFolder += 'Temp/'
elif FFunction == 3:
    imFolder += 'BiggerAcTii/'
elif FFunction == 4:
    imFolder += 'BiggerAcTii/'
elif FFunction == 5:
    imFolder += 'RailLong/'
elif FFunction == 6:   #long rail, two wheel rubs
    imFolder += '12mFromRFR/'

if myid==0:
    if os.path.isdir(imFolder):
        pass
    else:
        os.makedirs(imFolder)


#Set time step and grid step to be 10 steps per frequency and ten steps per wavelength respectively
#ts = 1 / frequency / 10    #time step
gs = (min(omegaL1, omegaT1) /12)    #grid step, omegaL2,omegaT2
ts = gs/((max(cl1,ct1))*(np.sqrt(3)))*0.95 #time step, cl2,ct2


#change to lower frequency but with dense grid
frequency = 16300

#Run for 3 seconds, what the laser can store:
#runtime = cycles / frequency #cycles / frequency 
#Tsteps = int(math.ceil(runtime / ts)) + 1       #total Time Steps

Tsteps = 6900   #calculated spereately for needed space to get reflections
runtime = Tsteps * ts

#number of grid points
gl1 = int(math.ceil(length1 / gs)) +1       #length 
gw1 = int(math.ceil(width1 / gs)) +1       #width
gh1 = int(math.ceil(height1 / gs)) +1       #height

if (myid == 0) :
    print('runtime, gs, ts, fr, gl, gw, gh, Tsteps: ', runtime, gs, ts, 1/ts, gl1, gw1, gh1, Tsteps)

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
matDensityAll=np.zeros((gl1,gw1,gh1))
matLambdaAll=np.zeros((gl1,gw1,gh1))
matMuAll=np.zeros((gl1,gw1,gh1))
matBCall=np.zeros((gl1,gw1,gh1))
signalLocation=np.zeros((gl1,gw1,gh1))

matDensityAll[:,:,:]=rho1
matLambdaAll[:,:,:]=lmbda1
matMuAll[:,:,:]=mu1
matBCall[:,:,:]=0

## for latter rail section, define the dimmmensions in terms of grid
HeadThickness = 0.05
WebThickness = 0.035
FootThickness = 0.03
HeadWidth = 0.102

relHeadThick = HeadThickness / height1
relWeb = WebThickness / width1
relFoot = FootThickness / height1
relHeadWidth = HeadWidth / width1

relStartHeadThick = 1 - relHeadThick
relStartWeb = 0.5 - (relWeb / 2.0)
relEndWeb = 0.5 + (relWeb / 2.0)
relStartHeadWidth = 0.5 - (relHeadWidth / 2.0)
relEndHeadWidth = 0.5 + (relHeadWidth / 2.0)


gridStartHead = round((gh1-3) * relStartHeadThick) + 1
gridStartWeb = round((gw1-3) * relStartWeb)  + 1
gridEndWeb = round((gw1-3) * relEndWeb)  + 1
gridEndFoot = round((gh1-3) * relFoot)  + 1
gridStartHeadWidth = round((gw1-3) * relStartHeadWidth)  + 1
gridEndHeadWidth = round((gw1-3)  * relEndHeadWidth)  + 1


#Make the Signal Location grid
if FFunction == 1:
    pnodes = max(int(whlayer / 2),3)
    contactLength = max(int(0.001 / gs),3)  #1 cm contact patch or 3 nodes, whichever is larger
    
    #starting at .25 down, to be between the first 2 ties
    WheelStartPoint = int(0.25 * gl1)
    
    signalLocation[WheelStartPoint:WheelStartPoint+contactLength,gridStartHeadWidth:gridEndHeadWidth, -3:] = 1
    
elif FFunction == 2:
     
    signalLocation[14:20,gridStartHeadWidth:gridStartHeadWidth+2,gridStartHead:zmax-2] = 1

    signalLocation[20,gridStartHeadWidth:gridStartHeadWidth+2,gridStartHead:zmax-2] = 0.5
    signalLocation[13,gridStartHeadWidth:gridStartHeadWidth+2,gridStartHead:zmax-2] = 0.5
    signalLocation[14:20,gridStartHeadWidth+2:gridStartHeadWidth+3,gridStartHead:zmax-2] = 0.5
    
elif FFunction == 3:
    signalLocation[2:4,:,:] = 1


elif FFunction == 4:
    start = 2*int(gh1/5)
    end = 3 * int(gh1/5)
    signalLocation[2:4,start:end,start:end] = 1

elif FFunction == 5:

    signalLocation[int(gl1/2)-5:int(gl1/2)+5,int(gw1/2)-5:int(gw1/2)+5,zmax-2] = 1
    signalLocation[int(gl1/2)-5:int(gl1/2)+5,int(gw1/2)-5:int(gw1/2)+5,zmax-1] = 0.5
    signalLocation[int(gl1/2)-5:int(gl1/2)+5,int(gw1/2)-5:int(gw1/2)+5,zmax-3] = 0.5
    

elif FFunction == 6:
    
    Wheel1Start = int(Wheel1Distance / gs)
    
    signalLocation[Wheel1Start:Wheel1Start+6,gridStartHeadWidth:gridStartHeadWidth+2,gridStartHead:zmax-2] = 1

    signalLocation[Wheel1Start+6,gridStartHeadWidth:gridStartHeadWidth+2,gridStartHead:zmax-2] = 0.5
    signalLocation[Wheel1Start-1,gridStartHeadWidth:gridStartHeadWidth+2,gridStartHead:zmax-2] = 0.5
    signalLocation[Wheel1Start:Wheel1Start+6,gridStartHeadWidth+2:gridStartHeadWidth+3,gridStartHead:zmax-2] = 0.5

    sep = int(1.360/gs) # Wheel 2 is centered 1.36 meters from wheel 1
    
    signalLocation[Wheel1Start+sep:Wheel1Start+6+sep,gridEndHeadWidth-2:gridEndHeadWidth,gridStartHead:zmax-2] = 1

    signalLocation[Wheel1Start+6+sep,gridEndHeadWidth-2:gridEndHeadWidth,gridStartHead:zmax-2] = 0.5
    signalLocation[Wheel1Start-1+sep,gridEndHeadWidth-2:gridEndHeadWidth,gridStartHead:zmax-2] = 0.5
    signalLocation[Wheel1Start+sep:Wheel1Start+6+sep,gridEndHeadWidth-3:gridEndHeadWidth-2,gridStartHead:zmax-2] = 0.5
    
    
if myid == 0:
    print('globs made, line 145')

#########
# FUnctions
def JBSU(x,y,z):
    try:
        if (matBCs[x,y,z] == 2 or matBCs[x-1,y,z] == 2 or matBCs[x,y-1,z] == 2 or matBCs[x,y,z-1] == 2):
            pass
        else:
            norm1=(1/gs)*(matLambda[x,y,z]+2*matMu[x,y,z])
            norm2=(1/gs)*(matLambda[x,y,z])

            ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
            sxx[x,y,z]=sxx[x,y,z]+ds*ts

            ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
            syy[x,y,z]=syy[x,y,z]+ds*ts

            ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
            szz[x,y,z]=szz[x,y,z]+ds*ts

        if (matBCs[x,y,z] == 2 or matBCs[x+1,y,z] == 2 or matBCs[x-1,y,z] == 2 or matBCs[x,y+1,z] == 2 
            or matBCs[x,y-1,z] == 2 or matBCs[x,y,z-1] == 2 or matBCs[x+1,y+1,z] == 2):
            pass
        else:
            shearDenomxy=(1/matMu[x,y,z])+(1/matMu[x+1,y,z])+(1/matMu[x,y+1,z])+(1/matMu[x+1,y+1,z])
            shearxy=4*(1/gs)*(1/shearDenomxy)
            ds=shearxy*(vx[x,y+1,z]-vx[x,y,z]+vy[x+1,y,z]-vy[x,y,z])
            sxy[x,y,z]=sxy[x,y,z]+ds*ts

        if (matBCs[x,y,z] == 2 or matBCs[x+1,y,z] == 2 or matBCs[x-1,y,z] == 2 or matBCs[x,y,z+1] == 2 
            or matBCs[x,y,z-1] == 2 or matBCs[x,y-1,z] == 2 or matBCs[x+1,y,z+1] == 2):
            pass
        else:
            shearDenomxz=(1/matMu[x,y,z])+(1/matMu[x+1,y,z])+(1/matMu[x,y,z+1])+(1/matMu[x+1,y,z+1])
            shearxz=4*(1/gs)*(1/shearDenomxz)
            ds=shearxz*(vx[x,y,z+1]-vx[x,y,z]+vz[x+1,y,z]-vz[x,y,z])
            sxz[x,y,z]=sxz[x,y,z]+ds*ts   

        if (matBCs[x,y,z] == 2 or matBCs[x,y,z+1] == 2 or matBCs[x,y,z-1] == 2 or matBCs[x,y+1,z] == 2 
            or matBCs[x,y-1,z] == 2 or matBCs[x-1,y,z] == 2 or matBCs[x,y+1,z+1] == 2):
            pass
        else:
            shearDenomyz=(1/matMu[x,y,z])+(1/matMu[x,y+1,z])+(1/matMu[x,y,z+1])+(1/matMu[x,y+1,z+1])
            shearyz=4*(1/gs)*(1/shearDenomyz)
            ds=shearyz*(vy[x,y,z+1]-vy[x,y,z]+vz[x,y+1,z]-vz[x,y,z])
            syz[x,y,z]=syz[x,y,z]+ds*ts
    except:
        print('Unrecognized BC stress', matBCs[x,y,z],x,y,z)

        
# %%

def JBUV(x,y,z):
    
    if matBCs[x,y,z] == 0: 
        dvxConst=2*(1/gs)*(1/(matDensity[x,y,z]+matDensity[x+1,y,z]))
        dv=dvxConst*( sxx[x+1,y,z]-sxx[x,y,z]
                     +sxy[x,y,z]-sxy[x,y-1,z]
                     +sxz[x,y,z]-sxz[x,y,z-1])
        vx[x,y,z]=vx[x,y,z]+dv*ts
    #x at 0
    elif (matBCs[x,y,z] ==2 or matBCs[x,y-1,z]==2 or matBCs[x,y,z-2]==2):
        pass #requires elements out of the direction
    elif matBCs[x+1,y,z] == 2:
        vx[x,y,z] += 2 * ts/gs * 1/(2 * matDensity[x,y,z]) * ((-2)*sxx[x,y,z])

    elif matBCs[x-1,y,z] ==2 :
        vx[x,y,z] += 2 * ts/gs * 1/(2 * matDensity[x,y,z]) * ((2)*sxx[x+1,y,z])

    else:
        dvxConst=2*(1/gs)*(1/(matDensity[x,y,z]+matDensity[x+1,y,z]))
        dv=dvxConst*( sxx[x+1,y,z]-sxx[x,y,z]
                     +sxy[x,y,z]-sxy[x,y-1,z]
                     +sxz[x,y,z]-sxz[x,y,z-1])
        vx[x,y,z]=vx[x,y,z]+dv*ts
    
    #Vy cases
    if matBCs[x,y,z] == 0: 
        dvyConst=2*(1/gs)*(1/(matDensity[x,y,z]+matDensity[x,y+1,z]))
        dv=dvyConst* ( sxy[x,y,z]-sxy[x-1,y,z]
                      +syy[x,y+1,z]-syy[x,y,z]
                      +syz[x,y,z]-syz[x,y,z-1])
        vy[x,y,z]=vy[x,y,z]+dv*ts
    #y = 0
    elif (matBCs[x,y,z] ==2 or matBCs[x-1,y,z] == 2 or matBCs[x,y,z-1] == 2):
        pass  #requires elements out of the direction
    elif matBCs[x,y+1,z] == 2:
        vy[x,y,z] += 2 * ts/gs * 1/(2 * matDensity[x,y,z]) * ((-2)*syy[x,y,z])
    elif matBCs[x,y-1,z] == 2:
        vy[x,y,z] += 2 * ts/gs * 1/(2 * matDensity[x,y,z]) * ((2)*syy[x,y+1,z])
    else:
        dvyConst=2*(1/gs)*(1/(matDensity[x,y,z]+matDensity[x,y+1,z]))
        dv=dvyConst* ( sxy[x,y,z]-sxy[x-1,y,z]
                      +syy[x,y+1,z]-syy[x,y,z]
                      +syz[x,y,z]-syz[x,y,z-1])
        vy[x,y,z]=vy[x,y,z]+dv*ts

    #Vz cases
    if matBCs[x,y,z] ==0:
        dvzConst=2*(1/gs)*(1/(matDensity[x,y,z]+matDensity[x,y,z+1]))
        dv=dvzConst*( sxz[x,y,z]-sxz[x-1,y,z]
                     +syz[x,y,z]-syz[x,y-1,z]
                     +szz[x,y,z+1]-szz[x,y,z])
        vz[x,y,z]=vz[x,y,z]+dv*ts
    #z at 0
    elif (matBCs[x,y,z] == 2 or matBCs[x-1,y,z] == 2 or matBCs[x,y-1,z]==2):
        pass
    elif matBCs[x,y,z+1] == 2:
        vz[x,y,z] += 2 * ts/gs *(1/(2 * matDensity[x,y,z])) * ((-2)*szz[x,y,z])
    elif matBCs[x,y,z-1] == 2:
        vz[x,y,z] += 2 * ts/gs *(1/(2 * matDensity[x,y,z])) * ((2)*szz[x,y,z+1])
    else:
        dvzConst=2*(1/gs)*(1/(matDensity[x,y,z]+matDensity[x,y,z+1]))
        dv=dvzConst*( sxz[x,y,z]-sxz[x-1,y,z]
                     +syz[x,y,z]-syz[x,y-1,z]
                     +szz[x,y,z+1]-szz[x,y,z])
        vz[x,y,z]=vz[x,y,z]+dv*ts

def setSimSpaceBC99(matBCs):
    
    matBCs[0,:,:]=99
    matBCs[xmax,:,:]=99
    matBCs[:,0,:]=99
    matBCs[:,ymax,:]=99
    matBCs[:,:,0]=99
    matBCs[:,:,zmax]=99
    
    return matBCs

def setSimSpaceBCs(matBCs):
    #Second Dimmension boundaries /y
    matBCs[:,0,:]=2
    matBCs[:,1,:]=1
    matBCs[:,ymax,:]=2
    matBCs[:,ymax-1,:]=2
    matBCs[:,ymax-2,:]=1

    #Third Dimmension Boundaries /z
    matBCs[:,:,0]=2
    matBCs[:,2:ymax-1,1]=1
    matBCs[:,:,zmax]=2
    matBCs[:,:,zmax-1]=2
    matBCs[:,2:ymax-1,zmax-2]=1
    
    #First Dimmension Boundaries /x
    #   handled different if this is going to be calculated by node
    #   others c does it different, but they split between nodes before calculating
    #   here we calculate the whole set and then parse
    matBCs[0,:,:]=2
    matBCs[1,2:ymax-1,1:zmax-1]=1
    matBCs[xmax,:,:]=2
    matBCs[xmax-1,:,:]=2
    matBCs[xmax-2,1:ymax-1,1:zmax-1]=1
    
    return matBCs

def setRailBCs(matBCs):
    #top of foot
    matBCs[:,1:gridStartWeb,gridEndFoot]=1
    matBCs[:,gridEndWeb:ymax-1,gridEndFoot]=1

    #Sides Web
    matBCs[:,gridStartWeb,gridEndFoot:gridStartHead] = 1
    matBCs[:,gridEndWeb,gridEndFoot:gridStartHead] =1

    #bottom Head
    matBCs[:,gridStartHeadWidth:gridStartWeb+1,gridStartHead] = 1
    matBCs[:,gridEndWeb:gridEndHeadWidth,gridStartHead] = 1

    #Sides HEad
    matBCs[:,gridStartHeadWidth,gridStartHead:zmax-1] = 1
    matBCs[:,gridEndHeadWidth,gridStartHead:zmax-1] = 1

    #air beside Web
    matBCs[:,1:gridStartWeb,gridEndFoot+1:gridStartHead] = 2
    matBCs[:,gridEndWeb+1:ymax,gridEndFoot+1:gridStartHead] = 2

    #air beside head
    matBCs[:,1:gridStartHeadWidth,gridStartHead:zmax] = 2
    matBCs[:,gridEndHeadWidth+1:ymax,gridStartHead:zmax] = 2

    
    return matBCs

def MakeFlaw(matBCs):
    MidPoint = int(gl1/2)
    StartTrans = int(gl1/5)*2
    EndTrans = int(gl1/5)*3
    
    TransToEnd = gl1-EndTrans
    MidTransToEnd = int(TransToEnd/2)+EndTrans
    QuarterTrans = int((EndTrans-StartTrans)/4)
    
    StartFlawX = MidTransToEnd - QuarterTrans
    EndFlawX = MidTransToEnd + QuarterTrans
    
    StartFlawY = MidPoint - QuarterTrans
    EndFlawY = MidPoint + QuarterTrans
    
    VertFlaw = int(gh1/8)
    VertStart = zmax - VertFlaw
    
    #main hole
    matBCs[StartFlawX:EndFlawX,StartFlawY:EndFlawY,VertStart:] = 2
    
    #edges
    matBCs[StartFlawX:EndFlawX,StartFlawY:EndFlawY,VertStart-1] = 1
    matBCs[StartFlawX-1,StartFlawY-1:EndFlawY+1,VertStart:zmax-2]=1
    matBCs[EndFlawX+1,StartFlawY-1:EndFlawY+1,VertStart:zmax-2]=1
    matBCs[StartFlawX-1:EndFlawX+1,StartFlawY-1,VertStart:zmax-2]=1
    matBCs[StartFlawX-1:EndFlawX+1,EndFlawY+1,VertStart:zmax-2]=1
    
    return matBCs
    

#matBCs = setSimSpaceBC99(matBCs)
matBCall = setSimSpaceBCs(matBCall)
    
if RailShape:
    #matDensity,matLambda,matMu = setAirCut(matDensity,matLambda,matMu)
    matBCall = setRailBCs(matBCall)
    #matBCs = addTies(matBCs,Ties)

#Add Flaw
if Flaw:
    matBCall = MakeFlaw(matBCall)

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
    print("line 369: glb inputx, local inputx id, local inputx:",inputx,inputid,inputlocx)


szzConst=2*ts/(gs*rho1)

amp=10000
decayRate= 0
sinConst=ts*amp/rho1

sinInputSignal=sinConst*np.sin(2*np.pi*frequency*timeVec)*np.exp(-decayRate*timeVec)
#sinInputSignal[int(.1*Tsteps+1):] = 0

#Create End Damping to make asorbing boundary
Absorber = np.ones((gl1,gw1,gh1))
StepAbsorption = 0.5
AbsorptionRange = 101

if Absorbing:
    for x in range(AbsorptionRange):
        Absorber[x:AbsorptionRange+1,:,:] *= StepAbsorption


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

matDensity = np.zeros((npx,gw1,gh1))
matLambda = np.zeros((npx,gw1,gh1))
matMu = np.zeros((npx,gw1,gh1))
matBCs = np.zeros((npx,gw1,gh1))
signalloc = np.zeros((npx,gw1,gh1))
mpiAbsorber=np.zeros((npx,gw1,gh1))

mpi_comm.Scatterv([matDensityAll,split,offset,MPI.DOUBLE], matDensity)
mpi_comm.Scatterv([matLambdaAll,split,offset,MPI.DOUBLE], matLambda)
mpi_comm.Scatterv([matMuAll,split,offset,MPI.DOUBLE], matMu)
mpi_comm.Scatterv([matBCall,split,offset,MPI.DOUBLE], matBCs)
mpi_comm.Scatterv([signalLocation[:,:,:],split,offset,MPI.DOUBLE], signalloc)
mpi_comm.Scatterv([Absorber,split,offset,MPI.DOUBLE], mpiAbsorber)


matDensity=distBox(matDensity,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
matLambda=distBox(matLambda,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
matMu=distBox(matMu,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
matBCs=distBox(matBCs,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
signalloc=distBox(signalloc,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
mpiAbsorber=distBox(mpiAbsorber,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        

#Now slab has local versions with ghosts of matProps
if (myid == 0) :
    print('split matprops to globs, scratttered parameters to processors, line 449')
    
    ## All signals at 4 nodes from end
    FromEnd = int(.05 / gs) # 5 cm from end, about where laser recording is done
    
    # Top of rail
    FSignalLocX= gl1-FromEnd
    FSignalLocY=int(gw1/2)
    FSignalLocZ=gh1-2

    ## End halfway up head
    BSignalLocX=gl1-FromEnd
    BSignalLocY=int(gw1/2)
    BSignalLocZ=gh1-int((gh1-gridStartHead)/2)

    ## left of head
    USignalLocX=gl1-FromEnd
    USignalLocY=gridStartHeadWidth
    USignalLocZ=gh1-int((gh1-gridStartHead)/2)

    ## right of head
    DSignalLocX=gl1-FromEnd
    DSignalLocY=gridEndHeadWidth
    DSignalLocZ=gh1-int((gh1-gridStartHead)/2)

    ## right of web
    RSignalLocX=gl1-FromEnd
    RSignalLocY=gridStartWeb
    RSignalLocZ=int(gh1/2)

    ## Left of Web
    LSignalLocX=gl1-FromEnd
    LSignalLocY=gridEndWeb
    LSignalLocZ=int(gh1/2)

    ## End of Web
    MSignalLocX=gl1-2
    MSignalLocY=int(gw1/2)
    MSignalLocZ=int(gh1/2)


    #signal locations going to be a quarter of the way in the middle from the 
    # Front, Back, Up side, Down side, Right, Left, and Middle Middle Middle
    FSignal=np.zeros((Tsteps,3))
    BSignal=np.zeros((Tsteps,3))
    USignal=np.zeros((Tsteps,3))
    DSignal=np.zeros((Tsteps,3))
    RSignal=np.zeros((Tsteps,3))
    LSignal=np.zeros((Tsteps,3))
    MSignal=np.zeros((Tsteps,3))
    
    # switch to 1000 from Tsteps for size, and then added a saving every 1000 to fit in memory
    Movements = np.zeros((gl1,gw1,gh1,SaveSize))
    MovementsX = np.zeros((gl1,gw1,gh1,SaveSize))
    MovementsY = np.zeros((gl1,gw1,gh1,SaveSize))
    MovementsZ = np.zeros((gl1,gw1,gh1,SaveSize))
    DisX = np.zeros((gl1,gw1,gh1))
    DisY = np.zeros((gl1,gw1,gh1))
    DisZ = np.zeros((gl1,gw1,gh1))
    
    
    Parameters = {"AirCut" : AirCut,
                  "RailShape":  RailShape,
                  "Flaw" : Flaw,
                  "AbsorptionUsed" : Absorbing,
                  "SavingToCSV" : CSVs,
                  "SavingToPickle" : Pickles,
                  "MakingImages" : Images,
                  "Length" : length1,
                  "Width" : width1,
                  "Height" : height1,
                  "Wheel1Distance" : Wheel1Distance,
                  "SaveFolder" : imFolder,
                  "RunTitle" : runName,
                  "TiesIncluded" : Ties,
                  "GridDesignFrequency" : frequency,
                  "InputSignalFrequency" : Signalfrequency,
                  "TimeStepFreqeuency" : 1 / ts,
                  "SimulationCycleLength" : cycles,
                  "ForcingFuctionNumber" : FFunction,
                  "PerWheelForce" : WheelLoad,
                  "PoisonsRatio" : pRatio1,
                  "YoungsModulous" : yModulus1,
                  "MaterialDensity" : rho1,
                  "LongitudinalWaveSpeed" : cl1,
                  "TransverseWaveSpeeed" : ct1,
                  "TimeStepsSimLength" : Tsteps,
                  "GridLengthNodes" : gl1,
                  "GridWidthNodes" : gw1,
                  "GridHeightNodes" : gh1,
                  "LargestXnode" : xmax,
                  "LargestYnode" : ymax,
                  "LargestZnode" : zmax,
                  "GridStep" : gs,
                  "TimeStep" : ts,
                  "RunTime" : runtime,
                  "SaveEveryXStep" : SaveSize,
                  "HeightStartHeadNode" : gridStartHead,
                  "WidthStartWebNode" : gridStartWeb,
                  "WidthEndWebNode" : gridEndWeb,
                  "HeightEndFootNode" : gridEndFoot,
                  "WidthStartHeadNode" : gridStartHeadWidth,
                  "WidthEndHeadNode" : gridEndHeadWidth,
                  "AbsorberLengthNodes" : AbsorptionRange,
                  "AbsorptionPerNode" : StepAbsorption
                 }
                  
    file=open(imFolder+'Parameters.p','wb')
    pickle.dump(Parameters,file)
    file.close()
    
    del Parameters
    
    MinMax = np.zeros((4,2))
    
    j=0
    
stime = time.time()


if (myid == 0 and CSVs == True):
    print('subs setup, line 1213.  About to start at ' + str(stime))
    writeFile = open(imFolder + 'LaserPoints.csv','a')
    writeFile.write("Time,topX, topY, topZ, endX, endY, endZ, rHeadX, rHeadY, rHeadZ, lHeadX, lHeadY, lHeadZ, rWebX, rWebY, rWebZ, lWebX, lWebY, lWebZ\n")
    AnimationData =  open(imFolder + 'Anima.csv','a')
    AnimationData.write("time,x,y,z,Energy,DisX,DisY,DisZ\n")

if Images == True:
    MidMatrix = np.zeros((gl1,Tsteps))
    MidMatrixX = np.zeros((gl1,Tsteps))
    MidMatrixY = np.zeros((gl1,Tsteps))
    MidMatrixZ = np.zeros((gl1,Tsteps))

    Movements = np.zeros((gl1,gw1,int(Tsteps/10)))
    
DisX = np.zeros((gl1,gw1,gh1))
DisY = np.zeros((gl1,gw1,gh1))
DisZ = np.zeros((gl1,gw1,gh1))



for t in range(0,Tsteps):
    if FFunction ==2:
        vz += signalloc * sinInputSignal[t]
    elif FFunction ==3:
        vx += signalloc * sinInputSignal[t]
    elif FFunction ==4:
        vx += signalloc * sinInputSignal[t]
    elif FFunction ==5:
        vz -= signalloc * sinInputSignal[t]
    elif FFunction ==6:
        vz -= signalloc * sinInputSignal[t]


    for x in range(1,npx+1):
        for y in range(0,ymax):
            for z in range(0,zmax):
                JBSU(x,y,z)


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


    for x in range(1,npx+1):
        for y in range(0,ymax):
            for z in range(0,zmax):
                JBUV(x,y,z)
    vx*=mpiAbsorber
    vy*=mpiAbsorber
    vz*=mpiAbsorber
        
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
        if CSVs == True:
            USignal[t]=[vxg[USignalLocX,USignalLocY,USignalLocZ],
                        vyg[USignalLocX,USignalLocY,USignalLocZ],
                        vzg[USignalLocX,USignalLocY,USignalLocZ]]
            DSignal[t]=[vxg[DSignalLocX,DSignalLocY,DSignalLocZ],
                        vyg[DSignalLocX,DSignalLocY,DSignalLocZ],
                        vzg[DSignalLocX,DSignalLocY,DSignalLocZ]]
            RSignal[t]=[vxg[RSignalLocX,RSignalLocY,RSignalLocZ],
                        vyg[RSignalLocX,RSignalLocY,RSignalLocZ],
                        vzg[RSignalLocX,RSignalLocY,RSignalLocZ]]
            LSignal[t]=[vxg[LSignalLocX,LSignalLocY,LSignalLocZ],
                        vyg[LSignalLocX,LSignalLocY,LSignalLocZ],
                        vzg[LSignalLocX,LSignalLocY,LSignalLocZ]]
            MSignal[t]=[vxg[MSignalLocX,MSignalLocY,MSignalLocZ],
                        vyg[MSignalLocX,MSignalLocY,MSignalLocZ],
                        vzg[MSignalLocX,MSignalLocY,MSignalLocZ]]
            FSignal[t]=[vxg[FSignalLocX,FSignalLocY,FSignalLocZ],
                        vyg[FSignalLocX,FSignalLocY,FSignalLocZ],
                        vzg[FSignalLocX,FSignalLocY,FSignalLocZ]]
            BSignal[t]=[vxg[BSignalLocX,BSignalLocY,BSignalLocZ],
                        vyg[BSignalLocX,BSignalLocY,BSignalLocZ],
                        vzg[BSignalLocX,BSignalLocY,BSignalLocZ]]

            DisX += vxg[:,:,:] * ts
            DisY += vyg[:,:,:] * ts
            DisZ += vzg[:,:,:] * ts
            DisX += vxg[:,:,:] * ts
            DisY += vyg[:,:,:] * ts
            DisZ += vzg[:,:,:] * ts
            Movements[:,:,:,t%SaveSize] = np.sqrt(DisX**2 + DisY**2 + DisZ**2)  # ad the modulous 1000 for big model and multi save
            MovementsX[:,:,:,t%SaveSize] = DisX
            MovementsY[:,:,:,t%SaveSize] = DisY
            MovementsZ[:,:,:,t%SaveSize] = DisZ

            Min = np.min(Movements[:,:,:,t%SaveSize])
            Max = np.max(Movements[:,:,:,t%SaveSize])
            if Min < MinMax[0,0]:
                MinMax[0,0] = Min
            if Min < MinMax[0,1]:
                MinMax[0,0] = Min

            Min = np.min(MovementsX[:,:,:,t%SaveSize])
            Max = np.max(MovementsX[:,:,:,t%SaveSize])
            if Min < MinMax[1,0]:
                MinMax[1,0] = Min
            if Min < MinMax[1,1]:
                MinMax[1,1] = Min

            Min = np.min(MovementsY[:,:,:,t%SaveSize])
            Max = np.max(MovementsY[:,:,:,t%SaveSize])
            if Min < MinMax[2,0]:
                MinMax[2,0] = Min
            if Min < MinMax[2,1]:
                MinMax[2,1] = Min

            Min = np.min(MovementsZ[:,:,:,t%SaveSize])
            Max = np.max(MovementsZ[:,:,:,t%SaveSize])
            if Min < MinMax[3,0]:
                MinMax[3,0] = Min
            if Min < MinMax[3,1]:
                MinMax[3,1] = Min

        if Images == True:
            MidMatrixX[:,t] = vxg[:,inputy,inputz]
            MidMatrixY[:,t] = vyg[:,inputy,inputz]
            MidMatrixZ[:,t] = vzg[:,inputy,inputz]

            MidMatrixX.append(vxg[:,inputy,inputz])
            MidMatrixY.append(vxg[:,inputy,inputz])
            MidMatrixZ.append(vxg[:,inputy,inputz])

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
           
            
            fig=plt.figure()
            plt.contourf(np.transpose(np.sqrt(DisX[:,gridStartHeadWidth+1,gridStartHead:]**2 + 
                                              DisY[:,gridStartHeadWidth+1,gridStartHead:]**2 + 
                                              DisZ[:,gridStartHeadWidth+1,gridStartHead:]**2)), cmap='seismic')
            plt.savefig(imFolder+'HeadStart/HS'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)   
           
            
            fig=plt.figure()
            plt.contourf(np.transpose(np.sqrt(DisX[:,gridEndHeadWidth-1,gridStartHead:]**2 + 
                                              DisY[:,gridEndHeadWidth-1,gridStartHead:]**2 + 
                                              DisZ[:,gridEndHeadWidth-1,gridStartHead:]**2)), cmap='seismic')
            plt.savefig(imFolder+'HeadEnd/HE'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)   
           
            
            fig=plt.figure()
            plt.contourf(np.transpose(np.sqrt(DisX[:,gridStartWeb+1,gridEndFoot:gridStartHead]**2 + 
                                              DisY[:,gridStartWeb+1,gridEndFoot:gridStartHead]**2 + 
                                              DisZ[:,gridStartWeb+1,gridEndFoot:gridStartHead]**2)), cmap='seismic')
            plt.savefig(imFolder+'WebStart/WS'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)   
           
            
            fig=plt.figure()
            plt.contourf(np.transpose(np.sqrt(DisX[:,gridEndWeb-1,gridEndFoot:gridStartHead]**2 + 
                                              DisY[:,gridEndWeb-1,gridEndFoot:gridStartHead]**2 + 
                                              DisZ[:,gridEndWeb-1,gridEndFoot:gridStartHead]**2)), cmap='seismic')
            plt.savefig(imFolder+'WebEnd/WE'+str(t).zfill(5)+'.png')
            plt.close(fig)   
        
        if CSVs == True:
        
            writeFile.write(str(t)+","
                             +str(FSignal[-1][0])+","+str(FSignal[-1][1])+","+str(FSignal[-1][2])+","
                             +str(BSignal[-1][0])+","+str(BSignal[-1][1])+","+str(BSignal[-1][2])+","
                             +str(USignal[-1][0])+","+str(USignal[-1][1])+","+str(USignal[-1][2])+","
                             +str(DSignal[-1][0])+","+str(DSignal[-1][1])+","+str(DSignal[-1][2])+","
                             +str(RSignal[-1][0])+","+str(RSignal[-1][1])+","+str(RSignal[-1][2])+","
                             +str(LSignal[-1][0])+","+str(LSignal[-1][1])+","+str(LSignal[-1][2])+","
                             +str(MSignal[-1][0])+","+str(MSignal[-1][1])+","+str(MSignal[-1][2])
                             +"\n"
                            )
        
            for x in range(np.shape(DisX)[0]):
                for y in range(np.shape(DisX)[1]):
                    for z in range(np.shape(DisX)[2]):
                        if matBCall[x,y,z] == 1:
                            AnimationData.write(str(t)+","+str(x)+","+str(y)+","+str(z)+","
                                                +str(np.sqrt(DisX[x,y,z]**2+DisY[x,y,z]**2+DisZ[x,y,z]**2))+","
                                                +str(DisX[x,y,z])+","+str(DisY[x,y,z])+","+str(DisZ[x,y,z])
                                                +"\n")
        
            
        if Pickles == True and t%SaveSize == 0:
            file=open(imFolder+'MovementsR2MM'+str(j).zfill(3)+'.p','wb')
            pickle.dump([Movements, MovementsX, MovementsY, MovementsZ],file)
            file.close()
            j+=1

            

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
        print(t,'/',Tsteps-1,'checksums vx, sxx:',ckvs,ckss, int((time.time()-stime)/60.0*100)/100)
    sys.stdout.flush()

if Images == True:
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



if myid ==0:
    #print(MidMatrix)
    
        
    if Pickels == True:
        file=open(imFolder+'MovementsR2MM'+str(j).zfill(3)+'.p','wb')
        pickle.dump([Movements[:,:,:,:Tsteps%SaveSize], MovementsX[:,:,:,:Tsteps%SaveSize], 
                     MovementsY[:,:,:,:Tsteps%SaveSize], MovementsZ[:,:,:,:Tsteps%SaveSize]],file) 
        #don't save the trailing bad data from the last
        file.close()

        file=open(imFolder+'MinMax.p','wb')
        pickle.dump(MinMax, file) 
        file.close()
    
    if CSVs == True:
        del MovementsX
        del MovementsY
        del MovementsZ

        del DisX
        del DisY
        del DisZ
    
        #MidMatrixX = np.matrix(MidMatrixX)
        #MidMatrixY = np.matrix(MidMatrixY)
        #MidMatrixZ = np.matrix(MidMatrixZ)
        writeFile.close()
    
        # %%
    if Images == True:
        import multiprocessing
        from joblib import Parallel, delayed
        num_jobs=30

        # %%
        def EnergyFig(t, xStart, xEnd, yStart, yEnd, zStart, zEnd, Folder, v, figH):

            fig = plt.figure(figsize=(6,figH), dpi=300)

            Image = Movements[xStart:xEnd, yStart:yEnd, zStart:zEnd,t].T

            x,y,z = np.shape(Image)
            if x ==1:
                I2 = np.squeeze(Image, axis=(0,))
            if y ==1:
                I2 = np.squeeze(Image, axis=(1,))
            if z ==1:
                I2 = np.squeeze(Image, axis=(2,))

            plt.contourf(I2, v, cmap=plt.cm.jet)
            plt.savefig(imFolder+Folder+'/Energy'+str(t).zfill(5)+'.png')
            plt.close(fig)

        # %%
        def AnimationBook(xStart, xEnd, yStart, yEnd, zStart, zEnd, Folder):

            if xStart - xEnd == 0.0:
                figH = 6 * (yEnd - yStart) / (zEnd-zStart)
                xEnd+=1
            elif yStart - yEnd == 0.0:
                figH = 6 * (zEnd - zStart) / (xEnd-xStart)
                yEnd+=1
            elif zStart - zEnd == 0.0:
                figH = 6 * (yEnd - yStart) / (xEnd-xStart)
                zEnd+=1
            else:
                figH = 0

            if figH==0:
                print("Error, no Dimmension is a plane",yStart-yEnd)
            else:
                EMin = np.min(Movements[xStart:xEnd, yStart:yEnd, zStart:zEnd,:])
                EMax = np.max(Movements[xStart:xEnd, yStart:yEnd, zStart:zEnd,:])
                v = np.linspace(EMin, EMax, 30, endpoint=True)[0:20]

                print('About to make frames for ',Folder)
                temp = Parallel(n_jobs=num_jobs)(delayed(EnergyFig)(t, xStart, xEnd, yStart, yEnd, zStart, zEnd, Folder,v,figH) for t in range(Tsteps))




        # %%
        AnimationBook(0,xmax,gridStartWeb+1,gridStartWeb+1,gridEndFoot,gridStartHead,"WebStart")
        # %%
        plt.close('all')
        AnimationBook(0,xmax,gridEndWeb-1,gridEndWeb-1,gridEndFoot,gridStartHead,"WebEnd")
        # %%
        plt.close('all')
        AnimationBook(0,xmax,gridEndHeadWidth-1,gridEndHeadWidth-1,gridStartHead,zmax,"HeadEnd")
        # %%
        plt.close('all')
        AnimationBook(0,xmax,gridStartHeadWidth+1,gridStartHeadWidth+1,gridStartHead,zmax,"HeadStart")
        # %%
        plt.close('all')
        AnimationBook(0,xmax,gridStartHeadWidth,gridEndHeadWidth,zmax-3,zmax-3,"TopSurface")
        # %%
        plt.close('all')
        AnimationBook(0,xmax,gridStartHeadWidth,gridStartHeadWidth,gridStartHead,zmax,"LeftSurface")
        # %%
        plt.close('all')
        AnimationBook(0,xmax,gridStartWeb,gridStartWeb,gridEndFoot,gridStartHead,"RightSurface")
        # %%
        plt.close('all')
        
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
       
    if CSVs == True:
        SaveData = np.concatenate((FSignal,BSignal), axis = 1)
        SaveData = np.concatenate((SaveData,USignal), axis =1)
        SaveData = np.concatenate((SaveData,DSignal), axis =1)
        SaveData = np.concatenate((SaveData,RSignal), axis =1)
        SaveData = np.concatenate((SaveData,LSignal), axis =1)
        SaveData = np.concatenate((SaveData,MSignal), axis =1)

        np.savetxt(runName+'LaserPoints.csv',SaveData,delimiter=", ")
    