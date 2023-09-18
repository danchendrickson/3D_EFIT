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
RailShape = True
Flaw = False

#Dimmesnsion of simulation space in meters
length1 = 25
width1 = 0.1524 # 0.1524
height1 = 0.1524

#Image Folder
imFolder = '/sciclone/scr10/dchendrickson01/EFIT/'
runName = 'DoubleRub'

#is the rail supported by 0, 1 or 2 ties
Ties = 0

#Choose ferquency to be used for excitment
frequency = 49720.0  #brute forced this number to be where simulation frequency 
#                     is 2,000,000, alowing for %10 to equal laser 200k same as actual
#            74574  is 3,000,000 hz running, and sample rate %15 is 200k same as actual, if we need more dense
Signalfrequency = 16300

cycles = 10000

figDPI = 600

#Forcing Function Location and type
# 1 for dropped wheel on top
# 2 for rubbing flange on side
# 3 for plane wave 

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
    imFolder += 'DoubleTest/'

    
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

#Run for 3 seconds, what the laser can store:
runtime = 3 #cycles / frequency 

Tsteps = int(math.ceil(runtime / ts)) + 1       #total Time Steps

#number of grid points
gl1 = int(math.ceil(length1 / gs)) +1       #length 
gw1 = int(math.ceil(width1 / gs)) +1       #width
gh1 = int(math.ceil(height1 / gs)) +1       #height

if (myid == 0) :
    print('runtime, gs, ts, fr, gl, gw, gh, Tsteps, imFolder : ', runtime, gs, ts, 1/ts, gl1, gw1, gh1, Tsteps, imFolder)

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
     
    signalLocation[14:20,gridStartHeadWidth:gridStartHeadWidth+2,gridStartHead:zmax-2] = 1

    signalLocation[20,gridStartHeadWidth:gridStartHeadWidth+2,gridStartHead:zmax-2] = 0.5
    signalLocation[13,gridStartHeadWidth:gridStartHeadWidth+2,gridStartHead:zmax-2] = 0.5
    signalLocation[14:20,gridStartHeadWidth+2:gridStartHeadWidth+3,gridStartHead:zmax-2] = 0.5

    sep = int(1.360/gs)
    
    signalLocation[14+sep:20+sep,gridEndHeadWidth-2:gridEndHeadWidth,gridStartHead:zmax-2] = 1

    signalLocation[20+sep,gridEndHeadWidth-2:gridEndHeadWidth,gridStartHead:zmax-2] = 0.5
    signalLocation[13+sep,gridEndHeadWidth-2:gridEndHeadWidth,gridStartHead:zmax-2] = 0.5
    signalLocation[14+sep:20+sep,gridEndHeadWidth-3:gridEndHeadWidth-2,gridStartHead:zmax-2] = 0.5
    
    
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
    elif matBCs[z,y,z-1] == 2:
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

mpi_comm.Scatterv([matDensityAll,split,offset,MPI.DOUBLE], matDensity)
mpi_comm.Scatterv([matLambdaAll,split,offset,MPI.DOUBLE], matLambda)
mpi_comm.Scatterv([matMuAll,split,offset,MPI.DOUBLE], matMu)
mpi_comm.Scatterv([matBCall,split,offset,MPI.DOUBLE], matBCs)
mpi_comm.Scatterv([signalLocation[:,:,:],split,offset,MPI.DOUBLE], signalloc)


matDensity=distBox(matDensity,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
matLambda=distBox(matLambda,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
matMu=distBox(matMu,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
matBCs=distBox(matBCs,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        
signalloc=distBox(signalloc,myid,gl1,gw1,gh1,npx,nprocs,mpi_comm)        

#Now slab has local versions with ghosts of matProps
if (myid == 0) :
    print('split matprops to globs, scratttered parameters to processors, line 449')
    
    ## All signals at 4 nodes from end
    FromEnd = 50
    # Top of rail
    FSignalLocX= gl1-FromEnd
    FSignalLocY=int(gw1/2)
    FSignalLocZ=gh1-4

    ## End halfway up head
    BSignalLocX=gl1-FromEnd
    BSignalLocY=int(gw1/2)
    BSignalLocZ=int((gh1-gridStartHead)/2+gridStartHead)

    ## left of head
    USignalLocX=gl1-FromEnd
    USignalLocY=gridStartHeadWidth+1
    USignalLocZ=int((gh1-gridStartHead)/2+gridStartHead)

    ## right of head
    DSignalLocX=gl1-FromEnd
    DSignalLocY=gridEndHeadWidth-1
    DSignalLocZ=int((gh1-gridStartHead)/2+gridStartHead)

    ## right of web
    RSignalLocX=gl1-FromEnd
    RSignalLocY=gridStartWeb+1
    RSignalLocZ=int((gridStartHead-gridEndFoot)/2+gridEndFoot)

    ## Left of Web
    LSignalLocX=gl1-FromEnd
    LSignalLocY=gridEndWeb-1
    LSignalLocZ=int((gridStartHead-gridEndFoot)/2+gridEndFoot)

    ## End of Web
    MSignalLocX=gl1-FromEnd
    MSignalLocY=int(gw1/2)
    MSignalLocZ=int((gridStartHead-gridEndFoot)/2+gridEndFoot)


    #signal locations going to be a quarter of the way in the middle from the 
    # Front, Back, Up side, Down side, Right, Left, and Middle Middle Middle
    FSignal=np.zeros((Tsteps,3))
    BSignal=np.zeros((Tsteps,3))
    USignal=np.zeros((Tsteps,3))
    DSignal=np.zeros((Tsteps,3))
    RSignal=np.zeros((Tsteps,3))
    LSignal=np.zeros((Tsteps,3))
    MSignal=np.zeros((Tsteps,3))

stime = time.time()

if (myid == 0 ):
    print('subs setup, line 1213.  About to start at ' + str(stime))

#MidMatrix = np.zeros((gl1,Tsteps))
#MidMatrixX = np.zeros((gl1,Tsteps))
#MidMatrixY = np.zeros((gl1,Tsteps))
#MidMatrixZ = np.zeros((gl1,Tsteps))

#Movements = np.zeros((gl1,gw1,int(Tsteps/10)))
DisX = np.zeros((gl1,gw1,gh1))
DisY = np.zeros((gl1,gw1,gh1))
DisZ = np.zeros((gl1,gw1,gh1))

writeFile = open(imFolder + 'LaserPoints.csv','a')

for t in range(0,Tsteps):
    if FFunction ==2:
        vz += signalloc * sinInputSignal[t]
    elif FFunction ==3:
        vx += signalloc * sinInputSignal[t]
    elif FFunction ==4:
        vx += signalloc * sinInputSignal[t]
    elif FFunction ==5:
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
        USignal[t]=[vxg[USignalLocX,USignalLocY,USignalLocZ],vyg[USignalLocX,USignalLocY,USignalLocZ],vzg[USignalLocX,USignalLocY,USignalLocZ]]
        DSignal[t]=[vxg[DSignalLocX,DSignalLocY,DSignalLocZ],vyg[DSignalLocX,DSignalLocY,DSignalLocZ],vzg[DSignalLocX,DSignalLocY,DSignalLocZ]]
        RSignal[t]=[vxg[RSignalLocX,RSignalLocY,RSignalLocZ],vyg[RSignalLocX,RSignalLocY,RSignalLocZ],vzg[RSignalLocX,RSignalLocY,RSignalLocZ]]
        LSignal[t]=[vxg[LSignalLocX,LSignalLocY,LSignalLocZ],vyg[LSignalLocX,LSignalLocY,LSignalLocZ],vzg[LSignalLocX,LSignalLocY,LSignalLocZ]]
        MSignal[t]=[vxg[MSignalLocX,MSignalLocY,MSignalLocZ],vyg[MSignalLocX,MSignalLocY,MSignalLocZ],vzg[MSignalLocX,MSignalLocY,MSignalLocZ]]
        FSignal[t]=[vxg[FSignalLocX,FSignalLocY,FSignalLocZ],vyg[FSignalLocX,FSignalLocY,FSignalLocZ],vzg[FSignalLocX,FSignalLocY,FSignalLocZ]]
        BSignal[t]=[vxg[BSignalLocX,BSignalLocY,BSignalLocZ],vyg[BSignalLocX,BSignalLocY,BSignalLocZ],vzg[BSignalLocX,BSignalLocY,BSignalLocZ]]

        #MidMatrixX[:,t] = vxg[:,inputy,inputz]
        #MidMatrixY[:,t] = vyg[:,inputy,inputz]
        #MidMatrixZ[:,t] = vzg[:,inputy,inputz]

        #MidMatrixX.append(vxg[:,inputy,inputz])
        #MidMatrixY.append(vxg[:,inputy,inputz])
        #MidMatrixZ.append(vxg[:,inputy,inputz])
        
        DisX += vxg[:,:,:] * ts
        DisY += vyg[:,:,:] * ts
        DisZ += vzg[:,:,:] * ts

        if t%10==0:
            #Movements[:,:,int(t/5)-1] = np.sqrt(DisX[:,:,zmax-3]**2 + DisY[:,:,zmax-3]**2 + DisZ[:,:,zmax-3]**2)
            writeFile.write(str(t)+','
                             +str(FSignal[-1][0])+','+str(FSignal[-1][1])+','+str(FSignal[-1][2])+','
                             +str(BSignal[-1][0])+','+str(BSignal[-1][1])+','+str(BSignal[-1][2])+','
                             +str(USignal[-1][0])+','+str(USignal[-1][1])+','+str(USignal[-1][2])+','
                             +str(DSignal[-1][0])+','+str(DSignal[-1][1])+','+str(DSignal[-1][2])+','
                             +str(RSignal[-1][0])+','+str(RSignal[-1][1])+','+str(RSignal[-1][2])+','
                             +str(LSignal[-1][0])+','+str(LSignal[-1][1])+','+str(LSignal[-1][2])+','
                             +str(MSignal[-1][0])+','+str(MSignal[-1][1])+','+str(MSignal[-1][2])+','
                             +'/n'
                            )
            
            
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
    writeFile.close()
    
    '''if np.shape(MidMatrixX)[0] == Tsteps:
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
       
    SaveData = np.concatenate((FSignal,BSignal), axis = 1)
    SaveData = np.concatenate((SaveData,USignal), axis =1)
    SaveData = np.concatenate((SaveData,DSignal), axis =1)
    SaveData = np.concatenate((SaveData,RSignal), axis =1)
    SaveData = np.concatenate((SaveData,LSignal), axis =1)
    SaveData = np.concatenate((SaveData,MSignal), axis =1)
    
    np.savetxt(runName+'LaserPoints.csv',SaveData,delimiter=", ")
    '''