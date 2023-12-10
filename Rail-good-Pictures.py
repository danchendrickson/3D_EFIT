# %% [markdown]
# # 3D EFIT Rail Code in MPI
# 
# ## Combined Zane structure, Eric MPI, My system

# %%
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
import matplotlib.animation as animation

# %%
MPIComm = Union[MPI.Intracomm, MPI.Intercomm]
mpi_comm = MPI.COMM_WORLD
myid = mpi_comm.Get_rank()                                                         
mpi_size = mpi_comm.Get_size()        
nprocs=mpi_size
myid = 0

# for overlapping slabs:  
# # points per proc along z = npz = gh1/nproc (+2 to ghost boundaries)
# glob_index = loc_index-1 + npz*myid
# loc_index = glob_index - npz*myid + 1
# myid given glob_index = glob_index/npz = ghloc-2

# set Constants
AirCut = False
RailShape = True
figDPI = 600

#Dimmesnsion of simulation space in meters
length1 = 3
width1 = 0.1524 # 0.1524
height1 = 0.1524

#Image Folder
imFolder = '/sciclone/scr10/dchendrickson01/EFIT/'
runName = 'RailRubLongg'

#is the rail supported by 0, 1 or 2 ties
Ties = 0
Flaw = False

cycles = 65

#Choose ferquency to be used for excitment
frequency = 16300
frequency = 27139
frequency = 60000

#Run for 4 Cycles:
runtime = cycles / frequency 

#Forcing Function Location and type
# 1 for dropped wheel on top
# 2 for rubbing flange on side
# 3 for plane wave
FFunction = 2

#MATERIAL 1 ((steel))
pRatio1 = 0.29                                    #poission's ratio in 
yModulus1 = 200 * (10**9)                           #youngs modulus in pascals
rho1 = 7800                                        #density in kg/m^3

#Image Folder
if FFunction == 1:
    imFolder += 'TopHit/'
elif FFunction == 2:
    imFolder += 'RailRubLongg/'
elif FFunction ==3:
    imFolder += 'Cube/'
elif FFunction ==4:
    imFolder += 'RailDense/'

WheelLoad = 173000 #crane force in Neutons

#CALCULATED PARAMETERS FROM INPUTS

mu1 = yModulus1/(2*(1+pRatio1))                    #second Lame Parameter
lmbda1 = 2 * mu1 * pRatio1 / (1 - 2 * pRatio1)     #first Lame Parameter

#Calculate speed of longitudinal and transverse waves in material 1
cl1 = np.sqrt((lmbda1 + 2* mu1)/rho1)
ct1 = np.sqrt(mu1/rho1)

#calculate wave lengths for material 1
omegaL1 = cl1 / frequency
omegaT1 = ct1 / frequency

# %%
#Set time step and grid step to be 10 steps per frequency and ten steps per wavelength respectively
#ts = 1 / frequency / 10    #time step
gs = (min(omegaL1, omegaT1) /13)    #grid step
ts = gs/((max(cl1,ct1))*(np.sqrt(3)))*0.93 #time step

#frequency = 16300

#Run for 4 Cycles:
runtime = cycles / frequency 

Tsteps = int(math.ceil(runtime / ts)) + 1       #total Time Steps

#number of grid points
gl1 = int(math.ceil(length1 / gs)) +1       #length 
gw1 = int(math.ceil(width1 / gs)) +1       #width
gh1 = int(math.ceil(height1 / gs)) +1       #height

print(gs, ts, gl1, gw1, gh1, Tsteps)

# %%
# Keep these as the global values
xmax=gl1-1
ymax=gw1-1
zmax=gh1-1

#due to double think BCs, adding one, but keeping the maxes
#gw1+=1
#gh1+=1

# %%
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
matDensity=np.zeros((gl1,gw1,gh1))
matLambda=np.zeros((gl1,gw1,gh1))
matMu=np.zeros((gl1,gw1,gh1))
matBCs=np.zeros((gl1,gw1,gh1))
signalLocation=np.zeros((gl1,gw1,gh1))

matDensity[:,:,:]=rho1
matLambda[:,:,:]=lmbda1
matMu[:,:,:]=mu1
matBCs[:,:,:]=0

#Make the Signal Location grid
if FFunction == 1:
    pnodes = int(whlayer / 2)
    contactLength = int(0.001 / gs)  #1 cm contact patch

    signalLocation[0:contactLength,gridStartHeadWidth:gridEndHeadWidth, -3:] = 1
    
elif FFunction == 2:
    pnodes = int(whlayer / 4)
    contactLength = max(int(0.004 / gs),3)  #4 cm contact patch

    signalLocation[0:contactLength,gridStartHeadWidth:gridStartHeadWidth+3, gridStartHead:] = 1
    print(FFunction, contactLength, np.sum(signalLocation))

    ## Find the share of the force per node for FF1
    
elif FFunction ==3:
    '''signalLocation[3:4,2:ymax-1,2:zmax-1] = 1
    signalLocation[3:4,2:ymax-1,2:zmax-1] = 1
    signalLocation[2:3,2:ymax-1,2:zmax-1] = 0.5
    signalLocation[4:5,2:ymax-1,2:zmax-1] = 0.5
    '''
    signalLocation[13:15,:,:] = 1
    signalLocation[12:13,:,:] = 0.5
    signalLocation[15:16,:,:] = 0.5
    
    
elif FFunction == 4:

    signalLocation[int(gl1/2)-5:int(gl1/2)+5,int(gw1/2)-5:int(gw1/2)+5,zmax-2] = 1
    signalLocation[int(gl1/2)-5:int(gl1/2)+5,int(gw1/2)-5:int(gw1/2)+5,zmax-1] = 0.5
    signalLocation[int(gl1/2)-5:int(gl1/2)+5,int(gw1/2)-5:int(gw1/2)+5,zmax-3] = 0.5
    
    
specificWheelLoad = WheelLoad / np.sum(signalLocation)


if myid == 0:
    print('globs made, line 145')

# %%
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


gridStartHead = round(gh1 * relStartHeadThick)
gridStartWeb = round(gw1 * relStartWeb)
gridEndWeb = round(gw1 * relEndWeb)
gridEndFoot = round(gh1 * relFoot)
gridStartHeadWidth = round(gw1 * relStartHeadWidth)
gridEndHeadWidth = round(gw1  * relEndHeadWidth)



# %%
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


# %%
def setSimSpaceBC99(matBCs):
    
    matBCs[0,:,:]=99
    matBCs[xmax,:,:]=99
    matBCs[:,0,:]=99
    matBCs[:,ymax,:]=99
    matBCs[:,:,0]=99
    matBCs[:,:,zmax]=99
    
    return matBCs


# %%
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
    

# %%
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
    

# %%
def setRailBCs(matBCs):
    #set the boundary conditions in material props4
    # Set the simulations pace boundarys
    # top
    
    #top of footing
    matBCs[:,:gridStartWeb,gridEndFoot]=1
    matBCs[:,ymax-gridStartWeb:,gridEndFoot]=1

    #sides of web
    matBCs[:,gridStartWeb,gridStartHead:gridEndFoot]=1
    matBCs[:,gridEndWeb,gridStartHead:gridEndFoot]=1

    #bottom of head
    matBCs[:,gridStartHeadWidth:gridStartWeb,gridStartHead]=1
    matBCs[:,gridEndWeb:gridEndHeadWidth,gridStartHead]=1

    #sides of head
    matBCs[:,gridStartHeadWidth,gridStartHead:]=1
    matBCs[:,gridEndHeadWidth,gridStartHead:]=1

    #zone 1 of air, left of web
    matBCs[:,:gridStartWeb,gridEndFoot:gridStartHead]=2
    # zone 2 of air left of head
    matBCs[:,:gridStartHeadWidth,:gridStartHead]=2
    # zone 3 of air, right of web
    matBCs[:,gridEndWeb:,gridEndFoot:gridStartHead]=2
    # zone 4 of air, right of head
    matBCs[:,gridEndHeadWidth:,gridStartHead:]=2

    
    return matBCs
    

# %%


# %%
#matBCs = setSimSpaceBC99(matBCs)
matBCs = setSimSpaceBCs(matBCs)
    
if RailShape:
    #matDensity,matLambda,matMu = setAirCut(matDensity,matLambda,matMu)
    matBCs = setRailBCs(matBCs)
    #matBCs = addTies(matBCs,Ties)

#Add Flaw
#    136 x 136 x15
if Flaw:
    matBCs = MakeFlaw(matBCs)


# %%
print(matBCs[9,:,:])

# %%
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
#sinInputSignal-=0.0000001
#sinInputSignal[int(1/12*Tsteps):]=0

fig = plt.figure()
plt.plot(sinInputSignal)
plt.show()

# %%
# MPI EJW Section #4 changes 

#initialize fields
vx=np.zeros((npx,gw1,gh1))
vy=np.zeros((npx,gw1,gh1))
vz=np.zeros((npx,gw1,gh1))

sxx=np.zeros((npx,gw1,gh1))
syy=np.zeros((npx,gw1,gh1))
szz=np.zeros((npx,gw1,gh1))
sxy=np.zeros((npx,gw1,gh1))
sxz=np.zeros((npx,gw1,gh1))
syz=np.zeros((npx,gw1,gh1))

#record the signal at a specified location
### ADD map function for this
signalLocx=int(gl1/2)
signalLocy=int(gw1/2)
signalLocz=int(gh1/2)

vxSignal=np.zeros(Tsteps)
vySignal=np.zeros(Tsteps)
vzSignal=np.zeros(Tsteps)

# %%

#record the signal at a specified location
### ADD map function for this
#SAME AS INPUTZ?

#manually setting x for long where wave may not propaagate full lengthh in time

FSignalLocX= 5 #int(gl1/4)
BSignalLocX=25 # int(3*gl1/4)
USignalLocX=10 #int(gl1/4)
DSignalLocX=10 #int(gl1/4)
RSignalLocX=10 #int(gl1/4)
LSignalLocX=10 #int(gl1/4)
MSignalLocX=10 #int(gl1/4)

FSignalLocY=int(gw1/2)
BSignalLocY=int(gw1/2)
USignalLocY=int(gw1/2)
DSignalLocY=int(gw1/2)
RSignalLocY=int(gw1/4)
LSignalLocY=int(3*gw1/4)
MSignalLocY=int(gw1/2)

FSignalLocZ=int(gh1/2)
BSignalLocZ=int(gh1/2)
USignalLocZ=int(3*gh1/4)
DSignalLocZ=int(gh1/4)
RSignalLocZ=int(gh1/2)
LSignalLocZ=int(gh1/2)
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

# %%
signalloc = np.zeros((npx,gw1,gh1))
signalloc=signalLocation[:,:,:]

stime = time.time()

# %%
# asdfasdf
MidMatrixX = np.zeros((gl1,Tsteps))
MidMatrixY = np.zeros((gl1,Tsteps))
MidMatrixZ = np.zeros((gl1,Tsteps))

Movements = np.zeros((gl1,gw1,gh1,Tsteps))
DisX = np.zeros((gl1,gw1,gh1))
DisY = np.zeros((gl1,gw1,gh1))
DisZ = np.zeros((gl1,gw1,gh1))


# %%
for t in range(Tsteps):
     
    if FFunction == 2:
        vz += signalloc * sinInputSignal[t]
        print(np.sum(vz))
    
    if FFunction ==3:
        vx += signalloc * sinInputSignal[t]

    if FFunction == 4:
        vz += signalloc * sinInputSignal[t]

    for x in range(0,xmax):
        for y in range(0,ymax):
            for z in range(0,zmax):
                JBSU(x,y,z)

    
    # cut boundaries off of arrays
    '''Remove for Jupyter single processor
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
    '''
    
    #if the forcing function is a stress
    if FFunction == 1:
        szz -= signalloc * specificWheelLoad

    for x in range(0,xmax):
        for y in range(0,ymax):
            for z in range(0,zmax):
                JBUV(x,y,z)

    MidMatrixX[:,t] = vx[:,MSignalLocY,MSignalLocZ]
    MidMatrixY[:,t] = vy[:,MSignalLocY,MSignalLocZ]
    MidMatrixZ[:,t] = vz[:,MSignalLocY,MSignalLocZ]
    
    vxg=vx #[1:npx+1,:,:]        
    vyg=vy#[1:npx+1,:,:]        
    vzg=vz#[1:npx+1,:,:]        

    
    if (myid == 0 ) :
        USignal[t]=[vxg[USignalLocX,USignalLocY,USignalLocZ],vyg[USignalLocX,USignalLocY,USignalLocZ],vzg[USignalLocX,USignalLocY,USignalLocZ]]
        DSignal[t]=[vxg[DSignalLocX,DSignalLocY,DSignalLocZ],vyg[DSignalLocX,DSignalLocY,DSignalLocZ],vzg[DSignalLocX,DSignalLocY,DSignalLocZ]]
        RSignal[t]=[vxg[RSignalLocX,RSignalLocY,RSignalLocZ],vyg[RSignalLocX,RSignalLocY,RSignalLocZ],vzg[RSignalLocX,RSignalLocY,RSignalLocZ]]
        LSignal[t]=[vxg[LSignalLocX,LSignalLocY,LSignalLocZ],vyg[LSignalLocX,LSignalLocY,LSignalLocZ],vzg[LSignalLocX,LSignalLocY,LSignalLocZ]]
        MSignal[t]=[vxg[MSignalLocX,MSignalLocY,MSignalLocZ],vyg[MSignalLocX,MSignalLocY,MSignalLocZ],vzg[MSignalLocX,MSignalLocY,MSignalLocZ]]
        FSignal[t]=[vxg[FSignalLocX,FSignalLocY,FSignalLocZ],vyg[FSignalLocX,FSignalLocY,FSignalLocZ],vzg[FSignalLocX,FSignalLocY,FSignalLocZ]]
        BSignal[t]=[vxg[BSignalLocX,BSignalLocY,BSignalLocZ],vyg[BSignalLocX,BSignalLocY,BSignalLocZ],vzg[BSignalLocX,BSignalLocY,BSignalLocZ]]
        
        DisX += vx[:,:,:] * ts
        DisY += vy[:,:,:] * ts
        DisZ += vz[:,:,:] * ts
        Movements[:,:,:,t] = np.sqrt(DisX**2 + DisY**2 + DisZ**2)
        
        if t%5==0:
        
            fig=plt.figure()
            plt.contourf(np.transpose(vxg[:,:,int(gh1/2)]), cmap='seismic')
            plt.savefig(imFolder+'Mid/vyMidHeightShear'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)
            
            fig=plt.figure()
            plt.contourf(np.transpose(vxg[:,int(gw1/2),:]), cmap='seismic')
            plt.savefig(imFolder + 'Vert/vyVertShear'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)    
            
            fig=plt.figure()
            plt.contourf(np.transpose(vxg[int(gl1/2),:,:]), cmap='seismic')
            plt.savefig(imFolder + 'Head/vyMidLenShear'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)  
    
            fig=plt.figure()
            plt.contourf(np.transpose(vxg[:,:,int(gh1/4)]), cmap='seismic')
            plt.savefig(imFolder+'zplane25/vy25Shear'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)
            
            fig=plt.figure()
            plt.contourf(np.transpose(vxg[:,:,int(3*gh1/4)]), cmap='seismic')
            plt.savefig(imFolder+'zplane75/vy75Shear'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)
            '''
            fig=plt.figure()
            plt.contourf(np.transpose(Movements[:,:,zmax-1,t]), cmap='seismic')
            plt.savefig(imFolder+'TopSurface/TopSurface'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)
            
            fig=plt.figure()
            plt.contourf(np.transpose(Movements[3:,ymax,:,t]), cmap='seismic')
            plt.savefig(imFolder+'RightSurface/RightSurface'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)
            
            fig=plt.figure()
            plt.contourf(np.transpose(Movements[3:,0,:,t]), cmap='seismic')
            plt.savefig(imFolder+'LeftSurface/LeftSurface'+str(t).zfill(5)+'.png')
            # SideRub vs TopHit for which case
            plt.close(fig)
            '''

    # Collect vx, sxx checksum contributions for printing
    vxt=vx[1:npx+1,:,:]
    sxxt=sxx[1:npx+1,:,:]

    #ckvs=np.array(0.0,'d')
    #ckss=np.array(0.0,'d')
    
    ckv=np.sum(np.absolute(vxg))
    cks=np.sum(np.absolute(sxx))
    #mpi_comm.Reduce(ckv,ckvs,op=MPI.SUM,root=0)
    #mpi_comm.Reduce(cks,ckss,op=MPI.SUM,root=0)

    if (myid == 0 ):
        print(t,'/',Tsteps-1,'check vx, sxx:',ckv,cks, int(((time.time()-stime)/60.0)*100)/100)
    sys.stdout.flush()

# %%
vxDisplacement = [0]
vyDisplacement = [0]
vzDisplacement = [0]

Times = np.linspace(0, len(MSignal), num=len(MSignal)+1)
Times *= ts

for i in range(len(MSignal)):
    vxDisplacement.append(vxDisplacement[i-1]+MSignal[i][0] * ts)
    vyDisplacement.append(vyDisplacement[i-1]+MSignal[i][1] * ts)
    vzDisplacement.append(vzDisplacement[i-1]+MSignal[i][2] * ts)

plt.clf()
plt.title('Middle Node')
plt.plot(Times,vxDisplacement,label='x')
plt.plot(Times,vyDisplacement,label='y')
plt.plot(Times,vzDisplacement,label='z')
plt.legend()
plt.savefig(imFolder+runName+'DisplaceMid.png')
#plt.show()

# %%
vxDisplacement = [0]
vyDisplacement = [0]
vzDisplacement = [0]

for i in range(len(FSignal)):
    vxDisplacement.append(vxDisplacement[i-1]+FSignal[i][0] * ts)
    vyDisplacement.append(vyDisplacement[i-1]+FSignal[i][1] * ts)
    vzDisplacement.append(vzDisplacement[i-1]+FSignal[i][2] * ts)

plt.clf()
plt.title('Middle XY Plane Quarter into Rod')
plt.plot(Times,vxDisplacement,label='x', linewidth=2)
plt.plot(Times,vyDisplacement,label='y', linewidth=4)
plt.plot(Times,vzDisplacement,label='z', linewidth=2)
plt.legend()
plt.savefig(imFolder+runName+'DisplaceFront.png')
plt.show()

# %%
vxDisplacement = [0]
vyDisplacement = [0]
vzDisplacement = [0]

for i in range(len(BSignal)):
    vxDisplacement.append(vxDisplacement[i-1]+BSignal[i][0] * ts)
    vyDisplacement.append(vyDisplacement[i-1]+BSignal[i][1] * ts)
    vzDisplacement.append(vzDisplacement[i-1]+BSignal[i][2] * ts)

plt.clf()
plt.title('Back Node')
plt.plot(Times,vxDisplacement)
plt.plot(Times,vyDisplacement)
plt.plot(Times,vzDisplacement)
plt.savefig(imFolder+runName+'DisplaceBack.png')

# %%
vxDisplacement = [0]
vyDisplacement = [0]
vzDisplacement = [0]

for i in range(len(RSignal)):
    vxDisplacement.append(vxDisplacement[i-1]+RSignal[i][0] * ts)
    vyDisplacement.append(vyDisplacement[i-1]+RSignal[i][1] * ts)
    vzDisplacement.append(vzDisplacement[i-1]+RSignal[i][2] * ts)

plt.clf()
plt.title('Right Node')
plt.plot(Times,vxDisplacement)
plt.plot(Times,vyDisplacement)
plt.plot(Times,vzDisplacement)
plt.savefig(imFolder+runName+'DisplaceRight.png')

# %%
vxDisplacement = [0]
vyDisplacement = [0]
vzDisplacement = [0]

for i in range(len(LSignal)):
    vxDisplacement.append(vxDisplacement[i-1]+LSignal[i][0] * ts)
    vyDisplacement.append(vyDisplacement[i-1]+LSignal[i][1] * ts)
    vzDisplacement.append(vzDisplacement[i-1]+LSignal[i][2] * ts)

plt.clf()
plt.title('Left Node')
plt.plot(Times,vxDisplacement)
plt.plot(Times,vyDisplacement)
plt.plot(Times,vzDisplacement)
plt.savefig(imFolder+runName+'DisplaceLeft.png')

# %%
vxDisplacement = [0]
vyDisplacement = [0]
vzDisplacement = [0]

for i in range(len(USignal)):
    vxDisplacement.append(vxDisplacement[i-1]+USignal[i][0] * ts)
    vyDisplacement.append(vyDisplacement[i-1]+USignal[i][1] * ts)
    vzDisplacement.append(vzDisplacement[i-1]+USignal[i][2] * ts)

plt.clf()
plt.title('Up Node')
plt.plot(Times,vxDisplacement)
plt.plot(Times,vyDisplacement)
plt.plot(Times,vzDisplacement)
plt.savefig(imFolder+runName+'DisplaceUp.png')

# %%
vxDisplacement = [0]
vyDisplacement = [0]
vzDisplacement = [0]

for i in range(len(DSignal)):
    vxDisplacement.append(vxDisplacement[i-1]+DSignal[i][0] * ts)
    vyDisplacement.append(vyDisplacement[i-1]+DSignal[i][1] * ts)
    vzDisplacement.append(vzDisplacement[i-1]+DSignal[i][2] * ts)

plt.clf()
plt.title('Down Node')
plt.plot(Times,vxDisplacement)
plt.plot(Times,vyDisplacement)
plt.plot(Times,vzDisplacement)
plt.savefig(imFolder+runName+'DisplaceDown.png')

# %%
#Data = [MSignal,USignal,DSignal,LSignal,RSignal,FSignal,BSignal]
print(np.shape(MSignal), np.shape(np.asarray(MSignal)))


# %%
fig = plt.figure(dpi=600)
endnode=500
plt.clf()
plt.title('Down Node')
plt.plot(Times[:endnode],vxDisplacement[:endnode],label='x')
plt.plot(Times[:endnode],vyDisplacement[:endnode],label='y')
plt.plot(Times[:endnode],vzDisplacement[:endnode],label='z')
plt.legend()
plt.savefig(imFolder+runName+'DisplaceMid2.png')
plt.show()

# %%
Cases = ['PlaneWave']    #'SideRub','TopHit', 'Cube', 'CubeS'
Views = ['Mid', 'Vert','Head','zplane25','zplane75']

# %%
Displacement = np.zeros(np.shape(MSignal))

for i in range(np.shape(MSignal)[1]):
    for j in range(np.shape(MSignal)[0]):
        if j == 0:
            Displacement[j,i]=MSignal[j,i]*ts
        else:
            Displacement[j,i]=Displacement[j-1,i]+MSignal[j,i]*ts
    
fig = plt.figure(dpi=600)
endnode=150
plt.clf()
plt.title('Middle Node')
plt.plot(Times[:endnode],Displacement[:endnode,0],label='x', linewidth = 3)
plt.plot(Times[:endnode],Displacement[:endnode,1],label='y', linewidth = 6)
plt.plot(Times[:endnode],Displacement[:endnode,2],label='z', linewidth = 3)
plt.legend()
plt.savefig(imFolder+runName+'DisplaceMid2.png')
plt.show()

# %%
MidDisplaceX = np.zeros(np.shape(MidMatrixX))
MidDisplaceY = np.zeros(np.shape(MidMatrixY))
MidDisplaceZ = np.zeros(np.shape(MidMatrixZ))


# %%
for i in range(np.shape(MidMatrixX)[0]):
    for j in range(np.shape(MidMatrixX)[1]):
        if j == 0:
            MidDisplaceX[i,j]=MidMatrixX[i,j]*ts
            MidDisplaceY[i,j]=MidMatrixY[i,j]*ts
            MidDisplaceZ[i,j]=MidMatrixZ[i,j]*ts
        else:
            MidDisplaceX[i,j]=MidDisplaceX[i,j-1]+MidMatrixX[i,j]*ts
            MidDisplaceY[i,j]=MidDisplaceY[i,j-1]+MidMatrixY[i,j]*ts
            MidDisplaceZ[i,j]=MidDisplaceZ[i,j-1]+MidMatrixZ[i,j]*ts
            

# %%
pts = 8
rng = int(gl1/pts)-1

fig = plt.figure(dpi=600, figsize=(6,4))
for i in range(pts):
    plt.plot(MidMatrixX[i*rng,:],label=str(i*rng))
plt.title('Velocity')
plt.legend()
plt.savefig(imFolder+runName+'MidVelocities.png')
plt.show()

# %%
fig = plt.figure(dpi=600, figsize=(6,4))
for i in range(pts):
    plt.plot(MidDisplaceX[i*rng,:],label=str(i*rng))
plt.legend()
plt.title('Displacement')
plt.savefig(imFolder+runName+'MidDisplacementsX.png')
plt.show()

# %%
fig = plt.figure(dpi=600, figsize=(6,4))
for i in range(pts):
    plt.plot(MidDisplaceY[i*rng,:],label=str(i*rng))
plt.legend()
plt.title('Displacement')
plt.savefig(imFolder+runName+'MidDisplacementsY.png')
plt.show()

# %%
fig = plt.figure(dpi=600, figsize=(6,4))
for i in range(pts):
    plt.plot(MidDisplaceZ[i*rng,:],label=str(i*rng))
plt.legend()
plt.title('Displacement')
plt.savefig(imFolder+runName+'MidDisplacementsZ.png')
plt.show()

# %%
import multiprocessing
from joblib import Parallel, delayed
num_jobs=19

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

        temp = Parallel(n_jobs=num_jobs)(delayed(EnergyFig)(t, xStart, xEnd, yStart, yEnd, zStart, zEnd, Folder,v,figH) for t in range(Tsteps))
        
        


# %%
plt.close('all')
AnimationBook(0,xmax,gridEndWeb-1,gridEndWeb-1,gridEndFoot,gridStartHead,"WebEnd")


# %%
plt.close('all')
AnimationBook(0,xmax,gridEndHeadWidth-1,gridEndHeadWidth-1,gridStartHead,zmax,"HeadEnd")


# %%
plt.close('all')
AnimationBook(0,xmax,gridStartHeadWidth,gridEndHeadWidth,zmax-3,zmax-3,"TopSurface")
plt.close('all')
AnimationBook(0,xmax,gridStartHeadWidth,gridStartHeadWidth,gridStartHead,zmax,"LeftSurface")
plt.close('all')
AnimationBook(0,xmax,gridStartWeb,gridStartWeb,gridEndFoot,gridStartHead,"RightSurface")
plt.close('all')



