import numpy as np
import matplotlib.pyplot as plt
import scipy as sp
import math
import time

# set Constants

#MATERIAL 1 (bronze)
pRatio1 = 0.34                                  #poission's ratio in 
yModulus1 = 103 * (10**9)                           #youngs modulus in pascals
rho1 = 8200                                        #density in kg/m^3

mu1 = yModulus1/(2*(1+pRatio1))                    #second Lame Parameter
lmbda1 = 2 * mu1 * pRatio1 / (1 - 2 * pRatio1)     #first Lame Parameter

#MATERIAL 1 (concrete)
#pRatio1 = 0.15                                    #poission's ratio in 
#yModulus1 = 24 * (10**9)                           #youngs modulus in pascals
#rho1 = 2400                                        #density in kg/m^3
#
#mu1 = yModulus1/(2*(1+pRatio1))                    #second Lame Parameter
#lmbda1 = 2 * mu1 * pRatio1 / (1 - 2 * pRatio1)     #first Lame Parameter


#Calculate speed of longitudinal and transverse waves in material 1
cl1 = np.sqrt((lmbda1 + 2* mu1)/rho1)
ct1 = np.sqrt(mu1/rho1)

print('material 1 wave speeds:' ,cl1,ct1)

#Choose ferquency to be used for excitment
frequency = 40000

#calculate wave lengths for material 1
omegaL1 = cl1 / frequency
omegaT1 = ct1 / frequency




#MATERIAL 2 (iron)
pRatio2=0.3
yModulus2= 200*(10**9)
rho2 = 7800       
mu2 = yModulus2/(2*(1+pRatio2))                    
lmbda2 = 2 * mu2 * pRatio2 / (1 - 2 * pRatio2)     

#Calculate speed of longitudinal and transverse waves in material 1
cl2= np.sqrt((lmbda2 + 2* mu2)/rho2)
ct2 = np.sqrt(mu2/rho2)

print('material 2 wave speeds:' ,cl2,ct2)

#calculate wavelengths in material 2
omegaL2 = cl2 / frequency
omegaT2 = ct2 / frequency


#dimensions of materials in meters
#the dimensions of material 1 should be greater than material 2
#diensions of material 1 

length1 = 0.60
width1 = 0.60
height1 = 0.150

#dimensions for material 2

length2= 0.1
width2= 0.1
height2= 0.15

#Run for 6 Cycles:
runtime = 8 / frequency 

#Set time step and grid step to be 10 steps per frequency and ten steps per wavelength respectively
#ts = 1 / frequency / 10    #time step
gs = (min(omegaL1, omegaT1,omegaL2,omegaT2) / 18)    #grid step
ts = gs/((max(cl1,ct1,cl2,ct2))*(np.sqrt(3)))*0.95 #time step

Tsteps = int(math.ceil(runtime / ts)) + 1       #total Time Steps

#number of grid points
gl1 = int(math.ceil(length1 / gs)) *2       #length 
gw1 = int(math.ceil(width1 / gs)) *2        #width
gh1 = int(math.ceil(height1 / gs)) *2       #height


#print(runtime, ts, gs, Tsteps, gl, gh)

print('runtime (s), time step size (s), total # of time steps:', runtime, ts, Tsteps)
print('grid step size, # of length pts, # of height pts, # of width pts:', gs, gl1, gw1,gh1)

#tensor to store material properties for each point
#0 index is density
#1 index is first Lame Parmaeter
#2 index is second lame parameter

matProps=np.zeros((3,gl1,gw1,gh1))
matProps[0,:,:,:]=rho1
matProps[1,:,:,:]=lmbda1
matProps[2,:,:,:]=mu1

#defining where the 2nd material lies in the grid
#for now we'll simulate a square rod of material 2 in the middle of a block 
#made up of material 1

gl2 = int(math.ceil(length2 / gs))       
gw2 = int(math.ceil(width2 / gs))         
gh2 = int(math.ceil(height2 / gs))

for x in range(gl1):
    for y in range(gw1):
        if x < y:
            matProps[0,x,y,:]=rho2
            matProps[1,x,y,:]=lmbda2
            matProps[2,x,y,:]=mu2

print('Total # of grid pts, # of material 1 pts, # of material 2 pts:', [gl1*gw1*gh1, gl1*gw1*gh1-8*gl2*gw2*gh2, 8*gl2*gw2*gh2])
print('material 2 indices:')
print('x:',[int(gl1/2)-gl2, int(gl1/2)+gl2])
print('y:',[int(gw1/2)-gw2, int(gw1/2)+gw2])
print('z:', [0, gh2*2])

#3D scatter plot of grid pts to show geometry 
#of the 2 materials

x=np.linspace(0,length1,gl1)
y=np.linspace(0,width1,gw1)
z=np.linspace(0,height1,gh1)

X, Y, Z = np.meshgrid(x, y, z)

U=matProps[0,:,:,:]

plt.figure()
ax=plt.axes(projection='3d')
fig=ax.scatter3D(X,Y,Z,c=U, alpha=0.02, marker='.', cmap='cool')
plt.savefig('efit_3d_Shape.png')

#define sine-exponential wave excitation

timeVec=np.linspace(0,runtime,Tsteps)

#radius
r=3
inputx=2
inputy=int(gw1/2)
inputz=int(gh1/2)


szzConst=2*ts/(gs*rho1)

amp=10000
decayRate= 50000
sinConst=ts*amp/rho1

sinInputSignal=sinConst*np.sin(2*np.pi*frequency*timeVec)*np.exp(-decayRate*timeVec)

plt.plot(timeVec,sinInputSignal)

plt.ioff()

#boundary values 
xmax=gl1-1
ymax=gw1-1
zmax=gh1-1

#initialize fields
vx=np.zeros((gl1,gw1,gh1))
vy=np.zeros((gl1,gw1,gh1))
vz=np.zeros((gl1,gw1,gh1))

sxx=np.zeros((gl1,gw1,gh1))
syy=np.zeros((gl1,gw1,gh1))
szz=np.zeros((gl1,gw1,gh1))
sxy=np.zeros((gl1,gw1,gh1))
sxz=np.zeros((gl1,gw1,gh1))
syz=np.zeros((gl1,gw1,gh1))

#record the signal at a specified location
signalLocx=int(gl1/2)
signalLocy=int(gw1/2)
signalLocz=int(gh1/2)

vxSignal=np.zeros(Tsteps)
vySignal=np.zeros(Tsteps)
vzSignal=np.zeros(Tsteps)

stime = time.time()

for t in range(0,Tsteps):
    
       
    #sin-exponential input
    vz[inputx,inputy,inputz]=vz[inputx,inputy,inputz]-szzConst*szz[inputx,inputy,inputz]+sinInputSignal[t]
    
    
    #update Stresses
    for x in range(gl1):
        for y in range(gw1):
            for z in range(gh1):
                
                #Calculate constants for stress equations
                norm1=(1/gs)*(matProps[1,x,y,z]+2*matProps[2,x,y,z])
                norm2=(1/gs)*(matProps[1,x,y,z])
                
                
                if x!=0 and x!=xmax and y!=0 and y!=ymax:
                    shearDenomxy=(1/matProps[2,x,y,z])+(1/matProps[2,x+1,y,z])+(1/matProps[2,x,y+1,z])+(1/matProps[2,x+1,y+1,z])
                    shearxy=4*(1/gs)*(1/shearDenomxy)
                    
                if x!=0 and x!=xmax and z!=0 and z!=zmax:
                    
                    shearDenomxz=(1/matProps[2,x,y,z])+(1/matProps[2,x+1,y,z])+(1/matProps[2,x,y,z+1])+(1/matProps[2,x+1,y,z+1])
                    shearxz=4*(1/gs)*(1/shearDenomxz)
                
                if y!=0 and y!=ymax and z!=0 and z!=zmax:
                    
                    shearDenomyz=(1/matProps[2,x,y,z])+(1/matProps[2,x,y+1,z])+(1/matProps[2,x,y,z+1])+(1/matProps[2,x,y+1,z+1])
                    shearyz=4*(1/gs)*(1/shearDenomyz)
                
                #FACES
                if x!=0 and x!=xmax and y!=0 and y!=ymax and z==0:
                    
                    
                    ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
                    sxx[x,y,z]=sxx[x,y,z]+ds*ts
                    
                    ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
                    syy[x,y,z]=syy[x,y,z]+ds*ts
                    
                    szz[x,y,z]=-szz[x,y,z+1]
                    
                    ds=shearxy*(vx[x,y+1,z]-vx[x,y,z]+vy[x+1,y,z]-vy[x,y,z])
                    sxy[x,y,z]=sxy[x,y,z]+ds*ts
                    del shearxy
                    
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                elif x==0 and y!=0 and y!=ymax and z!=0 and z!=zmax:
                    
                    
                    sxx[x,y,x]=-sxx[x+1,y,z]
                    
                    ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
                    syy[x,y,z]=syy[x,y,z]+ds*ts
                    
                    ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
                    szz[x,y,z]=szz[x,y,z]+ds*ts
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    
                    ds=shearyz*(vy[x,y,z+1]-vy[x,y,z]+vz[x,y+1,z]-vz[x,y,z])
                    syz[x,y,z]=syz[x,y,z]+ds*ts
                    del shearyz
                    
                elif x==xmax and y!=0 and y!=ymax and z!=0 and z!=zmax:
                    sxx[x,y,z]=-sxx[x-1,y,z]
                    
                    ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
                    syy[x,y,z]=syy[x,y,z]+ds*ts
                    
                    ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
                    szz[x,y,z]=szz[x,y,z]+ds*ts
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    
                    ds=shearyz*(vy[x,y,z+1]-vy[x,y,z]+vz[x,y+1,z]-vz[x,y,z])
                    syz[x,y,z]=syz[x,y,z]+ds*ts
                    del shearyz
                    
                elif x!=0 and x!=xmax and y==0 and z!=0 and z!=zmax:
                    
                    ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
                    sxx[x,y,z]=sxx[x,y,z]+ds*ts
                    
                    syy[x,y,z]=-syy[x,y+1,z]
                    
                    ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
                    szz[x,y,z]=szz[x,y,z]+ds*ts
                    
                    sxy[x,y,z]=0
                    
                    ds=shearxz*(vx[x,y,z+1]-vx[x,y,z]+vz[x+1,y,z]-vz[x,y,z])
                    sxz[x,y,z]=sxz[x,y,z]+ds*ts
                    del shearxz
                    
                    syz[x,y,z]=0
                    
                elif x!=0 and x!=xmax and y==ymax and z!=0 and z!=zmax:
                    
                    ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
                    sxx[x,y,z]=sxx[x,y,z]+ds*ts
                    
                    syy[x,y,z]=-syy[x,y-1,z]
                    
                    ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
                    szz[x,y,z]=szz[x,y,z]+ds*ts
                    
                    sxy[x,y,z]=0
                    
                    ds=shearxz*(vx[x,y,z+1]-vx[x,y,z]+vz[x+1,y,z]-vz[x,y,z])
                    sxz[x,y,z]=sxz[x,y,z]+ds*ts
                    del shearxz
                    
                    syz[x,y,z]=0
                    
                elif x!=0 and x!=xmax and y!=0 and y!=ymax and z==zmax:
                    
                    ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
                    sxx[x,y,z]=sxx[x,y,z]+ds*ts
                    
                    ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
                    syy[x,y,z]=syy[x,y,z]+ds*ts
                    
                    szz[x,y,z]=-szz[x,y,z-1]
                    
                    ds=shearxy*(vx[x,y+1,z]-vx[x,y,z]+vy[x+1,y,z]-vy[x,y,z])
                    sxy[x,y,z]=sxy[x,y,z]+ds*ts
                    del shearxy
                    
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                    
                #EDGES
                #bottom edges
                elif x==0 and y!=0 and y!=ymax and z==0:
                    
                    sxx[x,y,z]=-sxx[x+1,y,z]
                    
                    ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
                    syy[x,y,z]=syy[x,y,z]+ds*ts
                    
                    szz[x,y,z]=-szz[x,y,z+1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                elif x==xmax and y!=0 and y!=ymax and z==0:
                    
                    sxx[x,y,z]=-sxx[x-1,y,z]
                    
                    ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
                    syy[x,y,z]=syy[x,y,z]+ds*ts
                    
                    szz[x,y,z]-szz[x,y,z+1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                    
                elif x!=0 and x!=xmax and y==0 and z==0:
                    
                    ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
                    sxx[x,y,z]=sxx[x,y,z]+ds*ts
                    
                    syy[x,y,z]=-syy[x,y+1,z]
                    
                    szz[x,y,z]=-szz[x,y,z+1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                    
                elif x!=0 and x!=xmax and y==ymax and z==0:
                    
                    ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
                    sxx[x,y,z]=sxx[x,y,z]+ds*ts
                    
                    syy[x,y,z]=-syy[x,y-1,z]
                    
                    szz[x,y,z]=-szz[x,y,z+1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                #side edges
                elif x==0 and y==0 and z!=0 and z!=zmax:
                    
                    sxx[x,y,z]=-sxx[x+1,y,z]
                    
                    syy[x,y,z]=-syy[x,y+1,z]
                    
                    ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
                    szz[x,y,z]=szz[x,y,z]+ds*ts
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                    
                elif x==xmax and y==0 and z!=0 and z!=zmax:
                    
                    sxx[x,y,z]=-sxx[x-1,y,z]
                    
                    syy[x,y,z]=-syy[x,y+1,z]
                    
                    ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
                    szz[x,y,z]=szz[x,y,z]+ds*ts
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                elif x==0 and y==ymax and z!=0 and z!=zmax:
                    
                    sxx[x,y,z]=-sxx[x+1,y,z]
                    
                    syy[x,y,z]=-syy[x,y-1,z]
                    
                    ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
                    szz[x,y,z]=szz[x,y,z]+ds*ts
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                    
                elif x==xmax and y==ymax and z!=0 and z!=zmax:
                    
                    sxx[x,y,z]=-sxx[x-1,y,z]
                    
                    syy[x,y,z]=-syy[x,y-1,z]
                    
                    ds=norm1*(vz[x,y,z]-vz[x,y,z-1])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vy[x,y,z]-vy[x,y-1,z])
                    szz[x,y,z]=szz[x,y,z]+ds*ts
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                    
                #top edges
                elif x==0 and y!=0 and y!=ymax and z==zmax:
                    
                    sxx[x,y,z]=-sxx[x+1,y,z]
                    
                    ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
                    syy[x,y,z]=syy[x,y,z]+ds*ts
                    
                    szz[x,y,z]=-szz[x,y,z-1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                elif x==xmax and y!=0 and y!=ymax and z==zmax:
                    
                    sxx[x,y,z]=-sxx[x-1,y,z]
                    
                    ds=norm1*(vy[x,y,z]-vy[x,y-1,z])+norm2*(vx[x,y,z]-vx[x-1,y,z]+vz[x,y,z]-vz[x,y,z-1])
                    syy[x,y,z]=syy[x,y,z]+ds*ts
                    
                    szz[x,y,z]=-szz[x,y,z-1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                elif x!=0 and x!=xmax and y==0 and z==zmax:
                    
                    
                    ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
                    sxx[x,y,z]=sxx[x,y,z]+ds*ts
                    
                    syy[x,y,z]=-syy[x,y+1,z]
                    
                    szz[x,y,z]=-szz[x,y,z-1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                    
                elif x!=0 and x!=xmax and y==ymax and z==zmax:
                    
                    
                    ds=norm1*(vx[x,y,z]-vx[x-1,y,z])+norm2*(vy[x,y,z]-vy[x,y-1,z]+vz[x,y,z]-vz[x,y,z-1])
                    sxx[x,y,z]=sxx[x,y,z]+ds*ts
                    
                    syy[x,y,z]=-syy[x,y-1,z]
                    
                    szz[x,y,z]=-szz[x,y,z-1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                    
                #CORNERS
                
                elif x==0 and y==0 and z==0:
                    
                    sxx[x,y,z]=-sxx[x+1,y,z]
                    
                    syy[x,y,z]=-syy[x,y+1,z]
                    
                    szz[x,y,z]=-szz[x,y,z+1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                elif x==0 and y==0 and z==zmax:
                    
                    sxx[x,y,z]=-sxx[x+1,y,z]
                    
                    syy[x,y,z]=-syy[x,y+1,z]
                    
                    szz[x,y,z]=-szz[x,y,z-1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                elif x==0 and y==ymax and z==0:
                    
                    sxx[x,y,z]=-sxx[x+1,y,z]
                    
                    syy[x,y,z]=-syy[x,y-1,z]
                    
                    szz[x,y,z]=-szz[x,y,z+1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                elif x==0 and y==ymax and z==zmax:
                    
                    sxx[x,y,z]=-sxx[x+1,y,z]
                    
                    syy[x,y,z]=-syy[x,y-1,z]
                    
                    szz[x,y,z]=-szz[x,y,z-1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                elif x==xmax and y==0 and z==0:
                    
                    sxx[x,y,z]=-sxx[x-1,y,z]
                    
                    syy[x,y,z]=-syy[x,y+1,z]
                    
                    szz[x,y,z]=-szz[x,y,z+1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                elif x==xmax and y==0 and z==zmax:
                    
                    sxx[x,y,z]=-sxx[x-1,y,z]
                    
                    syy[x,y,z]=-syy[x,y+1,z]
                    
                    szz[x,y,z]=-szz[x,y,z-1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                elif x==xmax and y==ymax and z==0:
                    
                    sxx[x,y,z]=-sxx[x-1,y,z]
                    
                    syy[x,y,z]=-syy[x,y-1,z]
                    
                    szz[x,y,z]=-szz[x,y,z+1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                elif x==xmax and y==ymax and z==zmax:
                    
                    sxx[x,y,z]=-sxx[x-1,y,z]
                    
                    syy[x,y,z]=-syy[x,y-1,z]
                    
                    szz[x,y,z]=-szz[x,y,z-1]
                    
                    sxy[x,y,z]=0
                    sxz[x,y,z]=0
                    syz[x,y,z]=0
                    
                #NORMAL UPDATE
                else:
                    
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

                #delete variables for updates
                del norm1, norm2
                
                
    
    #update Velocities
    for x in range(gl1):
        for y in range(gw1):
            for z in range(gh1):
                
                #calculate constants for velocity
                if x!=xmax:
                    dvxConst=2*(1/gs)*(1/(matProps[0,x,y,z]+matProps[0,x+1,y,z]))
                    vminx=(2*ts)/(matProps[0,x+1,y,z]*gs)

                if y!=ymax:
                    dvyConst=2*(1/gs)*(1/(matProps[0,x,y,z]+matProps[0,x,y+1,z]))
                    vminy=(2*ts)/(matProps[0,x,y+1,z]*gs)

                if z!=zmax:
                    dvzConst=2*(1/gs)*(1/(matProps[0,x,y,z]+matProps[0,x,y,z+1]))
                    vminz=(2*ts)/(matProps[0,x,y,z+1]*gs)
                

                vmax=(2*ts)/(matProps[0,x,y,z]*gs)

                
                
                #FACES
                if x!=0 and x!=xmax and y!=0 and y!=ymax and z==0:
                    dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
                    vx[x,y,z]=vx[x,y,z]+dv*ts
                    
                    dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
                    vy[x,y,z]=vy[x,y,z]+dv*ts
                    
                    vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]
                    
                    
                elif x==0 and y!=0 and y!=ymax and z!=0 and z!=zmax:
                    
                    vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
                    
                    dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
                    vy[x,y,z]=vy[x,y,z]+dv*ts
                    
                    dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
                    vz[x,y,z]=vz[x,y,z]+dv*ts

                    
                elif x==xmax and y!=0 and y!=ymax and z!=0 and z!=zmax:
                    
                    vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
                    
                    dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
                    vy[x,y,z]=vy[x,y,z]+dv*ts
                    
                    dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
                    vz[x,y,z]=vz[x,y,z]+dv*ts
                    
                elif x!=0 and x!=xmax and y==0 and z!=0 and z!=zmax:
                    
                    dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
                    vx[x,y,z]=vx[x,y,z]+dv*ts
                    
                    vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
                    
                    dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
                    vz[x,y,z]=vz[x,y,z]+dv*ts
                    
                   
                elif x!=0 and x!=xmax and y==ymax and z!=0 and z!=zmax:
                    
                    dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
                    vx[x,y,z]=vx[x,y,z]+dv*ts
                    
                    vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
                    
                    dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
                    vz[x,y,z]=vz[x,y,z]+dv*ts                     
                    
                    
                elif x!=0 and x!=xmax and y!=0 and y!=ymax and z==zmax:
                    
                    dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
                    vx[x,y,z]=vx[x,y,z]+dv*ts
                    
                    
                    dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
                    vy[x,y,z]=vy[x,y,z]+dv*ts
                    
                    vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]

                    
                    
                #EDGES
                #bottom edges
                elif x==0 and y!=0 and y!=ymax and z==0:
                    
                    vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
                   
                    dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
                    vy[x,y,z]=vy[x,y,z]+dv*ts
                    
                    vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]
                    
                elif x==xmax and y!=0 and y!=ymax and z==0:
                    
                    vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
                    
                    dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
                    vy[x,y,z]=vy[x,y,z]+dv*ts
                    
                    vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]
                    
                elif x!=0 and x!=xmax and y==0 and z==0:
                    
                    dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
                    vx[x,y,z]=vx[x,y,z]+dv*ts
                    
                    vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
                    
                    vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]
                    
                elif x!=0 and x!=xmax and y==ymax and z==0:
                    
                    dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
                    vx[x,y,z]=vx[x,y,z]+dv*ts
                    
                    vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
                    
                    vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]
                    
                #side edges
                elif x==0 and y==0 and z!=0 and z!=zmax:
                    
                    vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
                    
                    vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
                    
                    dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
                    vz[x,y,z]=vz[x,y,z]+dv*ts

                elif x==xmax and y==0 and z!=0 and z!=zmax:
                    
                    vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
                    
                    vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
                    
                    dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
                    vz[x,y,z]=vz[x,y,z]+dv*ts
                    
                    
                elif x==0 and y==ymax and z!=0 and z!=zmax:
                    
                    vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
                    
                    vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
                    
                    dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
                    vz[x,y,z]=vz[x,y,z]+dv*ts
                    
                    
                    
                elif x==xmax and y==ymax and z!=0 and z!=zmax:
                    
                    vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
                    
                    vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
                    
                    dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
                    vz[x,y,z]=vz[x,y,z]+dv*ts
                    
                                        
                #top edges
                elif x==0 and y!=0 and y!=ymax and z==zmax:
                    
                    vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
                    
                    dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
                    vy[x,y,z]=vy[x,y,z]+dv*ts
                    
                    vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]
                    
                    
                elif x==xmax and y!=0 and y!=ymax and z==zmax:
                    
                    vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
                    
                    dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
                    vy[x,y,z]=vy[x,y,z]+dv*ts
                    
                    vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]
                    
                elif x!=0 and x!=xmax and y==0 and z==zmax:
                    
                    dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
                    vx[x,y,z]=vx[x,y,z]+dv*ts
                    
                    vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
                    
                    vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]
                    
                    
                elif x!=0 and x!=xmax and y==ymax and z==zmax:
                    
                    dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
                    vx[x,y,z]=vx[x,y,z]+dv*ts
                    
                    vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
                    
                    vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]
                    
                 
                #CORNERS
                elif x==0 and y==0 and z==0:
                    
                    vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
                    vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
                    vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]
                    
                elif x==0 and y==0 and z==zmax:
                    
                    vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
                    vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
                    vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]
                    
                elif x==0 and y==ymax and z==0:
                    
                    vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
                    vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
                    vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]
                    
                elif x==0 and y==ymax and z==zmax:
                    
                    vx[x,y,z]=vx[x,y,z]+vminx*sxx[x+1,y,z]
                    vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
                    vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]
                    
                elif x==xmax and y==0 and z==0:
                    
                    vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
                    vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
                    vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]
                    
                elif x==xmax and y==0 and z==zmax:
                    
                    vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
                    vy[x,y,z]=vy[x,y,z]+vminy*syy[x,y+1,z]
                    vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]
                    
                elif x==xmax and y==ymax and z==0:
                    
                    vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
                    vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
                    vz[x,y,z]=vz[x,y,z]+vminz*szz[x,y,z+1]
                    
                elif x==xmax and y==ymax and z==zmax:
                    
                    vx[x,y,z]=vx[x,y,z]-vmax*sxx[x,y,z]
                    vy[x,y,z]=vy[x,y,z]-vmax*syy[x,y,z]
                    vz[x,y,z]=vz[x,y,z]-vmax*szz[x,y,z]
                    
                else:
                    
                    dv=dvxConst*(sxx[x+1,y,z]-sxx[x,y,z]+sxy[x,y,z]-sxy[x,y-1,z]+sxz[x,y,z]-sxz[x,y,z-1])
                    vx[x,y,z]=vx[x,y,z]+dv*ts
    
                    dv=dvyConst*(sxy[x,y,z]-sxy[x-1,y,z]+syy[x,y+1,z]-syy[x,y,z]+syz[x,y,z]-syz[x,y,z-1])
                    vy[x,y,z]=vy[x,y,z]+dv*ts
    
                    dv=dvzConst*(sxz[x,y,z]-sxz[x-1,y,z]+syz[x,y,z]-syz[x,y-1,z]+szz[x,y,z+1]-szz[x,y,z])
                    vz[x,y,z]=vz[x,y,z]+dv*ts
            
                    del dvxConst, dvyConst, dvzConst
    
                
        
    
    #record signals
    vxSignal[t]=vx[signalLocx,signalLocy,signalLocz]
    vySignal[t]=vy[signalLocx,signalLocy,signalLocz]
    vzSignal[t]=vz[signalLocx,signalLocy,signalLocz]
    
    #save vx cut figure
    if t%10==0:
        fig=plt.figure()
        plt.contourf(np.transpose(vz[:,:,int(gh1/2)]), cmap='seismic')
        plt.title(str(t)+'/'+str(Tsteps-1)+' checksums vx, sxx:'+str(np.sum(np.absolute(vx)))+' '+str(np.sum(np.absolute(sxx)))+' time: '+str(time.time()-stime))
        plt.savefig('image/wvlts/zmidCut'+str(t)+'.png')
        plt.close(fig)
    
    
    print(t,'/',Tsteps-1,'checksums vx, sxx:',np.sum(np.absolute(vx)),np.sum(np.absolute(sxx)), time.time()-stime)

fig = plt.figure()
plt.contourf(np.transpose(vz[:,:,15]), cmap='seismic')
plt.savefig('result at 15.png')
plt.close(fig)


fig = plt.figure()
plt.plot(vxSignal)
plt.savefig('Xsignal Result.png')
plt.close(fig)

fig = plt.figure()
plt.plot(vySignal)
plt.savefig('Ysignal Result.png')
plt.close(fig)

fig = plt.figure()
plt.plot(vzSignal)
plt.savefig('Zsignal Result.png')
plt.close(fig)
