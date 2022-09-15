#Standard Header used on the projects

#first the major packages used for math and graphing
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler

from joblib import Parallel, delayed

#Custome graph format style sheet
#plt.style.use('Prospectus.mplstyle')

#If being run by a seperate file, use the seperate file's graph format and saving paramaeters
#otherwise set what is needed
if not 'Saving' in locals():
    Saving = True
if not 'Titles' in locals():
    Titles = True
if not 'Ledgends' in locals():
    Ledgends = True
if not 'FFormat' in locals():
    FFormat = '.eps'
if not 'location' in locals():
    #save location.  First one is for running on home PC, second for running on the work laptop.  May need to make a global change
    location = 'E:\\Documents\\Dan\\Code\\FigsAndPlots\\FigsAndPlotsDocument\\Figures\\'
    #location = 'C:\\Users\\dhendrickson\\Documents\\github\\FigsAndPlots\\FigsAndPlotsDocument\\Figures\\'

# Task Specific includes:

#import scipy.special as sp
import math
import matplotlib.animation as animation
import time
from numpy import inf
# Choose which EFIT_Class to use:
# import EFIT_Class as EFIT
# import EFIT_Class_StressLargerVelocity as EFIT
# import EFIT_Class_OrignalEqualGrid as EFIT
# import EFIT_Class_Parallel_EqualGrid as EFIT
# import EFIT_Class_VelocityLargerStress as EFIT
import EFIT_Class_OrignalEqualGrid_axisFlipper as EFIT
#import visvis as vv

# set Constants:
PoissonRatio = 0.3
YoungModulus = 20 * (10**9)
mu = 80 * (10**9)         #First Lame Parameter
lmbda = 2 * mu * PoissonRatio / (1 - 2 * PoissonRatio)     #second Lame Parameter
rho = 7800       #density kg/m^3

#Calculate speed of longitudinal and transverse waves
cl = np.sqrt((lmbda + 2* mu)/rho)
ct = np.sqrt(mu/rho)

#Choose ferquency to be used for excitment
frequency = 40000

#calculate wave length
omegal = cl / frequency
omegat = ct / frequency

#check max step size
dtmax = 1/ (max(cl,ct) * np.sqrt(3/(min(omegal,omegat)/10)))


# about 1foot (0.3m) of just the web of 175lbs rail 
# BeamLength = 0.3
# BeamHeight = 0.0762
# BeamWidth = 0.0381
BeamLength = 0.02
BeamHeight = 0.02
BeamWidth = 0.02

#Run for 6 Cycles:
runtime = 4.0 / frequency 

#Set time step and grid step to be 10 steps per frequency and ten steps per wavelength respectively
ts = 1 / frequency / 200  #  5e-08  #time step
gs = min(omegal, omegat) / 400    #grid step

Tsteps = int(math.ceil(runtime / ts)) + 1       #total Time STeps

gl = int(math.ceil(BeamLength / gs))         #number of grid points
gh = int(math.ceil(BeamHeight / gs)) 
gw = int(math.ceil(BeamWidth / gs)) 

print(dtmax,ts)

def makeAnimation(CenterZResults, title='Z'):
    y = np.linspace(0, BeamHeight, np.shape(CenterZResults[0][0])[1])
    x = np.linspace(0, BeamLength, np.shape(CenterZResults[0][0])[0])
    x,y = np.meshgrid(x,y)

    fig = plt.figure(figsize=(6.0,BeamHeight/BeamLength*6.0), dpi=100)
    ax = plt.axes(xlim=(0, BeamLength), ylim=(0, BeamHeight))  
    plt.ylabel(r'height')
    plt.xlabel(r'length')

    # animation function
    def animate(i): 
        z = np.matrix(CenterZResults[i][0][:,:]).T
        plt.title(str(i) + ' : ' + "{:.3e}".format(CenterZResults[i][1]))
        cont = plt.contourf(x, y, z, levels=5, cmap='gray') #,vmin=-100, vmax=100)
        #time.sleep(1)
        return cont  

    anim = animation.FuncAnimation(fig, animate, frames=np.shape(CenterZResults)[0])

    anim.save('animation'+title+'.gif')
    
# Inputs for forcing Function
Power = 10          #not sure on unit yet, just something for now
EmitterSize = 0.01  # meters, so 1 CM
Dimmension = 2      # 2 is in the z axis 
Direction = 1       # 1 is on top going down
CornerCut = 3       # how much of the corner is taken off of the square emitter
StartStep = 2
EndStep = Tsteps

#Run main function for time:
def loooop(looooooping): 
    
    lba=[]
    for i in range(3):
        lba.append(int('{0:08b}'.format(looooooping)[0:3+i][-1]))

    #[sx,sy,sz,vx,vy,vz] = lba #[1,1,1,1,1,1]
    [sx,sy,sz] = lba
    [vx,vy,vz]=[0,0,0]

    #Initialize EFIT Model
    Rail = EFIT.EFIT(gl, gw, gh, ts, gs, sx, sy, sz)

    #Set Material Properties consitant througout
    Rail.Gp[0,:,:,:] = rho  #constant Density
    Rail.Gp[1,:,:,:] = lmbda #Constant first Lamee parameter 
    Rail.Gp[2,:,:,:] = mu  #constant second Lamee parameter

    Rail.Gv[:,:,:,:]=0
    Rail.Gs[:,:,:,:,:]=0

    # CenterXResults = []
    # CenterYResults = []
    CenterZResults = []
    # AllVelocities=[]
    # VelocitiesX=[]
    # VelocitiesY=[]
    # VelocitiesZ=[]


    t=0
    CenterZResults.append((np.matrix(Rail.VelocityCut(1,2)),t))


    for i in range(Tsteps - 1):
        t = (i+1) * ts
    
        #Update Stresses at next half step:
        Rail.StepStresses()

        #update Velocity:
        Rail.StepVelocities()
        
        if i >= StartStep and i <= EndStep: 
            #Rail.ForcingFunctionImpulse(Power,EmitterSize,Dimmension,Direction, CornerCut)
            Rail.ForcingFunctionWave(t, frequency/2, 1,0.005,2,2)
            #ForcingFunctionWave(self, t, Hz = 40000, EP=100.0, size=0.04, Odim = 2, Dir=1):
            #
        #else:
        #    Rail.ForcingFunctionWave(t, frequency, 0)
        
        
        #print(str(i+1) + ' of ' + str(Tsteps-1) +' time steps. time is: '+ "{:.3e}".format(t)) #str(t))

        #  Save off each time step the currrent state for anmication later
        CenterZResults.append((np.matrix(Rail.VelocityCut(1,2)),t))
        # AllVelocities.append( Rail.VelocitySave())
        # VelocitiesX.append(Rail.VelocitySave(0))
        # VelocitiesY.append(Rail.VelocitySave(1))
        # VelocitiesZ.append(Rail.VelocitySave(2))
        
        # Store results mid process for latter animating
        if i % 10 == 9:
            print(str(i+1) + ' of ' + str(Tsteps-1) +' time steps. time is: '+ "{:.3e}".format(t)+' on loop '+str(looooooping)) #str(t))
        
        # Other data save out options at different time steps

    makeAnimation(CenterZResults, 'tries '+str(looooooping)) #+str(vx)+str(vy)+str(vz))

    return Rail


Results = Parallel(n_jobs=62)(delayed(loooop)(i) for i in range(8))

