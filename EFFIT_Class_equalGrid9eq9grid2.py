#import limited
from zlib import Z_DEFAULT_COMPRESSION
import numpy as np
#import visvis as vv
class EFIT:
    ts = 0
    ds = 0

    def __init__(self, Grid1, Grid2, Grid3, tStep, dStep,StressPlus = 1,VelocityPlus=1):
        #Initialize with the size of the grid in 3 dimmensions in number of nodes.  The distance step, and the time step


        #Grid shape is the 3 grid dimsions
        self.GridShape = (Grid1,Grid2,Grid3)

        #total gird size
        self.GridPoints = (Grid1)*(Grid2)*(Grid3)

        #MAybe to be used later, allows for some axis to be non right hand rule
        self.StPl = StressPlus
        self.VePl = VelocityPlus
        self.axis1 = 1
        self.axis2 = 1
        self.axis3 = 1

        #Lame parameters and density set as 0, will need to be reset in code after object is created
        self.lmbda=0
        self.mu=0
        self.rho=0
        
        #define empty grid for the 3 directions of velocity
        self.Gv1 = np.zeros(self.GridPoints,dtype="double").reshape(*self.GridShape)
        self.Gv2 = np.zeros(self.GridPoints,dtype="double").reshape(*self.GridShape)
        self.Gv3 = np.zeros(self.GridPoints,dtype="double").reshape(*self.GridShape)
        
        #define empty grid for the 3 directions of stress on 3 dimmensions of faces, not that the two cross dimmensions
        #are not included as the Stress 1-2 is equal to the Stress 2-1
        self.Gs11 = np.zeros(self.GridPoints,dtype="double").reshape(*self.GridShape)
        self.Gs22 = np.zeros(self.GridPoints,dtype="double").reshape(*self.GridShape)
        self.Gs33 = np.zeros(self.GridPoints,dtype="double").reshape(*self.GridShape)
        self.Gs12 = np.zeros(self.GridPoints,dtype="double").reshape(*self.GridShape)
        self.Gs23 = np.zeros(self.GridPoints,dtype="double").reshape(*self.GridShape)
        self.Gs31 = np.zeros(self.GridPoints,dtype="double").reshape(*self.GridShape)
        #define empty grid for the 3 scalar material properties at each node point.  Can honly hold scalar properties
        #Assumed properties are density, Lame 1, Lame 2
        #self.Gp = np.zeros(3*self.GridPoints,dtype="double").reshape(*self.GridShapeP)
        
        self.Max1 = Grid1 - 1
        self.Max2 = Grid2 - 1
        self.Max3 = Grid3 - 1

        self.ds = dStep
        self.ts = tStep

        self.ids = 1.0 / dStep
    
    def DeltaStress(self,x1,x2,x3):
        #Gets the change in the stresses per time at a certain coordinate juncture 
        #
        # Inputs: x1,x2,x3 coordinates of the cube in question.  Last time is assumed
        #
        # Outputs: 6 dimmensions of stress.  
        
        #Calculated stresses based on 3.55

        Ds11 =  ((self.ids) *
            ((self.lmbda+2*self.mu)*(self.Gv1[x1,x2,x3]-self.Gv1[x1-1,x2,x3]) +
                        self.lmbda*((self.Gv2[x1,x2,x3]-self.Gv2[x1,x2-1,x3]) +
                                    (self.Gv3[x1,x2,x3]-self.Gv3[x1,x2,x3-1]))
                )
            )
        Ds22 =  ((self.ids) *
            ((self.lmbda+2*self.mu)*(self.Gv2[x1,x2,x3]-self.Gv2[x1,x2-1,x3]) +
                        self.lmbda*((self.Gv3[x1,x2,x3]-self.Gv3[x1,x2,x3-1]) +
                                    (self.Gv1[x1,x2,x3]-self.Gv1[x1-1,x2,x3]))
                )
            )
        Ds33 =  ((self.ids) *
            ((self.lmbda+2*self.mu)*(self.Gv3[x1,x2,x3]-self.Gv3[x1,x2,x3-1]) +
                        self.lmbda*((self.Gv1[x1,x2,x3]-self.Gv1[x1-1,x2,x3]) +
                                    (self.Gv2[x1,x2,x3]-self.Gv2[x1,x2-1,x3]))
                )
            )
        Ds12 = (self.ids) * ( self.lmbda ) * (
                (self.Gv1[x1,x2+1,x3]-self.Gv1[x1,x2,x3]) + 
                (self.Gv2[x1+1,x2,x3]-self.Gv2[x1,x2,x3])
            )
        Ds23 = (self.ids) * ( self.lmbda ) * (
                (self.Gv2[x1,x2,x3+1]-self.Gv2[x1,x2,x3]) + 
                (self.Gv3[x1,x2+1,x3]-self.Gv3[x1,x2,x3])
            )
        Ds31 = (self.ids) * ( self.lmbda ) * (
                (self.Gv3[x1+1,x2,x3]-self.Gv3[x1,x2,x3]) + 
                (self.Gv1[x1,x2,x3+1]-self.Gv1[x1,x2,x3])
            )

        return Ds11, Ds22, Ds33, Ds12, Ds23, Ds31

    def DeltaVelocity(self, x1,x2,x3):
        #Gets the change in the velocity per time at a certain coordinate juncture 
        #
        # Inputs: x,y,z coordinates of the cube in question.  Last time is assumed
        #
        # Outputs: 3 dimmensions of velocity.  
        
        #Calculated velocity based on 3.54

        Dv1 = ((self.ids ) / (self.rho) ) * (
                (self.Gs11[x1+1,x2,x3] - self.Gs11[x1,x2,x3]) +
                (self.Gs12[x1,x2,x3] - self.Gs12[x1,x2-1,x3]) + 
                (self.Gs31[x1,x2,x3] - self.Gs31[x1,x2,x3-1])
            )            
        Dv2 = ((self.ids ) / (self.rho) ) * (
                (self.Gs12[x1,x2,x3] - self.Gs12[x1-1,x2,x3]) +
                (self.Gs22[x1,x2+1,x3] - self.Gs22[x1,x2,x3]) + 
                (self.Gs23[x1,x2,x3] - self.Gs23[x1,x2,x3-1])
            )            
        Dv3 = ((self.ids ) / (self.rho) ) * (
                (self.Gs31[x1,x2,x3] - self.Gs31[x1-1,x2,x3]) +
                (self.Gs23[x1,x2,x3] - self.Gs23[x1,x2-1,x3]) + 
                (self.Gs33[x1,x2,x3+1] - self.Gs33[x1,x2,x3])
            )
            

        return Dv1, Dv2, Dv3
    
    def UpdateStresses(self, x,y,z):
        #Updates velocity based off of previous velocities and delta velocity times time
        #
        # Inputs: Coordinates of cube in question.  Assumed last time step
        #
        # Output: updated self.Gs matrix
        
        Ds11, Ds22, Ds33, Ds12, Ds23, Ds31 = self.DeltaStress(x,y,z)

        self.Gs11[x,y,z] +=  Ds11 * self.ts * self.StPl * self.axis1
        self.Gs22[x,y,z] +=  Ds22 * self.ts * self.StPl * self.axis2
        self.Gs33[x,y,z] +=  Ds33 * self.ts * self.StPl * self.axis3
        self.Gs12[x,y,z] +=  Ds12 * self.ts * self.StPl * self.axis1 * self.axis2
        self.Gs23[x,y,z] +=  Ds23 * self.ts * self.StPl * self.axis2 * self.axis3
        self.Gs31[x,y,z] +=  Ds31 * self.ts * self.StPl * self.axis3 * self.axis1

        return self

    def UpdateVelocity(self,x,y,z):
        #Updates velocity based off of previous velocities and delta velocity times time
        #
        # Inputs: Coordinates of cube in question.  Assumed last time step
        #
        # Output: updated self.Gs matrix

        Dv1, Dv2, Dv3 = self.DeltaVelocity(x,y,z)

        self.Gv1[x,y,z] += Dv1 * self.ts * self.VePl * self.axis1
        self.Gv2[x,y,z] += Dv2 * self.ts * self.VePl * self.axis2
        self.Gv3[x,y,z] += Dv3 * self.ts * self.VePl * self.axis3

        return self
    
    def StepStresses(self,xx1=1,xx2=1,xx3=1):
        for i in range(self.Max1-1):
            x1 = i + xx1
            for j in range(self.Max2-1):
                x2 = j + xx2
                for k in range(self.Max3-1):
                    x3 = k + xx3
                    self.UpdateStresses(x1,x2,x3)
        '''
        self.Gs11[0,:,:] = -self.Gs11[1,:,:]
        self.Gs22[:,0,:] = -self.Gs22[:,1,:]
        self.Gs33[:,:,0] = -self.Gs33[:,:,1]
        self.Gs11[self.Max1,:,:] = -self.Gs11[self.Max1-1,:,:]
        self.Gs22[:,self.Max2,:] = -self.Gs22[:,self.Max2-1,:]
        self.Gs33[:,:,self.Max3] = -self.Gs33[:,:,self.Max3-1]
        self.Gs12[self.Max1,:,:]=0
        self.Gs23[self.Max1,:,:]=0
        self.Gs31[self.Max1,:,:]=0
        self.Gs12[:,self.Max2,:]=0
        self.Gs23[:,self.Max2,:]=0
        self.Gs31[:,self.Max2,:]=0
        self.Gs12[:,:,self.Max3]=0
        self.Gs23[:,:,self.Max3]=0
        self.Gs31[:,:,self.Max3]=0
        self.Gs12[0,:,:]=0
        self.Gs23[0,:,:]=0
        self.Gs31[0,:,:]=0
        self.Gs12[:,0,:]=0
        self.Gs23[:,0,:]=0
        self.Gs31[:,0,:]=0
        self.Gs12[:,:,0]=0
        self.Gs23[:,:,0]=0
        self.Gs31[:,:,0]=0
                '''

    def StepVelocities(self,xx1=1,xx2=1,xx3=1):
        for i in range(self.Max1-1):
            x1 = i + xx1
            for k in range(self.Max2-1):
                x2 = k + xx2
                for j in range(self.Max3-1):
                    x3 = j + xx3
                    self.UpdateVelocity(x1,x2,x3)
        #Boundary, floor is stationary
        self.Gv1[:,:,0]=0
        self.Gv2[:,:,0]=0
        self.Gv3[:,:,0]=0

    def ForcingFunctionWave(self, t, Hz = 40000, EP=100.0, size=0.04, Odim = 2, Dir=1):
        # Adds stresses from a force to the stress grid
        # Initially assumed a single force of a small plate sinosoidal ultrasound emitter.  More to be added later
        # 
        # Input:   t is the time
        #          Hz is the frequency of the signal
        #          EP is the force / density
        #          size is in meters the size of the emitter plate
        #          Odim is the dimmension the actio is workig on
        #          dir is for 1 max serface, 2 for middle of dimmension, other for min surface
        #
        # Outputs: no direct outputs, last time step stress is updated

        frequency = 1/Hz
        EmitterPreasure = EP / 2.0
        
        ##run for two periods and then stop:
        #if 2.0 / frequency < t:
            
        EmitterWidth = size / self.ds
        EmitterWidth = int(EmitterWidth)
        if EmitterWidth == 0: EmitterWidth = 2

        #emitter placed in middle of top face
        
        #StartX = int(self.MaxX / 2 - EmitterWidth / 2)
        #StartZ = int(self.MaxZ / 2 - EmitterWidth / 2)

        Temp = np.zeros((EmitterWidth,EmitterWidth))

        Temp[:,:] = np.sin(t/frequency*2*3.1415926) * EmitterPreasure

        if Odim == 0:
            Start0 = int((self.Max2 / 2) - (EmitterWidth / 2))
            Start1 = int((self.Max3 / 2) - (EmitterWidth / 2))
            Start2 = int((self.Max1 / 2) - (1))
            if Dir ==1:
                self.Gv1[self.Max1-1,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
                self.Gv1[self.Max1-2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
            elif Dir ==2:
                self.Gv1[Start2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
                self.Gv1[Start2-1,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
            else:
                self.Gv1[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
                self.Gv1[3,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
        elif Odim == 1:
            Start0 = int((self.Max1 / 2) - (EmitterWidth / 2))
            Start1 = int((self.Max3 / 2) - (EmitterWidth / 2))
            Start2 = int((self.Max2 / 2) - (1))
            if Dir ==1:
                self.Gv2[Start0:Start1+EmitterWidth,self.Max2-1,Start0:Start1+EmitterWidth] = Temp
                self.Gv2[Start0:Start1+EmitterWidth,self.Max2-2,Start0:Start1+EmitterWidth] = Temp
            elif Dir == 2:
                self.Gv2[Start0:Start1+EmitterWidth,Start2,Start0:Start1+EmitterWidth] = Temp
                self.Gv2[Start0:Start1+EmitterWidth,Start2-1,Start0:Start1+EmitterWidth] = Temp
            else:
                self.Gv2[Start0:Start1+EmitterWidth,1,Start0:Start1+EmitterWidth] = Temp
                self.Gv2[Start0:Start1+EmitterWidth,2,Start0:Start1+EmitterWidth] = Temp
        else:
            Start0 = int((self.Max1 / 2) - (EmitterWidth / 2))
            Start1 = int((self.Max2 / 2) - (EmitterWidth / 2))
            Start2 = int((self.Max3 / 2) - (1))
            if Dir ==1:
                self.Gv3[Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,self.Max3-1] = Temp
                self.Gv3[Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,self.Max3-2] = Temp
            elif Dir ==2:
                self.Gv3[Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,Start2] = Temp
                self.Gv3[Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,Start2-1] = Temp
            else:
                self.Gv3[Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,1] = Temp
                self.Gv3[Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,2] = Temp
        return np.sum(Temp)

    def ForcingFunctionImpulse(self, force = 100.0, emitter = 0.01, Odim = 2, Dir=1,CornerCut = 0):
        # Adds stresses from a force to the stress grid
        # Initially assumed a single force of a small plate sinosoidal ultrasound emitter.  More to be added later
        # 
        # Input is the time, impulse, square emitter size, what face it is on, and if it is top or bottom
        #
        # Outputs: no direct outputs, last time step stress is updated

        EmitterWidth = emitter * self.ids
        EmitterWidth = int(EmitterWidth)
        if EmitterWidth <= 3: EmitterWidth = 3

        Temp = np.zeros((EmitterWidth,EmitterWidth))
        
        EmmitterArea = EmitterWidth **2 - (4 * ((CornerCut*(CornerCut + 1)/2)))

        Temp[:,:] = -(force / EmmitterArea) / 2.0
        if CornerCut >= EmitterWidth:
            print('Corner Cut too large')
            CornerCut = 0
        if CornerCut > 0:
            CutMatrix = np.zeros((CornerCut,CornerCut))
            for j in range(CornerCut):
                for k in range(CornerCut):
                    if j < k: CutMatrix[j,k]=1
            Temp[EmitterWidth-CornerCut:EmitterWidth,0:CornerCut] *= CutMatrix
            Temp[EmitterWidth-CornerCut:EmitterWidth,EmitterWidth-CornerCut:EmitterWidth] *= np.flip(CutMatrix,1)
            Temp[0:CornerCut,0:CornerCut] *= np.flip(CutMatrix,0)
            Temp[0:CornerCut,EmitterWidth-CornerCut:EmitterWidth] *= np.flip(np.flip(CutMatrix,0),1)
        
        if Odim == 0:
            Start0 = int((self.MaxY / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxZ / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.GvX[self.MaxX,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp / self.Gp[0,self.MaxX,Start0,Start0]
            else:
                self.GvX[0,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp / self.Gp[0,0,self.MaxX,Start0,Start0]
        elif Odim == 1:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxZ / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.GvY[Start0:Start1+EmitterWidth,self.MaxY,Start0:Start1+EmitterWidth] = Temp / self.Gp[0,Start0,self.MaxY,Start0]
            else:
                self.GvY[Start0:Start1+EmitterWidth,0,Start0:Start1+EmitterWidth] = Temp / self.Gp[0,Start0,self.MaxY,Start0]
        else:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxY / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.GvZ[Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,self.MaxZ] = Temp / self.Gp[0,Start0,Start0,self.MaxZ]
            else:
                self.GvZ[Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,0] = Temp / self.Gp[0,Start0,Start0,self.MaxZ]

        EmitterWidth -= 2
        
        Temp = np.zeros((EmitterWidth,EmitterWidth))
        
        EmmitterArea = EmitterWidth **2 - (4 * ((CornerCut*(CornerCut + 1)/2)))

        Temp[:,:] = -(force / EmmitterArea) / 2.0
        if CornerCut >= EmitterWidth:
            print('Corner Cut too large')
            CornerCut = 0
        if CornerCut > 0:
            CutMatrix = np.zeros((CornerCut,CornerCut))
            for j in range(CornerCut):
                for k in range(CornerCut):
                    if j < k: CutMatrix[j,k]=1
            Temp[EmitterWidth-CornerCut:EmitterWidth,0:CornerCut] *= CutMatrix
            Temp[EmitterWidth-CornerCut:EmitterWidth,EmitterWidth-CornerCut:EmitterWidth] *= np.flip(CutMatrix,1)
            Temp[0:CornerCut,0:CornerCut] *= np.flip(CutMatrix,0)
            Temp[0:CornerCut,EmitterWidth-CornerCut:EmitterWidth] *= np.flip(np.flip(CutMatrix,0),1)

        if Odim == 0:
            Start0 = int((self.MaxY / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxZ / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.GvX[self.MaxX-1,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp / self.mu
            else:
                self.GvX[1,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp / self.mu
        elif Odim == 1:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxZ / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.GvY[1,Start0:Start1+EmitterWidth,self.MaxY-1,Start0:Start1+EmitterWidth] = Temp / self.mu
            else:
                self.GvY[2,Start0:Start1+EmitterWidth,1,Start0:Start1+EmitterWidth] = Temp / self.mu
        else:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxY / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.GvZ[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,self.MaxZ-1] = Temp / self.mu
            else:
                self.GvZ[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,1] = Temp / self.mu
 

        return self

    def VelocityCut(self, Dimm = 0, Component = -1, location = -1):
        #Gets a cut of the magnitude of velocity in a plane defined by the two dimmension, cut at the location given in location
        #
        # Inputs: Dimm is the dimmension that is orthoganl to the plane that will be reutrned.
        #         Location if given will be the cut allong the thrid dimmension that will be returned.  
        #               If not given, -1 is default, it will signal using the centroid by default
        #
        # Outputs: 2 dimesional matrix with the combined velocity in the 2 given dimensions

        Results = []

        if Dimm == 0:
            if location > self.Max1 or location < 0: location = -1
            if location == -1:
                location = int(self.Max1 / 2)
            Component0 = np.matrix((self.Max2, self.Max3))
            Component1 = np.matrix((self.Max2, self.Max3))
            Component2 = np.matrix((self.Max2, self.Max3))

            Component0 = self.Gv1[location,:,:]
            Component1 = self.Gv2[location,:,:]
            Component2 = self.Gv3[location,:,:]
            
        if Dimm == 1:
            if location > self.Max2 or location < 0: location = -1
            if location == -1:
                location = int(self.Max2 / 2)
            Component0 = np.matrix((self.Max1, self.Max3))
            Component1 = np.matrix((self.Max1, self.Max3))
            Component2 = np.matrix((self.Max1, self.Max3))

            Component0 = self.Gv1[:,location,:]
            Component1 = self.Gv2[:,location,:]
            Component2 = self.Gv3[:,location,:]
 
        if Dimm == 2:
            if location > self.Max3 or location < 0: location = -1
            if location == -1:
                location = int(self.Max3 / 2)
            Component0 = np.matrix((self.Max1, self.Max3))
            Component1 = np.matrix((self.Max1, self.Max3))
            Component2 = np.matrix((self.Max1, self.Max3))

            Component0 = self.Gv1[:,:,location]
            Component1 = self.Gv2[:,:,location]
            Component2 = self.Gv3[:,:,location]
            
        #Results = Component1 #+ Component2
        if Component == -1:
            Results = np.sqrt(Component0**2+Component1**2+Component2**2)
        elif Component == 0:
            Results = Component0
        elif Component == 1:
            Results = Component1
        elif Component == 2:
            Results = Component2
        else:
            print('Invalid component')
            Results = np.sqrt(Component0**2+Component1**2+Component2**2)


        return Results

    def StressCut(self, Dimm=0, Dimm1=0, Dimm2=0, location = -1):
        #Gets a cut of the magnitude of stress on the give face and direction, cut at the location given in location
        #
        # Inputs: Dimm is the dimmension to make the cut across
        #         Dimm1 is the dimmension that is orthoganl to the plane that will be reutrned.
        #         Dimm2 is the vector of the stress
        #         Location if given will be the cut allong the thrid dimmension that will be returned.  
        #               If not given, -1 is default, it will signal using the centroid by default
        #
        # Outputs: 2 dimesional matrix with the combined velocity in the 2 given dimensions

        Results = []

        if Dimm == 0:
            if location > self.MaxX or location < 0: location = -1
            if location == -1:
                location = int(self.MaxX / 2)
            Component1 = np.matrix((self.MaxY, self.MaxZ))
            
            Component1 = self.Gs[Dimm1,Dimm2,location,:,:]
            
            Results = Component1
            #Results = np.sqrt(np.add(np.multiply(Component1,Component1),np.multiply(Component2,Component2)))

        if Dimm == 1:
            if location > self.MaxY or location < 0: location = -1
            if location == -1:
                location = int(self.MaxY / 2)
            Component1 = np.matrix((self.MaxX, self.MaxZ))
            
            Component1 = self.Gs[Dimm1, Dimm2,:,location,:]
            
            Results = Component1
            # Results = np.sqrt(np.add(np.multiply(Component1,Component1),np.multiply(Component2,Component2)))
            
        if Dimm == 2:
            if location > self.MaxZ or location < 0: location = -1
            if location == -1:
                location = int(self.MaxZ / 2)
            Component1 = np.matrix((self.MaxX, self.MaxY))
            
            Component1 = self.Gs[Dimm1,Dimm2,:,:,location]
            
            Results = Component1
            #Results = np.sqrt(np.add(np.multiply(Component1,Component1),np.multiply(Component2,Component2)))
            

        return Results

    def VelocitySave(self, Dimm = -1):
        Component0 = self.Gv1[:,:,:]
        Component1 = self.Gv2[:,:,:]
        Component2 = self.Gv3[:,:,:]
        if Dimm == -1:
            Results = np.sqrt(Component0**2+Component1**2+Component2**2)
        elif Dimm == 0:
            Results = np.sqrt(Component0**2)
        elif Dimm == 1:
            Results = np.sqrt(Component1**2)
        elif Dimm == 2:
            Results = np.sqrt(Component2**2)
        else:
            print('Error, unknown Dimm', Dimm)
        
        return Results


    def CheckStressBoundary(self,x,y,z,Ds):
        #checks to see if a grid is at a boundary, and if so, adjusts boundary conditions appropriately
        #
        # Inputs: x,y,z coordinates of cube in question
        #
        # Outputs: Updated (if boundary) delta stress matrix

        #at front and back faces, stresses perpendicular to the face are 0:
        if x == 0 or x == self.MaxX-1:
            Ds[0,0]=0
            #Ds[0,1]=0
            #Ds[0,2]=0
        
        #at top face, stresses perpendicular to the face are 0:
        if y == self.MaxY-1:
            #Ds[1,0]=0
            Ds[1,1]=0
            #Ds[1,2]=0
        
        if y == 0:
            Ds[1,1]=-self.Gs[1,1,x,y+1,z]
                
        
        #at side faces, stresses perpendicular to the face are 0:
        if z == 0 or z == self.MaxZ-1:
            #Ds[2,0]=0
            #Ds[2,1]=0
            Ds[2,2]=0
        
        return Ds

    def CheckVelocityBoundary(self,x,y,z,Dv):
        #checks to see if a grid is at a boundary, and if so, adjusts boundary conditions appropriately
        #
        # Inputs: x,y,z coordinates of cube in question
        #
        # Outputs: Updated (if boundary) delta velocity vector
        
        #lower boundary is fixed, velocity is Zero
        if y == 0:
            Dv[1]=0
            if x <= 10 or x >= self.MaxX-10:
                Dv[0] = 0
                Dv[2] = 0

        if x == 0 or y == 0 or z == 0 or x == self.MaxX+1 or y == self.MaxY+1 or z == self.MaxZ+1:
            Dv[:]=0

        return Dv