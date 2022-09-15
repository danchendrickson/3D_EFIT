#import limited
from zlib import Z_DEFAULT_COMPRESSION
import numpy as np
#import visvis as vv
class EFIT:
    ts = 0
    ds = 0

    def __init__(self, xGrid, yGrid, zGrid, tStep, dStep,xdir=1,ydir=1,zdir=1):
        #Initialize with the size of the grid in X, Y, and Z number of nodes.  The distance step, and the time step


        #Velocity grid with 3 part vector at each node point, for 2 time steps
        self.GridShapeV = (3,xGrid,yGrid,zGrid)
        #Stresses with 3 part vector in each of 3 part phaces, at each nod epoint, for 2 time steps
        self.GridShapeS = (3,3,xGrid,yGrid,zGrid)
        #materials property gird.  Initially 3 properties needed: density, Lame 1, Lame 2
        self.GridShapeP = (3,xGrid,yGrid,zGrid)

        #total gird size
        self.GridPoints = (xGrid)*(yGrid)*(zGrid)
        
        #either 1 or -1, will be used in step calculations, allowing expiremnts with axis directions
        self.xDir = xdir
        self.yDir = ydir
        self.zDir = zdir
        
        #define empty grid for the 3 directions of velocity for 2 time steps
        self.Gv = np.zeros(3*self.GridPoints,dtype="double").reshape(*self.GridShapeV)
        #define empty grid for the 3 directions of stress on 3 dimmensions of faces for 2 time steps
        self.Gs = np.zeros(3*3*self.GridPoints,dtype="double").reshape(*self.GridShapeS)
        #define empty grid for the 3 scalar material properties at each node point.  Can honly hold scalar properties
        #Assumed properties are density, Lame 1, Lame 2
        self.Gp = np.zeros(3*self.GridPoints,dtype="double").reshape(*self.GridShapeP)
        
        self.MaxX = xGrid - 1
        self.MaxY = yGrid - 1
        self.MaxZ = zGrid - 1

        self.ds = dStep
        self.ts = tStep

        self.ids = 1.0 / dStep

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
    
    def DeltaStress(self,x,y,z):
        #Gets the change in the stresses per time at a certain coordinate juncture 
        #
        # Inputs: x,y,z coordinates of the cube in question.  Last time is assumed
        #
        # Outputs: 6 dimmensions of stress.  
        
        #Calculated stresses based on 3.55
        Ds = np.zeros((3,3))
        
        Lame1=self.Gp[1,x,y,z]
        Lame2=self.Gp[2,x,y,z]

        mm = [self.xDir, self.yDir,self.zDir]
        d1=[0,0,0]
        d2=[0,0,0]
        d3=[0,0,0]
        for i in range(3):
            d1[i]=-1*mm[i]
            d2[(i+1) % 3] = -1 * mm[(i+1) % 3]
            d3[(i+2) % 3] = -1 * mm[(i+2) % 3]
            uu=np.roll(mm,i)
            Ds[i,i] =  ((self.ids) *
                ((Lame1+2*Lame2)*(uu[0]*self.Gv[0,x,y,z]-uu[0]*self.Gv[0,x+d1[0],y+d1[1],z+d2[2]]) +
                    Lame1*((uu[1]*self.Gv[1,x,y,z]-uu[1]*self.Gv[1,x+d2[0],y+d2[1],z+d2[2]])
                           +(uu[2]*self.Gv[2,x,y,z]-uu[2]*self.Gv[2,x+d3[0],y+d3[1],z+d3[2]]))
                    )
                )

            for j in range(3):
                h1 = [0,0,0]
                h2=h1
                h1[i]=1*mm[i]
                h2[j]=1*mm[j]

                if i != j:
                    Ds[i,j] =  (self.ids) * (
                        4/( (1/self.Gp[2,x,y,z])
                            +(1/self.Gp[2,x+h1[0],y+h1[1],z+h1[2]])
                            +(1/self.Gp[2,x+h2[0],y+h2[1],z+h2[2]])
                            +(1/self.Gp[2,x+h1[0]+h2[0],y+h1[1]+h2[1],z+h1[2]+h2[2]])
                            )
                        ) * (
                            (self.Gv[i,x+h2[0],y+h2[1],z+h2[2]]-self.Gv[i,x,y,z])
                            + (self.Gv[j,x+h1[0],y+h1[1],z+h1[2]]-self.Gv[j,x,y,z])
                        )
                        
        

        return Ds

    def DeltaVelocity(self, x,y,z):
        #Gets the change in the velocity per time at a certain coordinate juncture 
        #
        # Inputs: x,y,z coordinates of the cube in question.  Last time is assumed
        #
        # Outputs: 3 dimmensions of velocity.  
        
        #Calculated velocity based on 3.54

        DV = np.zeros(3)

        for i in range(3):
            iS = [0,0,0]
            iM = [-1,-1,-1]
            dS = [self.xDir, self.yDir, self.zDir]
            iS[i]=1 * dS[i]
            iM[i] = 0
            iM = np.multiply(iM, dS)
            DV[0] = ((self.ids ) *
                (2 / (self.Gp[0,x,y,z]+self.Gp[0,x+iS[0],y+iS[1],z+iS[2]])) *
                (
                    (self.Gs[0,i,x+iS[0],y+iS[1],z+iS[2]] - self.Gs[0,i,x+iM[0],y+iM[1],z+iM[2]]) +
                    (self.Gs[1,i,x+iS[0],y+iS[1],z+iS[2]] - self.Gs[1,i,x+iM[0],y+iM[1],z+iM[2]]) + 
                    (self.Gs[2,i,x+iS[0],y+iS[1],z+iS[2]] - self.Gs[2,i,x+iM[0],y+iM[1],z+iM[2]])
                )
            )

        return DV
    
    def UpdateStresses(self, x,y,z):
        #Updates velocity based off of previous velocities and delta velocity times time
        #
        # Inputs: Coordinates of cube in question.  Assumed last time step
        #
        # Output: updated self.Gs matrix
        
        delS = self.DeltaStress(x,y,z)

        for i in range(3):
                for j in range(3):
                    self.Gs[j,i,x,y,z] +=  delS[j,i] * self.ts
        #self.Gs[:,:,x,y,z] += delS[:,:] * self.ts

        return self

    def UpdateVelocity(self,x,y,z):
        #Updates velocity based off of previous velocities and delta velocity times time
        #
        # Inputs: Coordinates of cube in question.  Assumed last time step
        #
        # Output: updated self.Gs matrix

        delV = self.DeltaVelocity(x,y,z)

        for i in range(3):
            self.Gv[i,x,y,z] += delV[i] * self.ts

        return self
    
    def StepStresses(self,xx=1,yy=1,zz=1):
        for i in range(self.MaxX-1):
            x = i + xx
            for j in range(self.MaxY-1):
                y = j + yy
                for k in range(self.MaxZ-1):
                    z = k + zz
                    self.UpdateStresses(x,y,z)
        
        self.Gs[0,0,0,:,:] = -self.Gs[0,0,1,:,:]
        self.Gs[1,1,:,0,:] = -self.Gs[1,1,:,1,:]
        self.Gs[2,2,:,:,0] = -self.Gs[2,2,:,:,1]
        self.Gs[0,0,self.MaxX,:,:] = -self.Gs[0,0,self.MaxX-1,:,:]
        self.Gs[1,1,:,self.MaxY,:] = -self.Gs[1,1,:,self.MaxY-1,:]
        self.Gs[2,2,:,:,self.MaxZ] = -self.Gs[2,2,:,:,self.MaxZ-1]
        self.Gs[0,1,self.MaxX-1,:,:]=0
        self.Gs[0,2,self.MaxX-1,:,:]=0
        self.Gs[1,0,self.MaxX-1,:,:]=0
        self.Gs[1,2,self.MaxX-1,:,:]=0
        self.Gs[1,0,:,self.MaxY-1,:]=0
        self.Gs[1,2,:,self.MaxY-1,:]=0
        self.Gs[0,1,:,self.MaxY-1,:]=0
        self.Gs[2,1,:,self.MaxY-1,:]=0
        self.Gs[2,0,:,:,self.MaxZ-1]=0
        self.Gs[2,1,:,:,self.MaxZ-1]=0
        self.Gs[0,2,:,:,self.MaxZ-1]=0
        self.Gs[1,2,:,:,self.MaxZ-1]=0
        self.Gs[0,1,0,:,:]=0
        self.Gs[0,2,0,:,:]=0
        self.Gs[1,0,0,:,:]=0
        self.Gs[1,2,0,:,:]=0
        self.Gs[1,0,:,0,:]=0
        self.Gs[1,2,:,0,:]=0
        self.Gs[0,1,:,0,:]=0
        self.Gs[2,1,:,0,:]=0
        self.Gs[2,0,:,:,0]=0
        self.Gs[2,1,:,:,0]=0
        self.Gs[0,2,:,:,0]=0
        self.Gs[1,2,:,:,0]=0
        

    def StepVelocities(self,xx=1,yy=1,zz=1):
        for i in range(self.MaxX-1):
            x = i + xx
            for k in range(self.MaxY-1):
                y = k + yy
                for j in range(self.MaxZ-1):
                    z = j + zz
                    self.UpdateVelocity(x,y,z)
        self.Gv[:,:,:,0]=0

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

        frequency = Hz
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

        Temp[:,:] = np.sin(frequency * t) * EmitterPreasure

        if Odim == 0:
            Start0 = int((self.MaxY / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxZ / 2) - (EmitterWidth / 2))
            Start2 = int((self.MaxX / 2) - (1))
            if Dir ==1:
                self.Gv[0,self.MaxX-1,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
                self.Gv[0,self.MaxX-2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
            elif Dir ==2:
                self.Gv[0,Start2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
                self.Gv[0,Start2-1,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
            else:
                self.Gv[0,2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
                self.Gv[0,3,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
        elif Odim == 1:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxZ / 2) - (EmitterWidth / 2))
            Start2 = int((self.MaxY / 2) - (1))
            if Dir ==1:
                self.Gv[1,Start0:Start1+EmitterWidth,self.MaxY-1,Start0:Start1+EmitterWidth] = Temp
                self.Gv[1,Start0:Start1+EmitterWidth,self.Maxy-2,Start0:Start1+EmitterWidth] = Temp
            elif Dir == 2:
                self.Gv[1,Start0:Start1+EmitterWidth,Start2,Start0:Start1+EmitterWidth] = Temp
                self.Gv[1,Start0:Start1+EmitterWidth,Start2-1,Start0:Start1+EmitterWidth] = Temp
            else:
                self.Gv[2,Start0:Start1+EmitterWidth,1,Start0:Start1+EmitterWidth] = Temp
                self.Gv[2,Start0:Start1+EmitterWidth,2,Start0:Start1+EmitterWidth] = Temp
        else:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxY / 2) - (EmitterWidth / 2))
            Start2 = int((self.MaxZ / 2) - (1))
            if Dir ==1:
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,self.MaxZ-1] = Temp
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,self.MaxZ-2] = Temp
            elif Dir ==2:
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,Start2] = Temp
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,Start2-1] = Temp
            else:
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,1] = Temp
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,2] = Temp
        return self

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
                self.Gv[0,self.MaxX,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp / self.Gp[0,self.MaxX,Start0,Start0]
            else:
                self.Gv[0,0,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp / self.Gp[0,0,self.MaxX,Start0,Start0]
        elif Odim == 1:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxZ / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.Gv[1,Start0:Start1+EmitterWidth,self.MaxY,Start0:Start1+EmitterWidth] = Temp / self.Gp[0,Start0,self.MaxY,Start0]
            else:
                self.Gv[2,Start0:Start1+EmitterWidth,0,Start0:Start1+EmitterWidth] = Temp / self.Gp[0,Start0,self.MaxY,Start0]
        else:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxY / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,self.MaxZ] = Temp / self.Gp[0,Start0,Start0,self.MaxZ]
            else:
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,0] = Temp / self.Gp[0,Start0,Start0,self.MaxZ]

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
                self.Gv[0,self.MaxX-1,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp / self.Gp[0,self.MaxX,Start0,Start0]
            else:
                self.Gv[0,1,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp / self.Gp[0,1,self.MaxX,Start0,Start0]
        elif Odim == 1:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxZ / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.Gv[1,Start0:Start1+EmitterWidth,self.MaxY-1,Start0:Start1+EmitterWidth] = Temp / self.Gp[0,Start0,self.MaxY,Start0]
            else:
                self.Gv[2,Start0:Start1+EmitterWidth,1,Start0:Start1+EmitterWidth] = Temp / self.Gp[0,Start0,self.MaxY,Start0]
        else:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxY / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,self.MaxZ-1] = Temp / self.Gp[0,Start0,Start0,self.MaxZ]
            else:
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,1] = Temp / self.Gp[0,Start0,Start0,self.MaxZ]
 

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
            if location > self.MaxX or location < 0: location = -1
            if location == -1:
                location = int(self.MaxX / 2)
            Component0 = np.matrix((self.MaxY, self.MaxZ))
            Component1 = np.matrix((self.MaxY, self.MaxZ))
            Component2 = np.matrix((self.MaxY, self.MaxZ))

            Component0 = self.Gv[0,location,:,:]
            Component1 = self.Gv[1,location,:,:]
            Component2 = self.Gv[2,location,:,:]
            
        if Dimm == 1:
            if location > self.MaxY or location < 0: location = -1
            if location == -1:
                location = int(self.MaxY / 2)
            Component0 = np.matrix((self.MaxX, self.MaxZ))
            Component1 = np.matrix((self.MaxX, self.MaxZ))
            Component2 = np.matrix((self.MaxX, self.MaxZ))

            Component0 = self.Gv[0,:,location,:]
            Component1 = self.Gv[1,:,location,:]
            Component2 = self.Gv[2,:,location,:]
 
        if Dimm == 2:
            if location > self.MaxZ or location < 0: location = -1
            if location == -1:
                location = int(self.MaxZ / 2)
            Component0 = np.matrix((self.MaxX, self.MaxY))
            Component1 = np.matrix((self.MaxX, self.MaxY))
            Component2 = np.matrix((self.MaxX, self.MaxY))

            Component0 = self.Gv[0,:,:,location]
            Component1 = self.Gv[1,:,:,location]
            Component2 = self.Gv[2,:,:,location]
            
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
        Component0 = self.Gv[0,:,:,:]
        Component1 = self.Gv[1,:,:,:]
        Component2 = self.Gv[2,:,:,:]
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
