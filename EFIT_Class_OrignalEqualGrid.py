import numpy as np
import visvis as vv
class EFIT:
    ts = 0
    ds = 0

    def __init__(self, xGrid, yGrid, zGrid, tStep, dStep):
        #Initialize with the size of the grid in X, Y, and Z number of nodes.  The distance step, and the time step


        #Velocity grid with 3 part vector at each node point, for 2 time steps
        self.GridShapeV = (3,xGrid,yGrid,zGrid)
        #Stresses with 3 part vector in each of 3 part phaces, at each nod epoint, for 2 time steps
        self.GridShapeS = (3,3,xGrid,yGrid,zGrid)
        #materials property gird.  Initially 3 properties needed: density, Lame 1, Lame 2
        self.GridShapeP = (3,xGrid,yGrid,zGrid)
     
        self.GridPoints = (xGrid)*(yGrid)*(zGrid)
        
        
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

        #experimental add for grid offsetting
        # x+=1
        # y+=1
        # z+=1
        #end experiment

        BCs = 0
        if x == self.MaxX: BCs+=1
        if x == 0: BCs+=1
        if y == self.MaxY: BCs+=1
        if y == 0: BCs+=1
        if z == self.MaxZ: BCs+=1
        if z == 0: BCs+=1

        if BCs == 0:
            Ds[0,0] =  ((self.ids) *
                ((Lame1+2*Lame2)*(self.Gv[0,x,y,z]-self.Gv[0,x-1,y,z]) +
                    Lame1*(self.Gv[1,x,y,z]-self.Gv[1,x,y-1,z]+self.Gv[2,x,y,z]-self.Gv[2,x,y,z-1])
                    )
                )
            Ds[1,1] =  ((self.ids) *
                ((Lame1+2*Lame2)*(self.Gv[1,x,y,z]-self.Gv[1,x,y-1,z]) +
                    Lame1*(self.Gv[0,x,y,z]-self.Gv[0,x-1,y,z]+self.Gv[2,x,y,z]-self.Gv[2,x,y,z-1])
                    )
                )
            Ds[2,2] =  ((self.ids) *
                ((Lame1+2*Lame2)*(self.Gv[2,x,y,z]-self.Gv[2,x,y,z-1]) +
                    Lame1*(self.Gv[0,x,y,z]-self.Gv[0,x-1,y,z]+self.Gv[1,x,y,z]-self.Gv[1,x,y-1,z])
                    )
                )
            Ds[0,1] =  (
                (self.ids) *
                (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x+1,y,z])+(1/self.Gp[2,x,y+1,z])+(1/self.Gp[2,x+1,y+1,z]))) *
                (self.Gv[0,x,y+1,z]-self.Gv[0,x,y,z] + self.Gv[1,x+1,y,z]-self.Gv[1,x,y,z] )
                )
            Ds[0,2] =  (
                (self.ids) *
                (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x+1,y,z])+(1/self.Gp[2,x,y,z+1])+(1/self.Gp[2,x+1,y,z+1]))) *
                (self.Gv[0,x,y,z+1]-self.Gv[0,x,y,z] +self.Gv[2,x+1,y,z]-self.Gv[2,x,y,z] )
                )
            Ds[1,2] =  (
                (self.ids) *
                (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x,y+1,z])+(1/self.Gp[2,x,y,z+1])+(1/self.Gp[2,x,y+1,z+1]))) *
                (self.Gv[1,x,y,z+1]-self.Gv[1,x,y,z] +self.Gv[2,x,y+1,z]-self.Gv[2,x,y,z] )
                )
        elif BCs == 1:
            if x == self.MaxX or x == 0: 
                Ds[1,2] = 0
                if x == 0:
                    Ds[0,1] =  (
                        (self.ids) *
                        (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x+1,y,z])+(1/self.Gp[2,x,y+1,z])+(1/self.Gp[2,x+1,y+1,z]))) *
                        (self.Gv[0,x,y+1,z]-self.Gv[0,x,y,z] + self.Gv[1,x+1,y,z]-self.Gv[1,x,y,z] )
                        )
                    Ds[0,2] =  (
                        (self.ids) *
                        (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x+1,y,z])+(1/self.Gp[2,x,y,z+1])+(1/self.Gp[2,x+1,y,z+1]))) *
                        (self.Gv[0,x,y,z+1]-self.Gv[0,x,y,z] +self.Gv[2,x+1,y,z]-self.Gv[2,x,y,z] )
                        )
                else:
                    Ds[0,1] =  (
                        (self.ids) *
                        (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x-1,y,z])+(1/self.Gp[2,x,y+1,z])+(1/self.Gp[2,x-1,y+1,z]))) *
                        (self.Gv[0,x,y+1,z]-self.Gv[0,x,y,z] + self.Gv[1,x-1,y,z]-self.Gv[1,x,y,z] )
                        )
                    Ds[0,2] =  (
                        (self.ids) *
                        (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x-1,y,z])+(1/self.Gp[2,x,y,z+1])+(1/self.Gp[2,x-1,y,z+1]))) *
                        (self.Gv[0,x,y,z+1]-self.Gv[0,x,y,z] +self.Gv[2,x-1,y,z]-self.Gv[2,x,y,z] )
                        )            
            elif y == self.MaxY or y ==0:
                Ds[0,2] = 0
                if y == 0:
                    Ds[0,1] =  (
                        (self.ids) *
                        (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x+1,y,z])+(1/self.Gp[2,x,y+1,z])+(1/self.Gp[2,x+1,y+1,z]))) *
                        (self.Gv[0,x,y+1,z]-self.Gv[0,x,y,z] + self.Gv[1,x+1,y,z]-self.Gv[1,x,y,z] )
                        )
                    Ds[1,2] =  (
                        (self.ids) *
                        (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x,y+1,z])+(1/self.Gp[2,x,y,z+1])+(1/self.Gp[2,x,y+1,z+1]))) *
                        (self.Gv[1,x,y,z+1]-self.Gv[1,x,y,z] +self.Gv[2,x,y+1,z]-self.Gv[2,x,y,z] )
                        )
                else:
                    Ds[0,1] =  (
                        (self.ids) *
                        (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x+1,y,z])+(1/self.Gp[2,x,y-1,z])+(1/self.Gp[2,x+1,y-1,z]))) *
                        (self.Gv[0,x,y-1,z]-self.Gv[0,x,y,z] + self.Gv[1,x+1,y,z]-self.Gv[1,x,y,z] )
                        )
                    Ds[1,2] =  (
                        (self.ids) *
                        (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x,y-1,z])+(1/self.Gp[2,x,y,z+1])+(1/self.Gp[2,x,y-1,z+1]))) *
                        (self.Gv[1,x,y,z+1]-self.Gv[1,x,y,z] +self.Gv[2,x,y-1,z]-self.Gv[2,x,y,z] )
                        )
            elif z == self.MaxZ or z == 0:
                Ds[0,1] = 0
                if z == 0:
                    Ds[0,2] =  (
                        (self.ids) *
                        (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x+1,y,z])+(1/self.Gp[2,x,y,z+1])+(1/self.Gp[2,x+1,y,z+1]))) *
                        (self.Gv[0,x,y,z+1]-self.Gv[0,x,y,z] +self.Gv[2,x+1,y,z]-self.Gv[2,x,y,z] )
                        )
                    Ds[1,2] =  (
                        (self.ids) *
                        (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x,y+1,z])+(1/self.Gp[2,x,y,z+1])+(1/self.Gp[2,x,y+1,z+1]))) *
                        (self.Gv[1,x,y,z+1]-self.Gv[1,x,y,z] +self.Gv[2,x,y+1,z]-self.Gv[2,x,y,z] )
                        )
                else:
                    Ds[0,2] =  (
                        (self.ids) *
                        (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x+1,y,z])+(1/self.Gp[2,x,y,z-1])+(1/self.Gp[2,x+1,y,z-1]))) *
                        (self.Gv[0,x,y,z-1]-self.Gv[0,x,y,z] +self.Gv[2,x+1,y,z]-self.Gv[2,x,y,z] )
                        )
                    Ds[1,2] =  (
                        (self.ids) *
                        (4/((1/self.Gp[2,x,y,z])+(1/self.Gp[2,x,y+1,z])+(1/self.Gp[2,x,y,z-1])+(1/self.Gp[2,x,y+1,z-1]))) *
                        (self.Gv[1,x,y,z-1]-self.Gv[1,x,y,z] +self.Gv[2,x,y+1,z]-self.Gv[2,x,y,z] )
                        )
            else:
                print('Stress Error BCs 1', x,y,z)
        elif BCs == 2:
            pass
        elif BCs == 3:
            pass
        else:
            print(x,y,z,'BCs Count')
        #Ds = self.CheckStressBoundary(x,y,z,Ds)
        Ds[2,0]=Ds[0,2]
        Ds[1,0]=Ds[0,1]      
        Ds[2,1]=Ds[1,2]

        return Ds

    def DeltaVelocity(self, x,y,z):
        #Gets the change in the velocity per time at a certain coordinate juncture 
        #
        # Inputs: x,y,z coordinates of the cube in question.  Last time is assumed
        #
        # Outputs: 3 dimmensions of velocity.  
        
        #Calculated velocity based on 3.54

        DV = np.zeros(3)

        BCs = 0
        if x == self.MaxX: BCs+=1
        if x == 0: BCs+=1
        if y == self.MaxY: BCs+=1
        if y == 0: BCs+=1
        if z == self.MaxZ: BCs+=1
        if z == 0: BCs+=1

        if BCs == 0:
                DV[0] = ((self.ids ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x+1,y,z])) *
                        (self.Gs[0,0,x+1,y,z] - self.Gs[0,0,x,y,z] + self.Gs[0,1,x,y,z] - self.Gs[0,1,x,y-1,z] + self.Gs[0,2,x,y,z] - self.Gs[0,2,x,y,z-1])
                        )
                DV[1] = ((self.ids ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x,y+1,z])) *
                        (self.Gs[0,1,x,y,z] - self.Gs[0,1,x-1,y,z] + self.Gs[1,1,x,y+1,z] - self.Gs[1,1,x,y,z] + self.Gs[1,2,x,y,z] - self.Gs[1,2,x,y,z-1])
                        )
                DV[2] = ((self.ids ) *
                        (2 / (self.Gp[0,x,y,z]+self.Gp[0,x,y,z+1])) *
                        (self.Gs[0,2,x,y,z] - self.Gs[0,2,x-1,y,z] + self.Gs[1,2,x,y,z] - self.Gs[1,2,x,y-1,z] + self.Gs[2,2,x,y,z+1] - self.Gs[2,2,x,y,z])
                        )
        elif BCs == 1:
            if x == 0:
                DV[0] = -2.0 * self.Gs[0,0,x+1,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[1] = -2.0 * self.Gs[1,1,x+1,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[2] = -2.0 * self.Gs[2,2,x+1,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
            elif y == 0:
                DV[0] = -2.0 * self.Gs[0,0,x,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[1] = -2.0 * self.Gs[1,1,x,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[2] = -2.0 * self.Gs[2,2,x,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
            elif z == 0:
                DV[0] = -2.0 * self.Gs[0,0,x,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[1] = -2.0 * self.Gs[1,1,x,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[2] = -2.0 * self.Gs[2,2,x,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
            elif x == self.MaxX or y == self.MaxY or z == self.MaxZ:
                DV[0] = 2.0 * self.Gs[0,0,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[1] = 2.0 * self.Gs[1,1,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[2] = 2.0 * self.Gs[2,2,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
            else:
                print('Error BCs 1',x,y,z)
        elif BCs == 2:
            if x == self.MaxX:
                if y == 0:
                    DV[0] = 2.0 * self.Gs[0,0,x,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[1] = 2.0 * self.Gs[1,1,x,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[2] = 2.0 * self.Gs[2,2,x,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                elif z == 0:
                    DV[0] = 2.0 * self.Gs[0,0,x,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[1] = 2.0 * self.Gs[1,1,x,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[2] = 2.0 * self.Gs[2,2,x,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                else:
                    DV[0] = 2.0 * self.Gs[0,0,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[1] = 2.0 * self.Gs[1,1,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[2] = 2.0 * self.Gs[2,2,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
            elif y == self.MaxY:
                if x == 0:
                    DV[0] = 2.0 * self.Gs[0,0,x+1,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[1] = 2.0 * self.Gs[1,1,x+1,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[2] = 2.0 * self.Gs[2,2,x+1,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                elif z == 0:
                    DV[0] = 2.0 * self.Gs[0,0,x,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[1] = 2.0 * self.Gs[1,1,x,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[2] = 2.0 * self.Gs[2,2,x,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                else:
                    DV[0] = 2.0 * self.Gs[0,0,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[1] = 2.0 * self.Gs[1,1,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[2] = 2.0 * self.Gs[2,2,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
            elif z == self.MaxZ:
                if y == 0:
                    DV[0] = 2.0 * self.Gs[0,0,x,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[1] = 2.0 * self.Gs[1,1,x,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[2] = 2.0 * self.Gs[2,2,x,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                elif x == 0:
                    DV[0] = 2.0 * self.Gs[0,0,x+1,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[1] = 2.0 * self.Gs[1,1,x+1,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[2] = 2.0 * self.Gs[2,2,x+1,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                else:
                    DV[0] = 2.0 * self.Gs[0,0,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[1] = 2.0 * self.Gs[1,1,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                    DV[2] = 2.0 * self.Gs[2,2,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
            elif x == 0 and y == 0:
                DV[0] = -2.0 * self.Gs[0,0,x+1,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[1] = -2.0 * self.Gs[1,1,x+1,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[2] = -2.0 * self.Gs[2,2,x+1,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
            elif y == 0 and z == 0:
                DV[0] = -2.0 * self.Gs[0,0,x,y+1,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[1] = -2.0 * self.Gs[1,1,x,y+1,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[2] = -2.0 * self.Gs[2,2,x,y+1,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
            elif x == 0 and z == 0:
                DV[0] = -2.0 * self.Gs[0,0,x+1,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[1] = -2.0 * self.Gs[1,1,x+1,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                DV[2] = -2.0 * self.Gs[2,2,x+1,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
            else:
                print('Error BCs 2',x,y,z)
        elif BCs == 3:
            if x == 0:
                if y == 0:
                    if z==0:
                        DV[0] = -2.0 * self.Gs[0,0,x+1,y+1,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[1] = -2.0 * self.Gs[1,1,x+1,y+1,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[2] = -2.0 * self.Gs[2,2,x+1,y+1,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                    else:
                        DV[0] = -2.0 * self.Gs[0,0,x+1,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[1] = -2.0 * self.Gs[1,1,x+1,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[2] = -2.0 * self.Gs[2,2,x+1,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                else:
                    if z == 0:
                        DV[0] = -2.0 * self.Gs[0,0,x+1,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[1] = -2.0 * self.Gs[1,1,x+1,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[2] = -2.0 * self.Gs[2,2,x+1,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                    else:
                        DV[0] = -2.0 * self.Gs[0,0,x+1,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[1] = -2.0 * self.Gs[1,1,x+1,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[2] = -2.0 * self.Gs[2,2,x+1,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
            else:
                if y == 0:
                    if z==0:
                        DV[0] = -2.0 * self.Gs[0,0,x,y+1,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[1] = -2.0 * self.Gs[1,1,x,y+1,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[2] = -2.0 * self.Gs[2,2,x,y+1,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                    else:
                        DV[0] = -2.0 * self.Gs[0,0,x,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[1] = -2.0 * self.Gs[1,1,x,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[2] = -2.0 * self.Gs[2,2,x,y+1,z] * (1/self.Gp[0,x,y,z]) * self.ids
                else:
                    if z == 0:
                        DV[0] = -2.0 * self.Gs[0,0,x,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[1] = -2.0 * self.Gs[1,1,x,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[2] = -2.0 * self.Gs[2,2,x,y,z+1] * (1/self.Gp[0,x,y,z]) * self.ids
                    else:
                        DV[0] = -2.0 * self.Gs[0,0,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[1] = -2.0 * self.Gs[1,1,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids
                        DV[2] = -2.0 * self.Gs[2,2,x,y,z] * (1/self.Gp[0,x,y,z]) * self.ids                                    
        else:
            DV[0]=0
            DV[1]=0
            DV[2]=0
            
        #DV = self.CheckVelocityBoundary(x,y,z,DV)

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
    
    def StepStresses(self):
        for x in range(self.MaxX+1):
            for y in range(self.MaxY+1):
                for z in range(self.MaxZ+1):
                    self.UpdateStresses(x,y,z)
        self.Gs[0,0,0,:,:] = -self.Gs[0,0,1,:,:]
        self.Gs[1,1,:,0,:] = -self.Gs[1,1,:,1,:]
        self.Gs[2,2,:,:,0] = -self.Gs[2,2,:,:,1]
        self.Gs[0,0,self.MaxX,:,:] = -self.Gs[0,0,self.MaxX-1,:,:]
        self.Gs[1,1,:,self.MaxY,:] = -self.Gs[1,1,:,self.MaxY-1,:]
        self.Gs[2,2,:,:,self.MaxZ] = -self.Gs[2,2,:,:,self.MaxZ-1]

    def StepVelocities(self):
        for x in range(self.MaxX+1):
            for y in range(self.MaxY+1):
                for z in range(self.MaxZ+1):
                    self.UpdateVelocity(x,y,z)
        self.Gv[:,:,:,0]=0

    def ForcingFunctionWave(self, t, Hz = 40000, EP=100.0, size=0.04, Odim = 2, Dir=1):
        # Adds stresses from a force to the stress grid
        # Initially assumed a single force of a small plate sinosoidal ultrasound emitter.  More to be added later
        # 
        # Input:   t is the time
        #          Hz is the frequency of the signal
        #          EP is the force / density
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
            if Dir ==1:
                self.Gv[0,self.MaxX-1,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
                self.Gv[0,self.MaxX,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
            else:
                self.Gv[0,0,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
                self.Gv[0,1,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
        elif Odim == 1:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxZ / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.Gv[1,Start0:Start1+EmitterWidth,self.MaxY-1,Start0:Start1+EmitterWidth] = Temp
                self.Gv[1,Start0:Start1+EmitterWidth,self.MaxY,Start0:Start1+EmitterWidth] = Temp
            else:
                self.Gv[2,Start0:Start1+EmitterWidth,1,Start0:Start1+EmitterWidth] = Temp
                self.Gv[2,Start0:Start1+EmitterWidth,0,Start0:Start1+EmitterWidth] = Temp
        else:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxY / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,self.MaxZ-1] = Temp
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,self.MaxZ] = Temp
            else:
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,1] = Temp
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,0] = Temp
        return self

    def ForcingFunctionImpulse(self, force = 100.0, emitter = 0.01, Odim = 2, Dir=1):
        # Adds stresses from a force to the stress grid
        # Initially assumed a single force of a small plate sinosoidal ultrasound emitter.  More to be added later
        # 
        # Input is the time, impulse, square emitter size, what face it is on, and if it is top or bottom
        #
        # Outputs: no direct outputs, last time step stress is updated

        EmitterWidth = emitter * self.ids
        EmitterWidth = int(EmitterWidth)
        if EmitterWidth == 0: EmitterWidth = 1

        
        Temp = np.zeros((EmitterWidth,EmitterWidth))

        Temp[:,:] = -force / 2.0

        if Odim == 0:
            Start0 = int((self.MaxY / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxZ / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.Gv[0,self.MaxX,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
                self.Gv[0,self.MaxX-1,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
            else:
                self.Gv[0,0,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
                self.Gv[0,1,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth] = Temp
        elif Odim == 1:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxZ / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.Gv[1,Start0:Start1+EmitterWidth,self.MaxY,Start0:Start1+EmitterWidth] = Temp
                self.Gv[1,Start0:Start1+EmitterWidth,self.MaxY-1,Start0:Start1+EmitterWidth] = Temp
            else:
                self.Gv[2,Start0:Start1+EmitterWidth,0,Start0:Start1+EmitterWidth] = Temp
                self.Gv[2,Start0:Start1+EmitterWidth,1,Start0:Start1+EmitterWidth] = Temp
        else:
            Start0 = int((self.MaxX / 2) - (EmitterWidth / 2))
            Start1 = int((self.MaxY / 2) - (EmitterWidth / 2))
            if Dir ==1:
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,self.MaxZ] = Temp
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,self.MaxZ-1] = Temp
            else:
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,1] = Temp
                self.Gv[2,Start0:Start1+EmitterWidth,Start0:Start1+EmitterWidth,0] = Temp


        return self

    def VelocityCut(self, Dimm, Component = -1, location = -1):
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
        if Dimm == -1:
            Component0 = self.Gv[0,:,:,:]
            Component1 = self.Gv[1,:,:,:]
            Component2 = self.Gv[2,:,:,:]
            
            #Results = Component1 #+ Component2
            Results = np.sqrt(Component0**2+Component1**2+Component2**2)
        elif Dimm == 0 or Dimm == 1 or Dimm == 2:
            Results = self.Gv[Dimm,:,:,:]
        else:
            print('Error, unknown Dimm', Dimm)
        
        return Results
