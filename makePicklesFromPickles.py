# %% [markdown]
# # Make Animation from pickels files

# %%
#Standard Header used on the projects

#first the major packages used for math and graphing
import numpy as np
import matplotlib.pyplot as plt
from cycler import cycler
import platform
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

# %%
import multiprocessing
from joblib import Parallel, delayed
num_jobs=30

# %%
## Task specific imports
import os as os
#import keras.utils as image
import glob
from PIL import Image, ImageDraw
import pickle

# %%
import math
import time

# %% [markdown]
# ## Choosing Platform
# Working is beinging conducted on several computers, and author needs to be able to run code on all without rewriting..  This segment of determines which computer is being used, and sets the directories accordingly.

# %%
HostName = platform.node()

if HostName == "Server":
    Computer = "Desktop"   
elif HostName[-6:] == 'wm.edu':
    Computer = "SciClone"
elif HostName == "SchoolLaptop":
    Computer = "LinLap"
elif HostName == "WTC-TAB-512":
    Computer = "PortLap"
else:
    Computer = "WinLap"
    
if Computer ==  "SciClone":
    rootfolder = '/sciclone/scr10/dchendrickson01/EFIT/'
    
else:
    asdfasdf

# %% [markdown]
# ## Load data

# %%
Case = '20m15XRfRq'
CasesAtaTime = 2
FilesAtTime = 5
ProcessPerFile = 5
Position = 0 #0 Energy, 1 X axis, 2 Y Axis, 3 X axis
skips = 2 #1 every data point, or skip what modulous

imFolder=rootfolder+Case+'/'
fileNames = glob.glob(imFolder+'Movements*.p')

Views=[]

# %%
fileName = imFolder+'Parameters.p'

file=open(fileName,'rb')
Parameters=pickle.load(file)

file.close()


# %%
xmax = Parameters["LargestXnode"]
ymax = Parameters["LargestYnode"]
zmax = Parameters["LargestZnode"]
gridStartWeb = Parameters["WidthStartWebNode"]
gridEndWeb = Parameters["WidthEndWebNode"]
gridEndFoot = Parameters["HeightEndFootNode"]
gridStartHead = Parameters["HeightStartHeadNode"]
gridStartHeadWidth = Parameters["WidthStartHeadNode"]
gridEndHeadWidth = Parameters["WidthEndHeadNode"]
DataBucketSize = Parameters["SaveEveryXStep"]

# %%
Parameters

# %% [markdown]
# ## Image making functions

# %%
def getFigData(fileName, Position, skips, xStart, xEnd, yStart, yEnd, zStart, zEnd):
    
    file=open(fileName,'rb')
    temp = pickle.load(file)
    file.close()
    
    Data = temp[Position][:,:,:,0::skips]
    
    del temp


    if xStart == xEnd:
        ReturnData = Data[xStart,:,:,:]
        ReturnData[:gridStartWeb-1,gridEndFoot+2:gridStartHead-2] = np.nan
        ReturnData[:gridStartHeadWidth-1,gridStartHead-2:] = np.nan
        ReturnData[gridEndWeb+2:,gridEndFoot+2:gridStartHead-1] = np.nan
        ReturnData[gridEndHeadWidth+4:,gridStartHead-1:] = np.nan
        
    elif yStart == yEnd:
        ReturnData = Data[:,yStart,:,:]
    elif zStart == zEnd:
        ReturnData = Data[:,:,zStart,:]
    else:
        print('Error no dimmension is a plane')
        ReturnData = []
        
        
    return ReturnData, int(fileName[-5:][:3])

# %%
def SimpleFig(Data, t, v, figW, figH, Folder):
    
    fig = plt.figure(figsize=(figH,figW), dpi=300)

    Data = Data.T
    
    #print(Data.shape)
    
    plt.contourf(Data, v, cmap=plt.cm.jet)
    plt.title(Folder+' '+str(t))
    plt.savefig(imFolder+Folder+'/DataFrame'+str(t).zfill(5)+'.png')
    
    fig.clf()
    plt.close(fig)
    del fig

    return 0

# %%
def CaseImage(View):
    imfolder = rootfolder + Case + '/'+View+'/'
    files = glob.glob(os.path.join(imfolder, '*.png'))
    files.sort()
    images = []

    for file in files:
        images.append(image.load_img(file))
    images[0].save(rootfolder + Case + '/Animated_'+View+'.gif',
        save_all=True, append_images=images[1:], optimize=False, duration=40, loop=0)
    return 0

# %%
def runCase(Case, skips):
    
    Position = Case[7]
    xStart = Case[0]
    xEnd = Case[1]
    yStart = Case[2]
    yEnd = Case[3]
    zStart = Case[4]
    zEnd = Case[5]
    Folder = Case[6]
    
    temp = Parallel(n_jobs=FilesAtTime)(delayed(getFigData)
                                        (fileName, Position, skips, xStart, xEnd, yStart, yEnd, zStart, zEnd) 
                                        for fileName in fileNames[:5])
    
    StackData = np.zeros((temp[0][0].shape[0],temp[0][0].shape[1],1))
    
    for group in temp:
        StackData = np.concatenate((StackData,group[0]),axis=2)
    del temp
    
    EMin = np.min(StackData)
    EMax = np.max(StackData)
    v = np.linspace(EMin, EMax, 30, endpoint=True)[0:20]
    
    if os.path.isdir(imFolder+Folder):
        pass
    else:
        os.makedirs(imFolder+Folder)
        
    figW = 6
    figH = figW * (StackData.shape[0]/StackData.shape[1])
    
    #print(StackData[:,:,1].T.shape,StackData[:,:,17].T.shape,StackData[:,:,99].T.shape)
    
    #Parallel(n_jobs=ProcessPerFile)(delayed(SimpleFig)
    #                                (StackData[:,:,i], i, v, figW, figH, Folder)
    #                                for i in range(StackData.shape[2]))
    for i in range(StackData.shape[2]):
        a = SimpleFig(StackData[:,:,i], i, v, figW, figH, Folder)
    
    del StackData
    #temp = CaseImage(Folder)

# %%
def getStackData(Case, skips):
    
    xStart = Case['xstart']
    xEnd = Case['xend']
    yStart = Case['ystart']
    yEnd = Case['yend']
    zStart = Case['zstart']
    zEnd = Case['zend']
    Folder = Case['name']
    Position = Case['Position']
    
    temp = Parallel(n_jobs=FilesAtTime)(delayed(getFigData)
                                        (fileName, Position, skips, xStart, xEnd, yStart, yEnd, zStart, zEnd) 
                                        for fileName in fileNames)
    
    StackData = np.zeros((temp[0][0].shape[0],temp[0][0].shape[1],1))
    
    for group in temp:
        StackData = np.concatenate((StackData,group[0]),axis=2)
    del temp
    
    return StackData

# %% [markdown]
# ## make all the frames of all the cases

# %%
DataCases = {
             'EndM4': {
                       'name'     : 'EndM4',
                       'xstart'   : xmax - 4,
                       'xend'     : xmax - 4,
                       'ystart'   : 0,
                       'yend'     : ymax,
                       'zstart'   : 0,
                       'zend'     : zmax,
                       'Position' : 0
                      },
             'EndM4y': {
                       'name'     : 'EndM4y',
                       'xstart'   : xmax - 4,
                       'xend'     : xmax - 4,
                       'ystart'   : 0,
                       'yend'     : ymax,
                       'zstart'   : 0,
                       'zend'     : zmax,
                       'Position' : 2
                      },
             'EndM2': {
                       'name'     : 'EndM2',
                       'xstart'   : xmax - 2,
                       'xend'     : xmax - 2,
                       'ystart'   : 0,
                       'yend'     : ymax,
                       'zstart'   : 0,
                       'zend'     : zmax,
                       'Position' : 0
                      },
             'EndM2y': {
                       'name'     : 'EndM2y',
                       'xstart'   : xmax - 2,
                       'xend'     : xmax - 2,
                       'ystart'   : 0,
                       'yend'     : ymax,
                       'zstart'   : 0,
                       'zend'     : zmax,
                       'Position' : 2
                      },
             'EndM2x': {
                       'name'     : 'EndM2x',
                       'xstart'   : xmax - 2,
                       'xend'     : xmax - 2,
                       'ystart'   : 0,
                       'yend'     : ymax,
                       'zstart'   : 0,
                       'zend'     : zmax,
                       'Position' : 1
                      },
             'WebStart': {
                       'name'     : 'WebStart',
                       'xstart'   : 0,
                       'xend'     : xmax ,
                       'ystart'   : gridStartWeb+1,
                       'yend'     : gridStartWeb+1,
                       'zstart'   : gridEndFoot,
                       'zend'     : gridStartHead,
                       'Position' : 1
                      },
             'WebEnd': {
                       'name'     : 'WebEnd',
                       'xstart'   : 0,
                       'xend'     : xmax ,
                       'ystart'   : gridEndWeb-1,
                       'yend'     : gridEndWeb-1,
                       'zstart'   : gridEndFoot,
                       'zend'     : gridStartHead,
                       'Position' : 1
                      },
             'TopSurface': {
                       'name'     : 'TopSurface',
                       'xstart'   : 0,
                       'xend'     : xmax ,
                       'ystart'   : gridStartHeadWidth-1,
                       'yend'     : gridEndHeadWidth+1,
                       'zstart'   : zmax-3,
                       'zend'     : zmax-3,
                       'Position' : 0
                      },
             'HeadEnd': {
                       'name'     : 'HeadEnd',
                       'xstart'   : 0,
                       'xend'     : xmax,
                       'ystart'   : gridEndHeadWidth-1,
                       'yend'     : gridEndHeadWidth-1,
                       'zstart'   : gridStartHead,
                       'zend'     : zmax,
                       'Position' : 2
                      },
             'LeftHead': {
                       'name'     : 'LeftHead',
                       'xstart'   : 0,
                       'xend'     : xmax ,
                       'ystart'   : gridStartHeadWidth,
                       'yend'     : gridStartHeadWidth,
                       'zstart'   : gridStartHead,
                       'zend'     : zmax,
                       'Position' : 1
                      },
             'RightHead': {
                       'name'     : 'RightHead',
                       'xstart'   : 0,
                       'xend'     : xmax ,
                       'ystart'   : gridEndHeadWidth,
                       'yend'     : gridEndHeadWidth,
                       'zstart'   : gridStartHead,
                       'zend'     : zmax,
                       'Position' : 1
                      },
             'EndM6': {
                       'xstart'   : xmax - 6,
                       'xend'     : xmax - 6,
                       'ystart'   : 0,
                       'yend'     : ymax,
                       'zstart'   : 0,
                       'zend'     : zmax,
                       'Position' : 0
                      },
             'HeadStart': {
                       'xstart'   : 0,
                       'xend'     : xmax ,
                       'ystart'   : gridStartHeadWidth+1,
                       'yend'     : gridStartHeadWidth+1,
                       'zstart'   : gridStartHead,
                       'zend'     : zmax,
                       'Position' : 0
                      },
             'MiddleHeadPlane': {
                       'xstart'   : 0,
                       'xend'     : xmax ,
                       'ystart'   : 0,
                       'yend'     : ymax,
                       'zstart'   : int((zmax - gridStartHead)/2)+gridStartHead,
                       'zend'     : int((zmax - gridStartHead)/2)+gridStartHead,
                       'Position' : 1
                      },
             'MiddleWebPlane': {
                       'xstart'   : 0,
                       'xend'     : xmax ,
                       'ystart'   : gridStartWeb,
                       'yend'     : gridEndWeb,
                       'zstart'   : int(zmax/2),
                       'zend'     : int(zmax/2),
                       'Position' : 1
                      },
             'MiddleVerticalPlane': {
                       'xstart'   : 0,
                       'xend'     : xmax ,
                       'ystart'   : int(ymax/2),
                       'yend'     : int(ymax/2),
                       'zstart'   : 0,
                       'zend'     : zmax,
                       'Position' : 1
                      }
    
}

# %%
stime = time.time()

for entry in DataCases:
    job = entry['name']

    Data = getStackData(DataCases[job],15)

    file=open(imFolder+'Data-'+job+'.p','wb')
    pickle.dump(Data,file)
    file.close()

    print(int((time.time()-stime)/60.0*100)/100)

