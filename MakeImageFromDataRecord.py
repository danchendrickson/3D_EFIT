import pandas as pd
import dask.dataframe as dd
import numpy as np

#
folder = '/scratch/EFIT/Double10m/'
file = 'Anima.csv'

# this is needed to call dask.compute
dataset= dd.read_csv(folder+file)

# for example take a subset
dataset = dataset[dataset['time']==12000]

maxDim = np.max([dataset['x'].max().compute(), dataset['y'].max().compute(), dataset['z'].max().compute()])
maxDim+=1
maxDim = int(maxDim)

data = np.empty((maxDim,int(maxDim/2),int(maxDim/2)))
color = np.empty((maxDim,int(maxDim/2),int(maxDim/2),4))

print('have data')

Max = dataset['DisX'].max().compute()
Min = dataset['DisX'].min().compute()

dataset['DisX'] += Min
dataset['DisX'] += 1


buckets = 15
bins = np.linspace(1,Min+Max+1,buckets+1)

import matplotlib.cm as colors

cmap = colors.get_cmap('rainbow')

colorBin = np.linspace(0,1,buckets)
Colors = []

for cb in colorBin:
    Colors.append(cmap(cb))

#print(rgba) # (0.99807766255210428, 0.99923106502084169, 0.74602077638401709, 1.0)

for part in dataset.partitions:
    for item in part.iterrows():
        data[item[1][1],item[1][2],item[1][3]] = item[1][5]
        color[item[1][1],item[1][2],item[1][3],:] = Colors[min(range(len(colorBin)), key=lambda i: abs(colorBin[i]-item[1][5]))]

print('image data made')

import matplotlib.pyplot as plt

# Define dimensions
Nx, Ny, Nz = maxDim, int(maxDim/2), int(maxDim/2)
#X, Y, Z = np.meshgrid(np.arange(Nx), np.arange(Ny), -np.arange(Nz))

# Create a figure with 3D ax
fig = plt.figure(figsize=(12,12), dpi=1000)
ax = fig.add_subplot(111, projection='3d')

print('fig set up')

data = data !=0

print('data made true/false')

ax.voxels(data,
          facecolors=color,
          linewidth=0.01)
#ax.set(xlabel='r', ylabel='g', zlabel='b')
#ax.set_aspect('equal')
# 

print('Axis set up')

# Show Figure
plt.savefig(folder+'singleimage2.png')
plt.show()
