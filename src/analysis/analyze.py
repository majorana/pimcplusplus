#!/usr/bin/python
import tables
from visual import *
#import Numeric
h5file = tables.openFile("../H2_bservables.h5");
data = h5file.root.Paths_2.Path.read()
names=h5file.root.Paths_2.SpeciesNames.read()

E = h5file.root.Energies_0.TotalEnergy.read()
Emean = sum(E)/len(E)
print Emean


path1=curve(color=color.yellow,radius=5)
path2=curve(color=color.blue,radius=5)
for pathNum in range(0,len(data)):
    rate(1)
    ball1 = sphere (pos=data[pathNum,2,0], radius=0.2, color=color.red)
    ball2 = sphere (pos=data[pathNum,3,0], radius=0.2, color=color.red)
    path1.visible=0
    path2.visible=0
    path1=curve(color=color.yellow,radius=0.0)
    path2=curve(color=color.blue,radius=0.0)
    for slice in range(0,len(data[pathNum][0])):
        path1.append(pos=data[pathNum][0][slice])
        path2.append(pos=data[pathNum][1][slice])
    scene.autoscale=0
            
                


