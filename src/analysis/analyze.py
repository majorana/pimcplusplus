#!/usr/bin/python
import tables
from visual import *
#import Numeric
h5file = tables.openFile("../Observables2.h5");
data = h5file.root.Paths_2.Path.read()
names=h5file.root.Paths_2.SpeciesNames.read()

path1=curve(color=color.yellow,radius=0.5)
path2=curve(color=color.blue,radius=0.5)
for pathNum in range(0,len(data)):
    rate(1)
    ball1 = sphere (pos=data[pathNum,2,0], radius=0.1, color=color.red)
    ball2 = sphere (pos=data[pathNum,3,0], radius=0.1, color=color.red)
    path1.visible=0
    path2.visible=0
    path1=curve(color=color.yellow)
    path2=curve(color=color.blue)
    for slice in range(0,len(data[pathNum][0])):
        path1.append(pos=data[pathNum][0][slice],radius=0.1)
        path2.append(pos=data[pathNum][1][slice],radius=0.1)
    scene.autoscale=0
            
                


