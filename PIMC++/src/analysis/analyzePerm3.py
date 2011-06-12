#!/usr/bin/python
import tables
import povexport
from visual import *
import stats
import sys
import Numeric

BoxSize=8.5
def PutInBox(v):
	for dim in range(0,len(v)):
		while(v[dim]>0.5*BoxSize):
			v[dim]=v[dim]-BoxSize
		while (v[dim]<-0.5*BoxSize):
			v[dim]=v[dim]+BoxSize
	return v


NDIM=3
def AllInBox(singlePath):
	singlePath[0]=PutInBox(singlePath[0])
	for dim in range(0,NDIM):
		for counter in range(1,len(singlePath)):
			mNew=-floor((singlePath[counter-1][dim]-singlePath[counter][dim])*1.0/BoxSize+0.5)
#			print singlePath[counter][dim]
#			print "hi"
#			print mNew*BoxSize
			singlePath[counter][dim]=singlePath[counter][dim]-mNew*BoxSize
	return singlePath
	
#h5file = tables.openFile("../TestPerm.h5")
h5file = tables.openFile(sys.argv[1])
data = h5file.root.PathDump_1.Path.read()

for sliceCounter in range(0,len(data)):
	for ptclNumber in range(0,len(data[0])):
		data[sliceCounter][ptclNumber]=AllInBox(data[sliceCounter][ptclNumber])
	
path1=[]
path2=[]
path3=[]
path4=[]
path5=[]
path6=[]
path7=[]
path8=[]

NumImages=0
for nx  in range(-NumImages,NumImages+1):
	for ny in range(-NumImages,NumImages+1):
		for nz in range(-NumImages,NumImages+1):
			path1.append(curve(color=color.yellow,radius=.1))
			path2.append(curve(color=color.blue,radius=.1))
			path3.append(curve(color=color.red,radius=.1))
			path4.append(curve(color=color.orange,radius=.1))
			path5.append(curve(color=color.green,radius=.1))
			path6.append(curve(color=color.cyan,radius=.1))
			path7.append(curve(color=color.magenta,radius=.1))
			path8.append(curve(color=color.white,radius=.1))

#timeball1=sphere(pos=data[0][0][0],radius=.1,color=color.green)
#timeball2=sphere(pos=data[0][0][0],radius=.1,color=color.orange)

ball1 = sphere (pos=data[0,0,0], radius=0.05, color=color.red)
ball2 = sphere (pos=data[0,1,0], radius=0.05, color=color.red)
ball3 = sphere (pos=data[0,2,0], radius=0.05, color=color.red)
ball4 = sphere (pos=data[0,3,0], radius=0.05, color=color.red)
ball5 = sphere (pos=data[0,4,0], radius=0.05, color=color.red)
ball6 = sphere (pos=data[0,5,0], radius=0.05, color=color.red)
ball7 = sphere (pos=data[0,6,0], radius=0.05, color=color.red)
ball8 = sphere (pos=data[0,7,0], radius=0.05, color=color.red)

mybox = box(pos=(0,0,0), length=8.5,
height=8.5, width=8.5)
mybox.opacity=0.13
scene.autoscale=1
pathNum=0
done=False 
ballslice=0
while (not(done)==True):
       if scene.mouse.clicked:
         m = scene.mouse.getclick()
         loc = m.pos
         print loc
       scene.autocenter
       s=scene.kb.getkey()
       if s=="left":
              pathNum=pathNum-1
              if pathNum<0:
                     pathNum=0
       elif s=="q":
              print "We are done now"
              done=True
       elif s=="right":
              pathNum=pathNum+1
              if pathNum>=len(data):
                     pathNum=len(data)-1
       elif s=="r":
       	     ballslice=ballslice+1
	     if ballslice>=len(data[pathNum][0]):
	     	ballslice=len(data[pathNum][0])-1
       elif s=="l":
       	     ballslice=ballslice-1
	     if ballslice<0:
	     	ballslice=0
       elif s=="b":
             mybox.visible=0
       elif s=="B":
             mybox.visible=1
       elif s=="C":
             scene.center=loc
       elif s=="p":
	       povexport.export(filename="export.pov")
       
       print ballslice
       print len(data[pathNum][0])
       ball1.visible=0
       ball2.visible=0
       ball3.visible=0
       ball4.visible=0
       ball5.visible=0
       ball6.visible=0
       ball7.visible=0
       ball8.visible=0
       ball1 = sphere (pos=(data[pathNum,0,ballslice]), radius=0.1, color=color.red)
       ball2 = sphere (pos=(data[pathNum,1,ballslice]), radius=0.1, color=color.red)
       ball3 = sphere (pos=(data[pathNum,2,ballslice]), radius=0.1, color=color.red)
       ball4 = sphere (pos=(data[pathNum,3,ballslice]), radius=0.1, color=color.red)
       ball5 = sphere (pos=(data[pathNum,4,ballslice]), radius=0.1, color=color.red)
       ball6 = sphere (pos=(data[pathNum,5,ballslice]), radius=0.1, color=color.red)
       ball7 = sphere (pos=(data[pathNum,6,ballslice]), radius=0.1, color=color.red)
       ball8 = sphere (pos=(data[pathNum,7,ballslice]), radius=0.1, color=color.red)
###       ball3 = sphere (pos=data[pathNum,2,0], radius=0.05, color=color.red)
       for nx  in range(0,2*NumImages+1):
	       for ny in range(0,2*NumImages+1):
		       for nz in range(0,2*NumImages+1):       
			       nzImage=nz-NumImages
			       nyImage=ny-NumImages
			       nxImage=nx-NumImages
			       if (1==1): #(nzImage==0 and nyImage==0) or (nzImage==0 and nxImage==0) or (nxImage==0 and nyImage==0):
				       path1[nx*9+ny*3+nz].visible=0
				       path2[nx*9+ny*3+nz].visible=0
				       path3[nx*9+ny*3+nz].visible=0
				       path4[nx*9+ny*3+nz].visible=0
				       path5[nx*9+ny*3+nz].visible=0
				       path6[nx*9+ny*3+nz].visible=0
				       path7[nx*9+ny*3+nz].visible=0
				       path8[nx*9+ny*3+nz].visible=0
				       path1[nx*9+ny*3+nz]=curve(color=color.yellow,radius=.05)
				       path2[nx*9+ny*3+nz]=curve(color=color.blue,radius=.05)
				       path3[nx*9+ny*3+nz]=curve(color=color.red,radius=.05)
				       path4[nx*9+ny*3+nz]=curve(color=color.orange,radius=.05)
				       path5[nx*9+ny*3+nz]=curve(color=color.green,radius=.05)
				       path6[nx*9+ny*3+nz]=curve(color=color.cyan,radius=.05)
				       path7[nx*9+ny*3+nz]=curve(color=color.magenta,radius=.05)
				       path8[nx*9+ny*3+nz]=curve(color=color.white,radius=.05)
##     timeball1=sphere(pos=data[pathNum][0][0],radius=.1,color=color.green)
##     timeball2=sphere(pos=data[pathNum][1][0],radius=.1,color=color.orange)
       for nx  in range(0,2*NumImages+1):
	       for ny in range(0,2*NumImages+1):
		       for nz in range(0,2*NumImages+1):       
			       for slice in range(0,len(data[pathNum][0])):
				       nzI=nz-NumImages
				       nyI=ny-NumImages
				       nxI=nx-NumImages
				       if (1==1): #(nzI==0 and nyI==0) or (nzI==0 and nxI==0) or (nxI==0 and nyI==0):
					       path1[nx*9+ny*3+nz].append(pos=(data[pathNum][0][slice]+[(nxI)*BoxSize,(nyI)*BoxSize,(nzI)*BoxSize]))
					       path2[nx*9+ny*3+nz].append(pos=(data[pathNum][1][slice]+[(nxI)*BoxSize,(nyI)*BoxSize,(nzI)*BoxSize]))
					       path3[nx*9+ny*3+nz].append(pos=(data[pathNum][2][slice]+[(nxI)*BoxSize,(nyI)*BoxSize,(nzI)*BoxSize]))
					       path4[nx*9+ny*3+nz].append(pos=(data[pathNum][3][slice]+[(nxI)*BoxSize,(nyI)*BoxSize,(nzI)*BoxSize]))
					       path5[nx*9+ny*3+nz].append(pos=(data[pathNum][4][slice]+[(nxI)*BoxSize,(nyI)*BoxSize,(nzI)*BoxSize]))
					       path6[nx*9+ny*3+nz].append(pos=(data[pathNum][5][slice]+[(nxI)*BoxSize,(nyI)*BoxSize,(nzI)*BoxSize]))
					       path7[nx*9+ny*3+nz].append(pos=(data[pathNum][6][slice]+[(nxI)*BoxSize,(nyI)*BoxSize,(nzI)*BoxSize]))
					       path8[nx*9+ny*3+nz].append(pos=(data[pathNum][7][slice]+[(nxI)*BoxSize,(nyI)*BoxSize,(nzI)*BoxSize]))
	       
###              path3.append(pos=data[pathNum][2][slice])
       scene.autoscale=1
            



