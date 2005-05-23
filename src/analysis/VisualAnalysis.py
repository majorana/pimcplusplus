from IO import *
import sys
import os
import stats
import numarray
from matplotlib.matlab import *
from HTMLgen import *
import povexport
from visual import *
from visual.controls import *

def GetPaths(infile):
     paths=infile.ReadVar("Path")
     return paths

def GetOpen(infile):
     theOpen=infile.ReadVar("OpenPtcl")
     theSlice=infile.ReadVar("OpenLinkSlice")
     theTail=infile.ReadVar("TailLocation")
     return theOpen

def InitVisualPaths(pathData):
     visualPath=[]
     visualBall=[]
#     print "Path Data is ",pathData
     numPaths=len(pathData[0])
     for pathNum in range(0,numPaths):
          visualPath.append(curve(color=color.green,radius=0.02))
          visualBall.append(sphere(pos=pathData[0,0,0], radius=0.05, color=color.red))
     return (visualPath,visualBall)




def ErasePaths(visualPath,visualBall):
     for pathNum in range(0,len(visualPath)):
          visualBall[pathNum].visible=0
          visualPath[pathNum].visible=0


box=zeros(3)+0.0
box[0]=24.0
box[1]=24.0
box[2]=24.0

def inBox(pathData,mcTime,pathNum,slice):
    myArray=zeros(len(pathData[mcTime][pathNum][slice]))+0.0
    for counter in range(0,len(pathData[mcTime][pathNum][slice])):
        myArray[counter]=pathData[mcTime][pathNum][slice][counter] % box[counter]
    return myArray

def Id(pathData,mcTime,pathNum,slice):
    return pathData[mcTime][pathNum][slice]

def PlotPaths(pathData,visualPath,visualBall,mcTime,ballTime):
     maxMCTime=len(pathData)
     if mcTime>=maxMCTime:
          mcTime=maxMCTime-1
     numPaths=len(pathData[mcTime])
     for pathNum in range(0,numPaths):
          if openData[mcTime]!=pathNum:
               visualPath[pathNum]=(curve(color=color.green,radius=0.02))
          else:
               visualPath[pathNum]=(curve(color=color.red,radius=0.02))
          visualBall[pathNum]=sphere(pos=(Id(pathData,mcTime,pathNum,ballTime)), radius=0.1, color=color.red)
#          visualBall[pathNum]=sphere(pos=(pathData[mcTime,pathNum,ballTime]), radius=2.1, color=color.red)
          numSlices=len(pathData[mcTime][pathNum])
          for slice in range(0,numSlices):
              visualPath[pathNum].append(pos=(Id(pathData,mcTime,pathNum,slice)))
#              visualPath[pathNum].append(pos=(pathData[mcTime][pathNum][slice]))




scene.exit=1
infile=IOSectionClass()
infile.OpenFile(sys.argv[1])
fileString=sys.argv[1]
dotLoc=string.rfind(fileString,'.')
dirName=fileString[0:dotLoc]
if not(os.access(dirName,os.F_OK)):
     os.mkdir(dirName)
os.chdir(dirName)
infile.OpenSection("Observables")
infile.OpenSection("PathDump")
pathData=GetPaths(infile)
openData=GetOpen(infile)
print openData
infile.CloseSection()
infile.CloseSection()
##c = controls() # Create controls window

# Create a button in the controls window:

##b = button( pos=(0,0), width=60, height=60, text='Quit', action=lambda: quit() )

def quit():
    print "I quit!"
    bedone=true
    
(visualPath,visualBall)=InitVisualPaths(pathData)


maxmcTime=len(pathData)-1

maxBallTime=len(pathData[0][0])-1
mcTime=0
ballTime=0
PlotPaths(pathData,visualPath,visualBall,mcTime,ballTime)
bedone=false


    

#    c.interact()
#    if (scene.kb.keys):
while (not(bedone)):
    if (scene.kb.keys):
        s=scene.kb.getkey()
        if s=="left":
            ErasePaths(visualPath,visualBall)
            mcTime=mcTime-1
            if mcTime<0:
                mcTime=0
        elif s=="right":
            ErasePaths(visualPath,visualBall)
            mcTime=mcTime+1
            if mcTime>maxmcTime:
                mcTime=maxmcTime
        elif s=="r":
            ErasePaths(visualPath,visualBall)
            ballTime=ballTime+1
            if ballTime>maxBallTime:
                ballTime=maxBallTime
        elif s=="l":
            ErasePaths(visualPath,visualBall)
            ballTime=ballTime-1
            if ballTime<0:
                ballTime=0
        elif s=="q":
            bedone=true
        print "The current mctime is ",mcTime
        PlotPaths(pathData,visualPath,visualBall,mcTime,ballTime)
scene.visible=false

## elif s=="r":
##         ballslice=ballslice+1
##         if ballslice>=len(data[pathNum][0]):
##             ballslice=len(data[pathNum][0])-1
##     elif s=="l":
##         ballslice=ballslice-1
##         if ballslice<0:
##             ballslice=0
##     elif s=="b":
##         mybox.visible=0
##     elif s=="B":
##         mybox.visible=1
##     elif s=="C":
##         scene.center=loc
##     elif s=="p":
##         povexport.export(filename="export.pov")
## elif s=="q":
##         print "We are done now"
##         done=True
       

 

