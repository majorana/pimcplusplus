from IO import *
import sys
import os
import stats
import numarray
from matplotlib.matlab import *
from HTMLgen import *
import povexport
from visual import *


def GetPaths(infile):
     paths=infile.ReadVar("Path")
     return paths

def InitVisualPaths(pathData):
     visualPath=[]
     visualBall=[]
     numPaths=len(pathData[0])
     for pathNum in range(0,numPaths):
          visualPath.append(curve(color=color.blue,radius=0.2))
          visualBall.append(sphere(pos=pathData[0,0,0], radius=0.05, color=color.red))
     return (visualPath,visualBall)




def ErasePaths(pathData,visualPath,visualBall,mcTime):
     maxMCTime=len(pathData)
     if mcTime>=maxMCTime:
          mcTime=maxMCTime-1
     numPaths=len(pathData[mcTime])
     for pathNum in range(0,numPaths):
          visualBall[pathNum].visible=0
          visualPath[pathNum].visible=0



def PlotPaths(pathData,visualPath,visualBall,mcTime):
     maxMCTime=len(pathData)
     if mcTime>=maxMCTime:
          mcTime=maxMCTime-1
     numPaths=len(pathData[mcTime])
     for pathNum in range(0,numPaths):
          visualPath[pathNum]=(curve(color=color.blue,radius=2.2))
          visualBall[pathNum].visible=0
          visualBall[pathNum]=sphere(pos=(pathData[mcTime,pathNum,0]), radius=2.1, color=color.red)
          visualBall[pathNum].visible=1
#          print pathNum,numPaths
          numSlices=len(pathData[mcTime][pathNum])
          for slice in range(0,numSlices):
              visualPath[pathNum].append(pos=(pathData[mcTime][pathNum][slice]))
              visualPath[pathNum].visible=1




infile=IOSectionClass()
infile.OpenFile(sys.argv[1])
fileString=sys.argv[1]
dotLoc=string.rfind(fileString,'.')
dirName=fileString[0:dotLoc]
if not(os.access(dirName,os.F_OK)):
     os.mkdir(dirName)
os.chdir(dirName)
infile.OpenSection("PathDump")
pathData=GetPaths(infile)
infile.CloseSection()
(visualPath,visualBall)=InitVisualPaths(pathData)
maxmcTime=len(pathData)-1



mcTime=0
PlotPaths(pathData,visualPath,visualBall,mcTime)
while (1==1):
    s=scene.kb.getkey()
    if s=="left":
        ErasePaths(pathData,visualPath,visualBall,mcTime)
        mcTime=mcTime-1
        if mcTime<0:
            mcTime=0
    elif s=="right":
        ErasePaths(pathData,visualPath,visualBall,mcTime)
        mcTime=mcTime+1
        if mcTime>maxmcTime:
            mcTime=maxmcTime
    elif s=="q":
        exit()
    print "The current mctime is ",mcTime
    PlotPaths(pathData,visualPath,visualBall,mcTime)

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
       

 

