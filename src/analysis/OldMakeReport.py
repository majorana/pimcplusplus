#!/usr/bin/python
from IO import *
import sys
import os
import stats
import numarray
from matplotlib.matlab import *
from HTMLgen import *
#import povexport
#from visual import *


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


def PlotPaths(pathData,visualPath,visualBall,mcTime):
     maxMCTime=len(pathData)
     if mcTime>=maxMCTime:
          mcTime=maxMCTime-1
     numPaths=len(pathData[mcTime])
     for pathNum in range(0,numPaths):
          visualPath[pathNum]=(curve(color=color.blue,radius=0.2))
          visualBall[pathNum].visible=0
          visualBall[pathNum]=sphere(pos=(pathData[0,pathNum,0]), radius=0.1, color=color.red)
          visualBall[pathNum].visible=1
          print pathNum,numPaths
          numSlices=len(pathData[mcTime][pathNum])
          for slice in range(0,numSlices):
               print pathNum,numPaths,numSlices,slice
               visualPath[pathNum].append(pos=(pathData[mcTime][pathNum][slice]))
     

def IsMonotonic (x):
     isMono = True
     for i in range(0,x.size()-2):
          isMono = isMono and (x[i+1] > x[i])
     return isMono
     

def ProduceCorrelationPicture(x,y,fileBase,hlabel,vlabel):
     clf()
#     infile.OpenSection("grid")
#     x=infile.ReadVar("Points")
#     infile.CloseSection()
     if (IsMonotonic(x)):
          plot(x, y)
     else:
          plot(x, y, 'o')
     h1=xlabel(hlabel)
     set(h1,"FontSize",20)
     v1=ylabel(vlabel)
     set(v1,"FontSize",20)
     labels = get(gca(), 'xticklabels')
     set(labels, 'fontsize', 16)
     labels = get(gca(), 'yticklabels')
     set(labels, 'fontsize', 16)
     savefig(fileBase+".png",dpi=60)
     savefig(fileBase+".ps")
     myImg=Image(fileBase+".png")
     return myImg


def ProcessCorrelationSection(infile,doc,currNum):
     sectionName=infile.GetName()
     doc.append(Heading(1,sectionName))
     hlabel=infile.ReadVar("xlabel")
     vlabel=infile.ReadVar("ylabel")
     data=infile.ReadVar("y")
     x=infile.ReadVar("x")
     description=infile.ReadVar("Description")
     doc.append(Heading(4,description))
     currNum=currNum+1
     baseName=sectionName+repr(currNum)
     myImg=ProduceCorrelationPicture(x, data[-1],baseName,hlabel,vlabel)
####     myImg=ProduceCorrelationPicture(x, data[-1000],baseName,hlabel,vlabel)
     doc.append(myImg)
##   Write ASCII data to a file
     asciiFileName = baseName + '.dat'
     asciiFile = open (asciiFileName, "w")
     n = len(x)
     for i in range(0,n):
##          asciiFile.write(repr(x[i]) + ' ' + repr(data[-1,i]) +'\n')
          asciiFile.write('%20.16e %20.16e\n' % (x[i], data[-1,i]))
     asciiFile.close()
     psFileName=baseName+'.ps'
     doc.append(BR())
     fileTable = Table()
     fileTable.body = [[Href(psFileName,'PostScript'), Href(asciiFileName,'ASCII data')]]
     fileTable.border=0
     fileTable.width='40%'
     fileTable.column1_align='center'
     fileTable.cell_align='center'
     doc.append(fileTable)
     return currNum
     

     
##      numVars=infile.CountVars()
##      for counter in range(0,numVars):
##           data=infile.ReadVar(counter)
##           if type(data)==numarray.numarraycore.NumArray:

##                varName=infile.GetVarName(counter)
##                doc.append(Name(sectionName+varName+repr(currNum)))
##                doc.append(Heading(2,varName))
##                myImg=ProduceCorrelationPicture(data[-1],varName+repr(currNum),'r',varName)

def ProduceTracePicture(data,fileBase,hlabel,vlabel):
    clf()
    x=fromfunction(lambda i:i,(len(data),))
    plot(x,data)
    h1=xlabel(hlabel)
    set(h1,"FontSize",20)
    v1=ylabel(vlabel)
    set(v1,"FontSize",20)
    labels = get(gca(), 'xticklabels')
    set(labels, 'fontsize', 16)
    labels = get(gca(), 'yticklabels')
    set(labels, 'fontsize', 16)
    savefig(fileBase+".png",dpi=60)
    savefig(fileBase+".ps")
    myImg=Image(fileBase+".png")
    return myImg



def MeanErrorString (mean, error):
     meanDigits = math.floor(math.log(abs(mean))/math.log(10))
     rightDigits = -math.floor(math.log(error)/math.log(10))+1
     formatstr = '%1.' + '%d' % rightDigits + 'f'
     meanstr  = formatstr % mean
     errorstr = formatstr % error
     return (meanstr, errorstr)

def ProcessScalarSection(infile,doc,currNum):
     sectionName=infile.GetName()
     doc.append(Heading(1,sectionName))
     toAddList=[]
    #put description
     myTable=Table()
     myTable.body=[['','Mean','Error','Variance', 'Kappa']]
     myTable.width='50%'
     numVars=infile.CountVars()
#     print "Num vars is ",numVars
     for counter in range(0,numVars):
          data=infile.ReadVar(counter)
          if type(data)==numarray.numarraycore.NumArray:
               currNum=currNum+1
               varName=infile.GetVarName(counter)
               baseName = varName+repr(currNum)
               toAddList.append(Name(sectionName+varName+repr(currNum)))
               toAddList.append(Heading(2,varName))
               myImg=ProduceTracePicture(data, baseName,'Blocks',varName)
               toAddList.append(myImg)

##   Write ASCII data to a file
               asciiFileName = baseName + '.dat'
               asciiFile = open (asciiFileName, "w")
               n = len(data)
               for i in range(0,n):
                    asciiFile.write('%20.16e\n' % data[i])
               asciiFile.close()
               psFileName=baseName+'.ps'
               doc.append(BR())
               fileTable = Table()
               fileTable.body = [[Href(psFileName,'PostScript'), Href(asciiFileName,'ASCII data')]]
               fileTable.border=0
               fileTable.width='40%'
               fileTable.column1_align='center'
               fileTable.cell_align='center'
               toAddList.append(fileTable)
               
               (mean,var,error,kappa)=stats.Stats(data)
               (meanstr, errorstr) = MeanErrorString (mean, error)
               myTable.body.append([Href("#"+sectionName+varName+repr(currNum),varName),\
                                    meanstr,errorstr, '%1.2e' % var ,'%1.2f' % kappa])
     doc.append(myTable)
     for counter in range(0,len(toAddList)):
          doc.append(toAddList[counter])
     return currNum


def ProcessSystemInfo(infile):
     tau=infile.ReadVar("tau")
     box=infile.ReadVar("Box")
     numTimeSlices=infile.ReadVar("NumTimeSlices")
     print numTimeSlices,tau
     beta=tau*numTimeSlices
     temp=1.0/beta
     systemTable=Table("System")
     systemTable.body=[["tau",repr(tau)]]
     systemTable.body.append(["# of Slices",repr(numTimeSlices)])
     systemTable.body.append(["beta", '%1.2f' % beta])
     systemTable.body.append(["temperature", '%1.4e' % temp])
     systemTable.body.append(["Box","[ "+'%1.2f' % box[0]+",  "+ '%1.2f' % box[1]+",  "+'%1.2f' % box[2]+" ]"])
     speciesTable=Table("Species")
     speciesTable.body=[]
     speciesTable.body.append(["Name","NumParticles","lambda","Type"])
     numSections=infile.CountSections2("Species")
     for spec in range(0,numSections):
          infile.OpenSection2("Species",spec)
          name=infile.ReadVar("Name")
          numPtcl=infile.ReadVar("NumParticles")
          lambdam = infile.ReadVar("lambda")
          type=infile.ReadVar("ParticleType")
          speciesTable.body.append([name,numPtcl,lambdam,type])
          infile.CloseSection()
     totalTable=Table()
     totalTable.body=[]
     totalTable.border=0
     totalTable.body.append([systemTable])
     totalTable.body.append([speciesTable])
     return totalTable
     
                       
                       


def ProcessRunInfo(infile):
     myTable=Table("Run Information")
     myTable.body=[]
     myTable.width='40%'
     numVars=infile.CountVars()
     for counter in range(0,numVars):
          data=infile.ReadVar(counter)
          varName=infile.GetVarName(counter)
          myTable.body.append([varName,data])
     return myTable



def ProcessTopTable(doc,infile):
     largeTable=Table()
     largeTable.border=0
     infile.OpenSection("RunInfo")
     runTable=ProcessRunInfo(infile)
     infile.CloseSection()
     infile.OpenSection("System")
     speciesTable=ProcessSystemInfo(infile)
     infile.CloseSection()
     largeTable.body.append([runTable,speciesTable])
     doc.append(largeTable)
     doc.append(HR())
     return doc

infile=IOSectionClass()
infile.OpenFile(sys.argv[1])
fileString=sys.argv[1]
dotLoc=string.rfind(fileString,'.')
dirName=fileString[0:dotLoc]
if not(os.access(dirName,os.F_OK)):
     os.mkdir(dirName)
os.chdir(dirName)
#infile.OpenSection("PathDump")
#pathData=GetPaths(infile)
#infile.CloseSection()
#(visualPath,visualBall)=InitVisualPaths(pathData)
#PlotPaths(pathData,visualPath,visualBall,0)

doc=SeriesDocument()
ProcessTopTable(doc,infile)


currNum=0
numSections=infile.CountSections()
print "The number of sections is ",numSections
for counter in range(0,numSections):
     infile.OpenSection(counter)
     print infile.GetName()
     myType=infile.ReadVar("Type")
     if myType=="Scalar":
          currNum=ProcessScalarSection(infile,doc,currNum)
          doc.append(HR())
     elif myType=="CorrelationFunction":
          currNum=ProcessCorrelationSection(infile,doc,currNum)
          doc.append(HR())
     infile.CloseSection()
doc.logo=""
doc.author="Ken and Bryan"
doc.email="bkclark@uiuc.edu"
doc.banner=("../analysis/pimcLogo.png")
doc.place_nav_buttons=0
doc.header()

doc.write("index.html")
