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


def GetPaths(infiles):
     paths=infiles.ReadVar("Path")
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
#     infiles.OpenSection("grid")
#     x=infiles.ReadVar("Points")
#     infiles.CloseSection()
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
     currAxis=axis()
     currAxis[0]=cutoff
     print "Curr axis is ",currAxis
     axis(currAxis)
     savefig(fileBase+".png",dpi=60)
     savefig(fileBase+".ps")
     myImg=Image(fileBase+".png")
     return myImg


def ProcessCorrelationSection(infiles,doc,currNum):
     sectionName=infiles.GetName()
     doc.append(Heading(1,sectionName))
     hlabel=infiles.ReadVar("xlabel")[0]
     vlabel=infiles.ReadVar("ylabel")[0]
     data=infiles.ReadVar("y")[0]
     x=infiles.ReadVar("x")[0]
     description=infiles.ReadVar("Description")[0]
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
     

     
##      numVars=infiles.CountVars()
##      for counter in range(0,numVars):
##           data=infiles.ReadVar(counter)
##           if type(data)==numarray.numarraycore.NumArray:

##                varName=infiles.GetVarName(counter)
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
     if (mean!=0.0):
          meanDigits = math.floor(math.log(abs(mean))/math.log(10))
     else:
          meanDigits=2
     if (error!=0.0):
          rightDigits = -math.floor(math.log(error)/math.log(10))+1
     else:
          rightDigits=2
     formatstr = '%1.' + '%d' % rightDigits + 'f'
     meanstr  = formatstr % mean
     errorstr = formatstr % error
     return (meanstr, errorstr)

def Avg (x):
     if x[0] == None:
          return None
     else:
          return sum(x)/len(x)

def WeightedAvg (means, errors):
     if (errors[0] != 0.0):
          weights = map (lambda x: 1.0/(x*x), errors)
          norm = 1.0/sum(weights)
          weights = map(lambda x: x*norm, weights)
          avg = 0.0
          error2 = 0.0
          for i in range (0,len(means)):
               avg = avg + means[i]*weights[i]
               error2 = error2 + weights[i]*errors[i]*errors[i]
          return (avg, math.sqrt(error2))
     else:
          return (Avg(means), 0.0)


def ProcessScalarSection(infiles,doc,currNum):
     sectionName=infiles.GetName()
     doc.append(Heading(1,sectionName))
     toAddList=[]
    #put description
     myTable=Table()
     myTable.body=[['','Mean','Error','Variance', 'Kappa']]
     myTable.width='50%'
     numVars=infiles.CountVars()
#     print "Num vars is ",numVars
     for counter in range(0,numVars):
          data = infiles.ReadVar(counter)
          if type(data[0])==numarray.numarraycore.NumArray:
               currNum=currNum+1
               varName=infiles.GetVarName(counter)
               baseName = varName+repr(currNum)
               toAddList.append(Name(sectionName+varName+repr(currNum)))
               toAddList.append(Heading(2,varName))
               myImg=ProduceTracePicture(data[0], baseName,'Blocks',varName)
               toAddList.append(myImg)

##   Write ASCII data to a file
               asciiFileName = baseName + '.dat'
               asciiFile = open (asciiFileName, "w")
               n = len(data)
               for i in range(0,n):
                    asciiFile.write('%20.16e\n' % data[0][i])
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
               meanlist = []
               varlist = []
               errorlist = []
               kappalist = []
               for d in data:
                    (mean,var,error,kappa)=stats.Stats(d)
                    meanlist.append(mean)
                    errorlist.append(error)
                    varlist.append(var)
                    kappalist.append(kappa)
               (mean,error) = WeightedAvg(meanlist, errorlist)
               print repr(mean) + "+/-" + repr(error)
                    
               (meanstr, errorstr) = MeanErrorString (mean, error)
               myTable.body.append([Href("#"+sectionName+varName+repr(currNum),varName),\
                                    meanstr,errorstr, '%1.2e' % var ,'%1.2f' % kappa])
     doc.append(myTable)
     for counter in range(0,len(toAddList)):
          doc.append(toAddList[counter])
     return currNum


def ProcessSystemInfo(infiles):
     tau=infiles.ReadVar("tau")[0]
     box=infiles.ReadVar("Box")[0]
     numTimeSlices=infiles.ReadVar("NumTimeSlices")[0]
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
     numSections=infiles.CountSections2("Species")
     for spec in range(0,numSections):
          infiles.OpenSection2("Species",spec)
          name=infiles.ReadVar("Name")[0]
          numPtcl=infiles.ReadVar("NumParticles")[0]
          lambdam = infiles.ReadVar("lambda")[0]
          type=infiles.ReadVar("ParticleType")[0]
          speciesTable.body.append([name,numPtcl,lambdam,type])
          infiles.CloseSection()
     totalTable=Table()
     totalTable.body=[]
     totalTable.border=0
     totalTable.body.append([systemTable])
     totalTable.body.append([speciesTable])
     return totalTable
     
                       
                       


def ProcessRunInfo(infiles):
     myTable=Table("Run Information")
     myTable.body=[]
     myTable.width='40%'
     numVars=infiles.CountVars()
     for counter in range(0,numVars):
          data=infiles.ReadVar(counter)[0]
          varName=infiles.GetVarName(counter)
          myTable.body.append([varName,data])
     return myTable



def ProcessTopTable(doc,infiles):
     largeTable=Table()
     largeTable.border=0
     infiles.OpenSection("RunInfo")
     runTable=ProcessRunInfo(infiles)
     infiles.CloseSection()
     infiles.OpenSection("System")
     speciesTable=ProcessSystemInfo(infiles)
     infiles.CloseSection()
     largeTable.body.append([runTable,speciesTable])
     doc.append(largeTable)
     doc.append(HR())
     return doc

def ProcessMove(doc,infiles):
     doc.append(Heading(1,"Moves"))
     myTable=Table()
     myTable.body=[['Moves','Acceptance']]
     myTable.width='50%'
     infiles.OpenSection("Moves")
     numMoves=infiles.CountSections()
     for i in range (0, numMoves):
          infiles.OpenSection(i)
          name=infiles.GetName()
          ar = Avg(infiles.ReadVar("AcceptRatio"))
          if (ar!=None):
               totAccept=ar[0]
               numAccept=0
               for counter in range(1,len(ar)):
                    totAccept=totAccept+ar[counter]
                    numAccept=numAccept+1
               myTable.body.append([name+" ",totAccept/numAccept])
          else:
               myTable.body.append([name,"No acceptance available"])
          infiles.CloseSection()
     infiles.CloseSection() # "Moves"
     doc.append(myTable)



basename = sys.argv[1]
infiles = IOSectionClassList()
infiles.OpenFiles(basename);
      
print 'Found ' +repr(infiles.len()) + ' output files.'

dirName=basename 
cutoff=None
if (os.access(dirName+".pref",os.F_OK)):
     print dirName+".pref"
     prefFile.OpenFile(dirName+".pref")
     print prefFile.ReadVar("cutoff")
     cutoff=prefFile.ReadVar("cutoff")
     prefFile.CloseFile() 
if cutoff==None:
     cutoff=0
if not(os.access(dirName,os.F_OK)):
     os.mkdir(dirName)
os.chdir(dirName)

doc=SeriesDocument()
ProcessTopTable(doc,infiles)
ProcessMove(doc,infiles)

currNum=0
infiles.OpenSection("Observables")
numSections=infiles.CountSections()
#print "The number of sections is ",numSections
for counter in range(0,numSections):
     infiles.OpenSection(counter)
     print infiles.GetName()
     myType=infiles.ReadVar("Type")[0]
     print "myType = " + myType
     if myType=="Scalar":
          currNum=ProcessScalarSection(infiles,doc,currNum)
          doc.append(HR())
     elif myType=="CorrelationFunction":
          currNum=ProcessCorrelationSection(infiles,doc,currNum)
          doc.append(HR())
     infiles.CloseSection()
infiles.CloseSection() # "Observables"


doc.logo=""
doc.author="Ken and Bryan"
doc.email="bkclark@uiuc.edu"
doc.banner=("../analysis/pimcLogo.png")
doc.place_nav_buttons=0
doc.header()

doc.write("index.html")
