import numarray
from IO import *
from Tables import *
from HTMLgen import *
from HTMLPlots import *
import stats
import math

def numericMap(func,myArray):
    newArray=zeros(len(myArray))+0.0
    for counter in range(0,len(myArray)):
        newArray[counter]=func(myArray[counter])
    return newArray

def WriteAsciiFile (asciiFileName,x,y):
     asciiFile = open (asciiFileName, "w")
     n = len(x)
     for i in range(0,n):
##          asciiFile.write(repr(x[i]) + ' ' + repr(data[-1,i]) +'\n')
          asciiFile.write('%20.16e %20.16e\n' % (x[i], y[i]))
     asciiFile.close()
     return

def GetImportanceFunction():
    infile = IOSectionClass()
    infile.OpenFile("pimc.in");
    infile.OpenSection("Action")
    infile.OpenSection("OpenLoop")
    ImportanceSample=infile.ReadVar("ImportanceSample")
    Xis=infile.ReadVar("Xis")
    if ImportanceSample=="Exponential" and Xis=="Distance":
        print "Now processing exponential importance sample"
        a=infile.ReadVar("a")
        alpha=infile.ReadVar("alpha")
        s=infile.ReadVar("s")
        return lambda x: -math.log(math.exp(-alpha*(x-s)*(x-s)))
    else:
        print "Code not yet written for said importance function"
        

def Processnofr(infiles,summaryDoc,detailedDoc,StartCut):
    #acquire data about the correlation section
#    species1=infiles.ReadVar("Species1")[0]
#    species2=infiles.ReadVar("Species2")[0]
    hlabel="r"    #infiles.ReadVar("xlabel")[0]
    vlabel="g(r)" #infiles.ReadVar("ylabel")[0]
    data=infiles.ReadVar("y")
    data = map (lambda y: y[StartCut:],data)
    if (data==[None]):
         return currNum
    x=infiles.ReadVar("x")[0]
    if (x==None):
         return currNum
    description=infiles.ReadVar("Description")[0]

    numProcs=len(data)
    numDataPoints=len(data[0][0])
    meanArrays=zeros((numProcs,numDataPoints))+0.0
    errorArrays=zeros((numProcs,numDataPoints))+0.0
    for proc in range(0,numProcs):
       for pos in range(0,numDataPoints):
           (meanArrays[proc,pos],var,errorArrays[proc,pos],kappa)=stats.Stats(data[proc][:,pos])
    totalMean=zeros(numDataPoints)+0.0
    totalError=zeros(numDataPoints)+0.0
    for pos in range(0,numDataPoints):
        (totalMean[pos],totalError[pos])=stats.WeightedAvg(meanArrays[:,pos],errorArrays[:,pos])
    baseName="nofr"
    
##Produce Image
    print len(x),len(totalMean),len(totalError)
    myImgTable=ProduceCorrelationPictureWithErrors(x, totalMean,totalError,baseName,hlabel,vlabel)
    importanceFunction=GetImportanceFunction()
    y=totalMean.copy()
    y1=totalMean+totalError
    y2=totalMean-totalError
    for counter in range(0,len(y)):
        if y[counter]!=0:
            y[counter]=math.log(y[counter]*math.exp(importanceFunction(x[counter])))/math.log(10)
#            y[counter]=math.log(y[counter]*math.exp(importanceFunction(x[counter]))/(x[counter]*x[counter]))/math.log(10)
        else:
            y[counter]=0
    for counter in range(0,len(y1)):
        if y1[counter]!=0:
            y1[counter]=math.log(y1[counter]*math.exp(importanceFunction(x[counter])))/math.log(10)
#            y1[counter]=math.log(y1[counter]*math.exp(importanceFunction(x[counter]))/(x[counter]*x[counter]))/math.log(10)
        else:
            y1[counter]=0

    for counter in range(0,len(y2)):
        if y2[counter]!=0:
            y2[counter]=math.log(y2[counter]*math.exp(importanceFunction(x[counter])))/math.log(10)
#            y2[counter]=math.log(y2[counter]*math.exp(importanceFunction(x[counter]))/(x[counter]*x[counter]))/math.log(10)
        else:
            y2[counter]=0
    yError=(y1-y2)/2
    baseName="nofrCorrected"
    correctedImgTable=ProduceCorrelationPictureWithErrors(x, y,yError,baseName,hlabel,vlabel)
    sectionName=infiles.GetName()
    description=infiles.ReadVar("Description")[0]
    summaryDoc.append(Heading(1,sectionName))
    summaryDoc.append(Heading(4,description))
    summaryDoc.append(myImgTable)
    summaryDoc.append(correctedImgTable)
    return 0
