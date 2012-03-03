#import numarray
import numpy
from IO import *
from Tables import *
from HTMLgen import *
from HTMLPlots import *
from GraphDraw import *
import stats
import pickle
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
        

def Processnofr(infiles,summaryDoc,detailedDoc,StartCut,box):
    print "Processing nofr"
    #acquire data about the correlation section
#    species1=infiles.ReadVar("Species1")[0]
#    species2=infiles.ReadVar("Species2")[0]
    hlabel="r"    #infiles.ReadVar("xlabel")[0]
    vlabel="n(r)" #infiles.ReadVar("ylabel")[0]
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
    print len(x), len(totalMean), len(totalError)
    first=0
    last=len(totalMean)
    for i in range(0,last):
      if totalMean[i]!=0:
        first=i
        break
    for i in range(0,last):
      if totalMean[i]<0.00001:
        totalMean[i]=0.00001
        print "Warning: Zero Value, setting to ", totalMean[i]
    x=x[first:last]
    totalMean=totalMean[first:last]
    totalError=totalError[first:last]
    print "First non-zero element: ", first, "(", x[0], ",", totalMean[0], ")"
    y=totalMean.copy()/totalMean[0]
    f = open('../nofrlineDump', 'a')
    pickle.dump(x.tolist(),f)
    pickle.dump(y.tolist(),f)
    pickle.dump(totalError.tolist(),f)
    f.close()
    myImgTable=ProduceErrorCorrelationPicture(x, y, totalError, baseName,hlabel,vlabel,box[0]/2)
    myImgTable2=ProduceCorrelationPictureLogLog(x, y, baseName+"loglog",hlabel,vlabel,box[0]/2)
    importanceFunction=GetImportanceFunction()
    y=totalMean.copy()
    y1=totalMean+totalError
    y2=totalMean-totalError
    for counter in range(0,len(y)):
        if y[counter]!=0:
            y[counter]=math.log(y[counter])/math.log(10)
#            y[counter]=math.log(y[counter]*math.exp(importanceFunction(x[counter]))/(x[counter]*x[counter]))/math.log(10)
        else:
            y[counter]=0
            print "Warning: Zero Value"
    for counter in range(0,len(y1)):
        if y1[counter]!=0:
            y1[counter]=math.log(y1[counter])/math.log(10)
#            y1[counter]=math.log(y1[counter]*math.exp(importanceFunction(x[counter]))/(x[counter]*x[counter]))/math.log(10)
        else:
            y1[counter]=0

    for counter in range(0,len(y2)):
        if y2[counter]!=0:
            y2[counter]=math.log(y2[counter])/math.log(10)
#            y2[counter]=math.log(y2[counter]*math.exp(importanceFunction(x[counter]))/(x[counter]*x[counter]))/math.log(10)
        else:
            y2[counter]=0
    yError=(y1-y2)/2
    baseName="nofrCorrected"
    correctedImgTable=ProduceErrorCorrelationPicture(x, y, yError, baseName,hlabel,vlabel,box[0]/2)
    #correctedImgTable2=ProduceCorrelationPictureLogLog(x, y, yError, baseName+"loglog",hlabel,vlabel,box[0]/2)
    sectionName=infiles.GetName()
    description=infiles.ReadVar("Description")[0]
    summaryDoc.append(Heading(1,sectionName))
    summaryDoc.append(Heading(4,description))
    summaryDoc.append(myImgTable)
    summaryDoc.append(myImgTable2)
    summaryDoc.append(correctedImgTable)
    #summaryDoc.append(correctedImgTable2)
    return 0
