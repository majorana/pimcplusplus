import numarray
import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
import numarray.fft
import stats
import math

def PrintArray(x,beta):
     for counter in range(0,len(x)/2):
          i=x[counter]
          print counter*2*3.14159/beta,i.real #math.sqrt(i.real*i.real+i.imag*i.imag)


def ProcessJosephson(infiles,summaryDoc,detailedDoc,StartCut,tau,numTimeSlices):
    hlabel="r"    #infiles.ReadVar("xlabel")[0]
    vlabel="g(r)" #infiles.ReadVar("ylabel")[0]
    data=infiles.ReadVar("PhiK")
#    data = map (lambda y: y[StartCut:],data)
    currNum=0
    beta=tau*numTimeSlices
    print "My beta is ",beta
    if (data==[None]):
         return currNum
    description=infiles.ReadVar("Description")[0]
    data=sum(data[0])/len(data[0])
    fftData=fft(data)
    for i in range(0,len(fftData)):
         fftData[i]=(0.0+i)*tau*fftData[i]/beta
    PrintArray(fftData,beta)
#    toPrint=fft(data)
#    for counter in range(0,len(data)):
#         print toPrint[counter].real
#    print "FFT DONE!",beta
#    for counter in range(0,len(data)):
#         print data[counter]
    


##     numProcs=len(data)
##     numDataPoints=len(data[0][0])
##     meanArrays=zeros((numProcs,numDataPoints))+0.0
##     errorArrays=zeros((numProcs,numDataPoints))+0.0
##     for proc in range(0,numProcs):
##        for pos in range(0,numDataPoints):
##            (meanArrays[proc,pos],var,errorArrays[proc,pos],kappa)=stats.Stats(data[proc][:,pos])
##     totalMean=zeros(numDataPoints)+0.0
##     totalError=zeros(numDataPoints)+0.0
##     for pos in range(0,numDataPoints):
##         (totalMean[pos],totalError[pos])=stats.WeightedAvg(meanArrays[:,pos],errorArrays[:,pos])
##     baseName="nofr"
    
## ##Produce Image
##     print len(x),len(totalMean),len(totalError)
##     myImgTable=ProduceCorrelationPictureWithErrors(x, totalMean,totalError,baseName,hlabel,vlabel)
    

##     sectionName=infiles.GetName()
##     description=infiles.ReadVar("Description")[0]
##     summaryDoc.append(Heading(1,sectionName))
##     summaryDoc.append(Heading(4,description))
##     summaryDoc.append(myImgTable)
    return 0
