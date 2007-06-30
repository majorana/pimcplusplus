#import numarray
import numpy
import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
from GraphDraw import *
import stats

def WriteAsciiFile_1 (asciiFileName,x,y):
     asciiFile = open (asciiFileName, "w")
     n = len(x)
     for i in range(0,n):
##          asciiFile.write(repr(x[i]) + ' ' + repr(data[-1,i]) +'\n')
          asciiFile.write('%20.16e %20.16e\n' % (x[i], y[i]))
     asciiFile.close()
     return

def WriteAsciiFile (asciiFileName,x,y,z):
     asciiFile = open (asciiFileName, "w")
     n = len(x)
     for i in range(0,n):
##          asciiFile.write(repr(x[i]) + ' ' + repr(data[-1,i]) +'\n')
          asciiFile.write('%20.16e %20.16e %20.16e\n' % (x[i], y[i], z[i]))
     asciiFile.close()
     return




def ProcessStructureFactor(infiles,summaryDoc,detailedDoc,StartCut):
     #acquire data about the structure factor
     currNum=0
     sectionName=infiles.GetName()
     hlabel=infiles.ReadVar("xlabel")[0]
     vlabel=infiles.ReadVar("ylabel")[0]
     data=infiles.ReadVar("y")
     if (data==[None]):
          return currNum
     numProcs = len(data)
     (_,numBins) = shape(data[0])
     mShape = (numProcs, numBins)
     meanArray=zeros(mShape)+0.0
     errorArray=zeros(mShape)+0.0
     proc = 0
     for d in data:
         for bin in range(0,numBins):
             (mean,var,error,kappa)=stats.Stats(d[StartCut:-1,bin])
##             (mean,var,error,kappa)=stats.Stats(d[StartCut:StartCut+1,bin])
##             (mean,var,error,kappa)=stats.Stats(d[0:50,bin])
             meanArray[proc,bin]=mean
             errorArray[proc,bin]=error
         proc = proc + 1
     mean  = zeros(numBins) + 0.0
     error = zeros(numBins) + 0.0
     for bin in range(0,numBins):
#         (mean[bin],error[bin])=stats.WeightedAvg(meanArray[:,bin],errorArray[:,bin])
         (mean[bin],error[bin])=stats.UnweightedAvg(meanArray[:,bin],errorArray[:,bin])
     y=mean
     x=infiles.ReadVar("x")[0]
     kVecs=infiles.ReadVar("kVecs")[0]
     print "KVECS Coming"
     print kVecs
     print "KVECS DONE"
     if (x==None):
          return currNum
     toSort=[]
     for counter in range(0,len(x)):
          toSort.append((x[counter],y[counter],error[counter]))
     toSort.sort()
     for counter in range(0,len(x)):
          x[counter]=toSort[counter][0]
          y[counter]=toSort[counter][1]
          error[counter]=toSort[counter][2]

     xNew=[]
     yNew=[]
     yErrorNew=[]
     counter=1
#     xNew.append(x[0])
#     yNew.append(y[0])
     while (counter<len(x)):
          totalY=y[counter]
          totalYerr=error[counter]*error[counter]
          numY=1
          counter=counter+1
          while (counter<len(x) and x[counter]-x[counter-1]<1e-10):
               totalY=totalY+y[counter]
               totalYerr=totalYerr+error[counter]*error[counter]
               numY=numY+1
               counter=counter+1
          xNew.append(x[counter-1])
          yNew.append(totalY/(numY+0.0))
          yErrorNew.append(sqrt(totalYerr)/(numY+0.0))
          
     description=infiles.ReadVar("Description")[0]
     
     currNum=currNum+1
     baseName=sectionName+repr(currNum)+'.'+repr(StartCut)

##Produce Image
     myImg=ProduceCorrelationPicture(x, y,baseName,hlabel,vlabel)
     myImg_avg=ProduceCorrelationPicture(xNew,yNew,baseName+"_avg",hlabel,vlabel)
##Produce Ascii file
     asciiFileName = baseName + '.dat'
     asciiFileName_avg = baseName + '_avg.dat'
     WriteAsciiFile(asciiFileName,x,y,error)
     WriteAsciiFile(asciiFileName_avg,xNew,yNew,yErrorNew)
     xkVecs=[]
     ykVecs=[]
     for i in range(0,len(kVecs)):
          xkVecs.append(kVecs[i][0])
          ykVecs.append(kVecs[i][1])
     WriteAsciiFile("SFVecs",xkVecs,ykVecs,y)
     psFileName=baseName+'.ps'
     psFileName_avg=baseName+'_avg.ps'

##create table with file names in it
     fileTable=BuildTable(1,2)
     fileTable.body = [[Href(psFileName,'PostScript'), Href(asciiFileName,'ASCII data')]]
     fileTable_avg=BuildTable(1,2)
     fileTable_avg.body = [[Href(psFileName_avg,'PostScript'), Href(asciiFileName_avg,'ASCII data')]]
     
##Write things to summaryDocument
     summaryDoc.append(Heading(1,sectionName))
     summaryDoc.append(Heading(4,description))
     summaryDoc.append(myImg)
     summaryDoc.append(BR())
     summaryDoc.append(fileTable)
     summaryDoc.append(myImg_avg)
     summaryDoc.append(BR())
     summaryDoc.append(fileTable_avg)
     return currNum
