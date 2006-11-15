import numarray
import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
from GraphDraw import *
import stats

def WriteAsciiFile (asciiFileName,x,y):
     asciiFile = open (asciiFileName, "w")
     n = len(x)
     for i in range(0,n):
##          asciiFile.write(repr(x[i]) + ' ' + repr(data[-1,i]) +'\n')
          asciiFile.write('%20.16e %20.16e\n' % (x[i], y[i]))
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
             meanArray[proc,bin]=mean
             errorArray[proc,bin]=error
         proc = proc + 1
     mean  = zeros(numBins) + 0.0
     error = zeros(numBins) + 0.0
     for bin in range(0,numBins):
         (mean[bin],error[bin])=stats.WeightedAvg(meanArray[:,bin],errorArray[:,bin])
     y=mean
     x=infiles.ReadVar("x")[0]
     if (x==None):
          return currNum
     toSort=[]
     for counter in range(0,len(x)):
          toSort.append((x[counter],y[counter]))
     toSort.sort()
     for counter in range(0,len(x)):
          x[counter]=toSort[counter][0]
          y[counter]=toSort[counter][1]

     xNew=[]
     yNew=[]
     counter=1
#     xNew.append(x[0])
#     yNew.append(y[0])
     while (counter<len(x)):
          totalY=y[counter]
          numY=1
          counter=counter+1
          while (counter<len(x) and x[counter]-x[counter-1]<1e-10):
               totalY=totalY+y[counter]
               numY=numY+1
               counter=counter+1
          xNew.append(x[counter-1])
          yNew.append(totalY/(numY+0.0))
     description=infiles.ReadVar("Description")[0]
     
     currNum=currNum+1
     baseName=sectionName+repr(currNum)

##Produce Image
     myImg=ProduceCorrelationPicture(x, y,baseName,hlabel,vlabel)
     myImg_avg=ProduceCorrelationPicture(xNew,yNew,baseName+"_avg",hlabel,vlabel)
##Produce Ascii file
     asciiFileName = baseName + '.dat'
     asciiFileName_avg = baseName + '_avg.dat'
     WriteAsciiFile(asciiFileName,x,y)
     WriteAsciiFile(asciiFileName_avg,xNew,yNew)
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
