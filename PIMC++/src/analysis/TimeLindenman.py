#import numarray
import numpy
import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
from GraphDraw import *
import stats
    
def Process(infiles,summaryDoc,detailedDoc,StartCut):
    if StartCut==None:
        StartCut=0
    #acquire data about the correlation section
    baseName = "TimeLindenment" 
    hlabel="T"    #infiles.ReadVar("xlabel")[0]
    vlabel="Gamma_t" #infiles.ReadVar("ylabel")[0]
    data=infiles.ReadVar("TimeDisp")
    numSteps=infiles.ReadVar("NumStep")
    currNum = 0
    if data==None or data==[None]:
         return currNum
    s = shape(data[0])
    numProcs = len(data)
    (_,numBins) = shape(data[0])
    mShape = (numProcs, numBins)
    meanArray=zeros(mShape)+0.0
    errorArray=zeros(mShape)+0.0
    proc = 0
    for (d,s) in zip(data,numSteps): #loop over procs?
        for bin in range(0,numBins):
            totalSteps=sum(s[StartCut:-1,bin])
            dataSize=sum(s[StartCut:-1,bin]>0.5) # total amount of data
            print "Datasize is ",dataSize,totalSteps
            newData=zeros(dataSize)+0.0
            newDataCount=0
            for i in range(0,len(d[StartCut:-1,bin])):
                if s[StartCut+i,bin]>0.5: #i.e there is data here
                    newData[newDataCount]=d[StartCut+i,bin]/totalSteps
                    newDataCount=newDataCount+1
            (mean,var,error,kappa)=stats.Stats(newData)
            mean=mean*dataSize
            var=var*dataSize
            error=error*dataSize
            meanArray[proc,bin]=mean
            errorArray[proc,bin]=error
        proc = proc + 1
    mean  = zeros(numBins) + 0.0
    error = zeros(numBins) + 0.0
    for bin in range(0,numBins):
        (mean[bin],error[bin])=stats.UnweightedAvg(meanArray[:,bin],errorArray[:,bin])
        
#    x=infiles.ReadVar("x")[0]
    x=range(0,numBins)
    fillx = x
    revx= x[::-1]
    clf()
    fillx=concatenate((x,revx)) #reverse(x))
    y1=mean[1:-1]+error[1:-1]  #don't output the first point
    y2=mean[1:-1]-error[1:-1]  # don't output the first point
    filly=concatenate((y1,y2[::-1]))
    fill (fillx, filly,hold=True)
#    hold on
    plot (x, mean,'r',hold=True)
    xlabel ('r')
    ylabel ('g(r)')
    imageName = baseName + ".png"
    epsName = baseName +".eps"
    savefig(imageName,dpi=60)
    savefig(epsName,dpi=60)
    if (x==None):
         return currNum
    description=infiles.ReadVar("Description")[0]

    currNum=currNum+1

##Produce Image
    #myImg=ProduceCorrelationPicture(x, meanArray,baseName,hlabel,vlabel)
##Produce Ascii file
    asciiName = baseName + '.dat'
    asciiFile = open (asciiName, "w")
    n = len(x)
    for i in range(1,n): #don't write the first point
          asciiFile.write('%20.16e %20.16e %20.16e\n' \
                          % (x[i], mean[i], error[i]))
    asciiFile.close()
    
    #WriteAsciiFile(asciiFileName,x,y)

##create table with file names in it
    fileTable=BuildTable(1,2)
    fileTable.body = [[Href(epsName,'PostScript'), Href(asciiName,'ASCII data')]]
##Write things to document
    img = Image(imageName);
    summaryDoc.append(Heading(2,"Hexatic" ))
    summaryDoc.append(Heading(4,description))
    summaryDoc.append(img)
    summaryDoc.append(BR())
    summaryDoc.append(fileTable)
    print "I'm out of  hexatic"
    return currNum
