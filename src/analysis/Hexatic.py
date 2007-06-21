import numarray
import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
from GraphDraw import *
import stats
    
def ProcessHexatic(infiles,summaryDoc,detailedDoc,StartCut):
    print "I'm in side hexatic"
    #acquire data about the correlation section
    baseName = "Hexatic" 
    hlabel="r"    #infiles.ReadVar("xlabel")[0]
    vlabel="g_6(r)" #infiles.ReadVar("ylabel")[0]
    data=infiles.ReadVar("HexaticReal")
    currNum = 0
    if data==None or data==[None]:
         return currNum
    s = shape(data[0])
    numProcs = len(data)
    (_,numBins) = shape(data[0])
##    #hack
##    numBin=int(numBins/3)
##    #end hack
    mShape = (numProcs, numBins)
    meanArray=zeros(mShape)+0.0
    errorArray=zeros(mShape)+0.0
    proc = 0
    for d in data:
        for bin in range(0,numBins): 
            (mean,var,error,kappa)=stats.Stats(d[StartCut:-1,bin])
#            (mean,var,error,kappa)=stats.Stats(d[StartCut:-1,bin])
            meanArray[proc,bin]=mean
            errorArray[proc,bin]=error
        proc = proc + 1
    mean  = zeros(numBins) + 0.0
    error = zeros(numBins) + 0.0
    for bin in range(0,numBins):
        (mean[bin],error[bin])=stats.UnweightedAvg(meanArray[:,bin],errorArray[:,bin])
        
    x=infiles.ReadVar("x")[0]
#    x=range(0,numBins)
    fillx = x
    revx= x[::-1]
    clf()
    fillx=concatenate((x,revx)) #reverse(x))
    for bin in range(0,numBins):
        deltax=2*3.14159*x[bin]*(x[1]-x[0])
        if (bin!=numBins-1):
            deltax=3.14159*(x[bin+1]**2-x[bin]**2)
            #42.0=box lenght, 562=ptcl number
        boxSize=21.121225
        nump=142
        mean[bin]=mean[bin]/(deltax)*boxSize*boxSize/(nump*nump)
#        mean[bin]=mean[bin]/(deltax)*42.0*42.0/(562.0*562.0)
        
#    filly = meanArray+errorArray
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
