import numarray
import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
import stats

def ProcessPairCorrelation(infiles,summaryDoc,detailedDoc,StartCut):
    #acquire data about the correlation section
    species1=infiles.ReadVar("Species1")[0]
    species2=infiles.ReadVar("Species2")[0]
    hlabel="r"    #infiles.ReadVar("xlabel")[0]
    vlabel="g(r)" #infiles.ReadVar("ylabel")[0]
    data=infiles.ReadVar("y")
    if (data==[None]):
         return currNum
#    data = map (lambda y: y[StartCut:-1],data)
    data=sum(data)/len(data)
    print shape(data)
    meanArray=zeros(len(data[0]))+0.0
    errorArray=zeros(len(data[0]))+0.0
    for col in range(0,len(data[0])):
        (mean,var,error,kappa)=stats.Stats(data[:][col])
        meanArray[col]=mean
        errorArray[col]=error
    x=infiles.ReadVar("x")[0]
    fillx = x
    revx= x[::-1]
    clf()
    fillx=concatenate((x,revx)) #reverse(x))
#    filly = meanArray+errorArray
    y1=meanArray+errorArray
    y2=meanArray-errorArray
    filly=concatenate((y1,y2[::-1]))
    print shape(fillx),shape(filly)
    fill (fillx, filly,hold=True)
#    hold on
    plot (x, meanArray,'r',hold=True)
    savefig("test.png",dpi=45)
    if (x==None):
         return currNum
    description=infiles.ReadVar("Description")[0]

    currNum=currNum+1
    baseName=sectionName+repr(currNum)

##Produce Image
    myImg=ProduceCorrelationPicture(x, y,baseName,hlabel,vlabel)
##Produce Ascii file
    asciiFileName = baseName + '.dat'
    WriteAsciiFile(asciiFileName,x,y)
    psFileName=baseName+'.ps'

##create table with file names in it
    fileTable=BuildTable()
    fileTable.body = [[Href(psFileName,'PostScript'), Href(asciiFileName,'ASCII data')]]
##Write things to document
    doc.append(Heading(1,sectionName))
    doc.append(Heading(4,description))
    doc.append(myImg)
    doc.append(BR())
    doc.append(fileTable)
    return currNum
