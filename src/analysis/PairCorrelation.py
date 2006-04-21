import numarray
import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
from GraphDraw import *
import stats

def ProcessPairCorrelation(infiles,summaryDoc,detailedDoc,StartCut):
    #acquire data about the correlation section
    species1=infiles.ReadVar("Species1")[0]
    species2=infiles.ReadVar("Species2")[0]
    baseName = "%s-%s_PC" % (species1, species2)
    hlabel="r"    #infiles.ReadVar("xlabel")[0]
    vlabel="g(r)" #infiles.ReadVar("ylabel")[0]
    data=infiles.ReadVar("y")
    currNum = 0
    if (data==[None]):
         return currNum
#    data = map (lambda y: y[StartCut:-1],data)
    data=sum(data)/len(data)
    meanArray=zeros(len(data[0]))+0.0
    errorArray=zeros(len(data[0]))+0.0
    for col in range(0,len(data[0])):
        (mean,var,error,kappa)=stats.Stats(data[:,col])
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
    fill (fillx, filly,hold=True)
#    hold on
    plot (x, meanArray,'r',hold=True)
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
    #WriteAsciiFile(asciiFileName,x,y)

##create table with file names in it
    fileTable=BuildTable(1,2)
    fileTable.body = [[Href(epsName,'PostScript'), Href(asciiName,'ASCII data')]]
##Write things to document
    img = Image(imageName);
    summaryDoc.append(Heading(2,"%s-%s pair correlation" \
                              % (species1, species2)))
    summaryDoc.append(Heading(4,description))
    summaryDoc.append(img)
    summaryDoc.append(BR())
    summaryDoc.append(fileTable)
    return currNum
