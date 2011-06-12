from Tables import *
from HTMLgen import *
from HTMLPlots import *
##import numarray
import numpy
import stats


def ProcessWindingNumber(infiles,summaryDoc,detailedDoc,StartCut):
    #acquire data about the correlation section
#    species1=infiles.ReadVar("Species1")[0]
#    species2=infiles.ReadVar("Species2")[0]
#    baseName = "%s-%s_PC" % (species1, species2)
    baseName="WindingNumber"
    hlabel="r"    #infiles.ReadVar("xlabel")[0]
    vlabel="g(r)" #infiles.ReadVar("ylabel")[0]
    data=infiles.ReadVar("WindingNumber")
    currNum = 0
    if (data==[None]):
         return currNum
    s = shape(data[0])
    numProcs = len(data)
    (_,numBins) = shape(data[0])
    mShape = (numProcs, numBins)
    meanArray=zeros(mShape)+0.0
    errorArray=zeros(mShape)+0.0
    proc = 0
    for d in data:
        for bin in range(0,numBins):
            (mean,var,error,kappa)=stats.Stats(d[:,bin])
            meanArray[proc,bin]=mean
            errorArray[proc,bin]=error
        proc = proc + 1
    mean  = zeros(numBins) + 0.0
    error = zeros(numBins) + 0.0
    for bin in range(0,numBins):
        (mean[bin],error[bin])=stats.UnweightedAvg(meanArray[:,bin],errorArray[:,bin])
    print mean,error
#    x=infiles.ReadVar("x")[0]
    x=[1,2,3]
    fillx = x
    revx= x[::-1]
#    clf()
    fillx=concatenate((x,revx)) #reverse(x))
    filly = meanArray+errorArray
    y1=mean+error
    y2=mean-error
    filly=concatenate((y1,y2[::-1]))
    fill (fillx, filly,hold=True)
#    hold on
    plot (x, mean,'r',hold=True)
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
    summaryDoc.append(Heading(2,"WindingNumber"))
    summaryDoc.append(Heading(4,description))
    summaryDoc.append(img)
    summaryDoc.append(BR())
    summaryDoc.append(fileTable)
    return currNum




## def ProcessWindingNumber(infiles,summaryDoc,detailedDoc,StartCut):
##     meanList=[]
##     errorList=[]
##     varList=[]
##     kappaList=[]
##     variabList=[]
##     numProcs=0
##     N=infiles.CountVars()
##     for i in range(0,N):
##         data = infiles.ReadVar(i)
##         if (type(data[0])==numarray.numarraycore.NumArray):
##             variabList.append(i)
##             numProcs=len(data)
##     scalarVarTable=BuildTable(numProcs+1,len(variabList)+1)
##     scalarVarTable.body[0][0]="Proc"
##     scalarTracePageHTMLList=[]
##     summaryTable=BuildTable(len(variabList)+1,5)
##     summaryTable.body[0]=["Winding Number","Mean","Error","Variance","Kappa"]
##     row = 0
##     for varCounter in variabList:
##         row=row+1
##         meanList=[]
##         errorList=[]
##         varList=[]
##         kappaList=[]
##         varName = infiles.GetVarName(varCounter)
##         data = infiles.ReadVar(varCounter)
##         numProcs=len(data)
##         baseName=varName+"WindingNumber"
##         scalarVarTable.body[0][row]=varName
##         scalarTracePageHTMLList.append(BuildScalarTracePage(data,baseName,varName,StartCut))
##         for proc in range(0,numProcs):
##             (mean,var,error,kappa)=stats.Stats(data[proc][StartCut:-1])
##             meanList.append(mean)
##             errorList.append(error)
##             varList.append(var)
##             kappaList.append(kappa)
##             (meanstr, errorstr) = stats.MeanErrorString (mean, error)
##             scalarVarTable.body[proc+1][0]=repr(proc)
##             scalarVarTable.body[proc+1][row]=meanstr + "+/-" + errorstr
##         (totalMean,totalError)=stats.WeightedAvg(meanList,errorList)
##         totalVar=sum(varList)/len(varList)
##         totalKappa=sum(kappaList)/len(kappaList)
##         (totalMeanStr,totalErrorStr)=stats.MeanErrorString(totalMean,totalError)
##         summaryTable.body[row]=[varName,totalMeanStr,totalErrorStr,"%1.2e" % totalVar,"%1.2f" % totalKappa]

##     summaryDoc.append(summaryTable)
##     for pageName in scalarTracePageHTMLList:
##         myFrame=IFrame()
##         myFrame.src=pageName
##         myFrame.width="100%"
##         myFrame.height="375"
##         detailedDoc.append(myFrame)
##     detailedDoc.append(summaryTable)
##     return infiles,summaryDoc,detailedDoc
