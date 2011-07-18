from Tables import *
from HTMLgen import *
from HTMLPlots import *
import numpy
import stats


def ProcessSuperfluidFraction(infiles,summaryDoc,detailedDoc,StartCut,temp):
    data=infiles.ReadVar("SuperfluidFraction")
    currNum = 0
    if (data==[None]):
         return currNum
    s = shape(data[0])
    numProcs = len(data)
    (_,numBins) = shape(data[0])
    dim = numBins
    # print s, numProcs, numBins
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
    meanSF=average(mean)
    errorSF=sqrt(sum([x*x for x in error]))

#    print "Superfluid",meanSF,errorSF
    SFTable=BuildTable(3,dim+2)
    if (dim==1):
      SFTable.body[0]=[" ","x","total"]
    elif (dim==2):
      SFTable.body[0]=[" ","x", "y","total"]
    elif (dim==3):
      SFTable.body[0]=[" ","x", "y","z","total"]
    else:
      print "Unrecognized dimension."
    SFTable.body[1][0]="Mean"
    SFTable.body[2][0]="Error"
    for i in range(1,dim+1):
        SFTable.body[1][i]= '%1.3f' % mean[i-1]
        SFTable.body[2][i]= '%1.3f' % error[i-1]
    SFTable.body[1][dim+1] = '%1.3f' % meanSF
    SFTable.body[2][dim+1] = '%1.3f' % errorSF

    clf()
    baseName="SuperfluidFraction"
    vlabel="$\rho_s$"
##    matplotlib.rc('text', usetex=True)

##    v1=ylabel(r"$\rho_s$")
    v1=ylabel("rho_s")
    setp(v1,"FontSize",20)
    labels = get(gca(), 'yticklabels')
    setp(labels, 'fontsize', 16)

    if (dim==1):
      barLoc = [0.1]
      errorLoc = [0.5]
      labels = ['x']
    elif (dim==2):
      barLoc = [0.1, 1.1]
      errLoc = [0.5, 1.5]
      labels = ['x','y']
    elif (dim==3):
      barLoc = [0.1, 1.1, 2.1]
      errLoc = [0.5, 1.5, 2.5]
      labels = ['x','y','z']
    else:
      print "Unrecognized Dimension."

    bar (barLoc, mean, yerr=error, ecolor='r')
    hold (True);
    errorbar (errLoc, mean, error, fmt='r.', capsize=15.0, ms=00.0);
    
    ticks = errLoc;
    xticks(ticks, labels)
    xticklabels = getp(gca(), 'xticklabels')
    setp (xticklabels, fontsize=18)
    hold (False);
    
    # errorbar([1,2,3],mean,error,fmt='b.')
    #axis([0, 3,0,1.5])

    DataDump = open ('SuperfluidFraction.dat', 'w')
    
    for row in range(0, numProcs):
      totalMean = 0
      totalError = 0
      for col in range(0, len(barLoc)):
        if (numProcs > 1):
          totalMean += mean[row][col]
          totalError += error[row][col]
        else:
          totalMean += mean[col]
          totalError += error[col]
      totalMean = totalMean/len(barLoc)
      totalError = totalError/len(barLoc)
      DataDump.write ('%1.6f ' % temp)
      DataDump.write ('%1.6f ' % totalMean)
      DataDump.write ('%1.6f ' % totalError)
      DataDump.write ('\n')
    DataDump.close()

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
    summaryDoc.append(Heading(2,"Superfluidfraction"))
    summaryDoc.append(Heading(dim+2,description))
    summaryDoc.append(SFTable)
    summaryDoc.append(img)
    summaryDoc.append(BR())
    summaryDoc.append(fileTable)
    return infiles,summaryDoc,detailedDoc





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
