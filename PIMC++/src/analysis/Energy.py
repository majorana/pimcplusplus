import numpy

import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
import stats

def ProcessEnergy(infiles,summaryDoc,detailedDoc,StartCut):
    summaryDoc.append(Heading(2,"Energy"))
    description=infiles.ReadVar("Description")[0]
    summaryDoc.append(Heading(4,description))
    N=infiles.CountVars()
    varList = []
    numProcs=0
    for i in range(0,N):
        data = infiles.ReadVar(i)
        if (type(data[0])==numpy.ndarray):
            varList.append(i)
            numProcs=len(data)
    scalarVarTable=BuildTable(numProcs+1,len(varList)+1)
    scalarVarTable.body[0][0]="Proc"
    scalarTracePageHTMLList=[]
    summaryTable=BuildTable(len(varList)+1,5)
    summaryTable.body[0]=["Energy","Mean","Error","Variance","Kappa"]

    
    Etable = 1.0*numpy.zeros((numProcs, 2*len(varList)))
    col=0
    for var in varList:
        data = infiles.ReadVar(var)
        for proc in range(0,numProcs):
            (mean,var, error,kappa) = stats.Stats(data[proc][StartCut:-1])
            Etable[proc,2*col]   = mean
            Etable[proc,2*col+1] = error
        col = col + 1
    Efile = open ('Energies.dat', 'w')
    for row in range(0, numProcs):
        for col in range(0,2*len(varList)):
            Efile.write ('%1.6f ' % Etable[row][col])
        Efile.write ('\n')
    Efile.close()
            
            
    row = 0            
    for varCounter in varList:
        row = row + 1
        meanList=[]
        errorList=[]
        varList=[]
        kappaList=[]
        data = infiles.ReadVar(varCounter)
        varName = infiles.GetVarName(varCounter)
        baseName=varName+"Energy"
        if (len(data[0]) > 1 and numProcs<50):
            scalarTracePageHTMLList.append(BuildScalarTracePage(data,baseName,varName,StartCut))
        scalarVarTable.body[0][row]=varName
        for proc in range(0,numProcs):
            (mean,var,error,kappa)=stats.Stats(data[proc][StartCut:-1])
            meanList.append(mean)
            errorList.append(error)
            varList.append(var)
            kappaList.append(kappa)
            (meanstr, errorstr) = stats.MeanErrorString (mean, error)
            scalarVarTable.body[proc+1][0]=repr(proc)
            scalarVarTable.body[proc+1][row]=meanstr + "+/-" + errorstr
        (totalMean,totalError)=stats.UnweightedAvg(meanList,errorList)
        totalVar=sum(varList)/len(varList)
        totalKappa=sum(kappaList)/len(kappaList)
        (totalMeanStr,totalErrorStr)=stats.MeanErrorString(totalMean,totalError)
        # If the variance is essentially zero, the autocorrelation time is undefined
        if (totalVar < 1.0e-20):
            kappaString = "N/A"
        else:
            kappaString = '%1.2f' % totalKappa
        summaryTable.body[row]=[varName,totalMeanStr,totalErrorStr, '%1.2e' % totalVar, kappaString]
    summaryDoc.append(summaryTable)
    detailedDoc.append(scalarVarTable)
    for pageName in scalarTracePageHTMLList:
        myFrame=IFrame()
        myFrame.src=pageName
        myFrame.width="100%"
        myFrame.height="375"
        detailedDoc.append(myFrame)
    return infiles,summaryDoc,detailedDoc


