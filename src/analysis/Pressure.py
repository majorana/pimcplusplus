from Tables import *
from HTMLgen import *
from HTMLPlots import *
import numarray
import stats

def ProcessPressure(infiles,summaryDoc,detailedDoc,StartCut):
#    TotalData=infiles.ReadVar("Total")
    meanList=[]
    errorList=[]
    varList=[]
    kappaList=[]
    variabList=[]
    N=infiles.CountVars()
    for i in range(0,N):
        data = infiles.ReadVar(i)
        if (type(data[0])==numarray.numarraycore.NumArray):
            variabList.append(i)
            numProcs=len(data)
    scalarVarTable=BuildTable(numProcs+1,len(variabList)+1)
    scalarVarTable.body[0][0]="Proc"
    scalarTracePageHTMLList=[]
    summaryTable=BuildTable(len(variabList)+1,5)
    summaryTable.body[0]=["Pressure","Mean","Error","Variance","Kappa"]
    row = 0
    for varCounter in variabList:
        row=row+1
        meanList=[]
        errorList=[]
        varList=[]
        kappaList=[]
        varName = infiles.GetVarName(varCounter)
        data = infiles.ReadVar(varCounter)
        numProcs=len(data)
        baseName=varName+"Pressure"
        scalarVarTable.body[0][row]=varName
        scalarTracePageHTMLList.append(BuildScalarTracePage(data,baseName,varName,StartCut))
        for proc in range(0,numProcs):
            (mean,var,error,kappa)=stats.Stats(data[proc][StartCut:-1])
            meanList.append(mean)
            errorList.append(error)
            varList.append(var)
            kappaList.append(kappa)
            (meanstr, errorstr) = stats.MeanErrorString (mean, error)
            scalarVarTable.body[proc+1][0]=repr(proc)
            scalarVarTable.body[proc+1][row]=meanstr + "+/-" + errorstr
        (totalMean,totalError)=stats.WeightedAvg(meanList,errorList)
        totalVar=sum(varList)/len(varList)
        totalKappa=sum(kappaList)/len(kappaList)
        (totalMeanStr,totalErrorStr)=stats.MeanErrorString(totalMean,totalError)
        summaryTable.body[row]=[varName,totalMeanStr,totalErrorStr,totalVar,totalKappa]


    summaryDoc.append(summaryTable)
    for pageName in scalarTracePageHTMLList:
        myFrame=IFrame()
        myFrame.src=pageName
        myFrame.width="100%"
        myFrame.height="375"
        detailedDoc.append(myFrame)
    detailedDoc.append(summaryTable)
    return infiles,summaryDoc,detailedDoc
