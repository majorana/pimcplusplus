from Tables import *
from HTMLgen import *
from HTMLPlots import *
import stats

def ProcessPressure(infiles,summaryDoc,detailedDoc,StartCut):
#    TotalData=infiles.ReadVar("Total")
    meanList=[]
    errorList=[]
    varList=[]
    kappaList=[]
    data = infiles.ReadVar("Total")
    numProcs=len(data)
    varName = "TotalPressure"
    baseName=varName
    scalarVarTable=BuildTable(numProcs+1,len(varList)+1)
    scalarVarTable.body[0][0]="Proc"
    scalarTracePageHTMLList=[]
    summaryTable=BuildTable(len(varList)+1,5)
    summaryTable.body[0]=["Pressure","Mean","Error","Variance","Kappa"]
    row = 0
    scalarVarTable.body[0][row]=varName
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
    (totalMean,totalError)=stats.WeightedAvg(meanList,errorList)
    totalVar=sum(varList)/len(varList)
    totalKappa=sum(kappaList)/len(kappaList)
    (totalMeanStr,totalErrorStr)=stats.MeanErrorString(totalMean,totalError)
    summaryTable.body.append([varName,totalMeanStr,totalErrorStr,totalVar,totalKappa])
    summaryDoc.append(summaryTable)
    for pageName in scalarTracePageHTMLList:
        myFrame=IFrame()
        myFrame.src=pageName
        myFrame.width="100%"
        myFrame.height="375"
        detailedDoc.append(myFrame)
    detailedDoc.append(summaryTable)
    return infiles,summaryDoc,detailedDoc
