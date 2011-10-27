##import numarray
import numpy
import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
import stats


def ProcessDisplaceMove(infiles,summaryDoc,detailedDoc,StartCut):
    totalAcceptRatioVec=infiles.ReadVar("AcceptRatio")
    for proc in range(0,len(totalAcceptRatioVec)):
        totalAcceptRatioVec[proc]=sum(totalAcceptRatioVec[proc])/len(totalAcceptRatioVec[proc])
    totalAcceptRatio=sum(totalAcceptRatioVec)/len(totalAcceptRatioVec)
    AcceptTable=BuildTable(1,2)
    AcceptTable.body[0][0]="Displace"
    AcceptTable.body[0][1]=totalAcceptRatio    
    summaryDoc.append(AcceptTable)

def ProcessBisectionBlock(infiles,summaryDoc,detailedDoc,StartCut):
    totalAcceptRatioVec=infiles.ReadVar("AcceptRatio")
    for proc in range(0,len(totalAcceptRatioVec)):
        totalAcceptRatioVec[proc]=sum(totalAcceptRatioVec[proc])/len(totalAcceptRatioVec[proc])
    totalAcceptRatio=sum(totalAcceptRatioVec)/len(totalAcceptRatioVec)
    numStages=infiles.CountSections()
    AcceptTable=BuildTable(2,numStages+2)
    AcceptTable.body[0][0]="Move"
    AcceptTable.body[0][1]="Total Acceptance Ratio"
    AcceptTable.body[1][0]="Bisection Block"
    AcceptTable.body[1][1]=totalAcceptRatio
    for stage in range(0,numStages):
        infiles.OpenSection(stage)
        acceptedPerms=infiles.ReadVar("AcceptRatio")
        for proc in range(0,len(acceptedPerms)):
                acceptedPerms[proc]=sum(acceptedPerms[proc])/len(acceptedPerms[proc])
        acceptedPerm=sum(acceptedPerms)/len(acceptedPerms)
        AcceptTable.body[0][stage+2]="Stage "+str(stage)
        AcceptTable.body[1][stage+2]=acceptedPerm
        print acceptedPerm
        infiles.CloseSection()
    summaryDoc.append(AcceptTable)
#    
#    
#    for stage in range(0,numStages):
#        infiles.OpenSection(stage)
#        
#        #        stageName = infiles.GetName()



##     for secNum in range(0, numSections):
##         infiles.OpenSection(secNum)


    
##     N=infiles.CountVars()
##     varList = []
##     numProcs=0
##     for i in range(0,N):
##         data = infiles.ReadVar(i)
##         if (type(data[0])==numarray.numarraycore.NumArray):
##             varList.append(i)
##             numProcs=len(data)
##     scalarVarTable=BuildTable(numProcs+1,len(varList)+1)
##     scalarVarTable.body[0][0]="Proc"
##     scalarTracePageHTMLList=[]
##     summaryTable=BuildTable(len(varList)+1,5)
##     summaryTable.body[0]=["Energy","Mean","Error","Variance","Kappa"]
