import numarray
import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
import stats




def ProcessBisectionBlock(infiles,summaryDoc,detailedDoc,StartCut):
    totalAcceptRatioVec=infiles.ReadVar("AcceptRatio")
    for proc in range(0,len(totalAcceptRatioVec)):
        totalAcceptRatioVec[proc]=totalAcceptRatioVec[proc][len(totalAcceptRatioVec[proc])-1]
    totalAcceptRatio=sum(totalAcceptRatioVec)/len(totalAcceptRatioVec)
    numStages=infiles.CountSections()
    infiles.OpenSection(0)
    acceptedPerms=infiles.ReadVar("Acceptance Ratio")
    print "Accepted perms: ",acceptedPerms
    print "Length",len(acceptedPerms[0])
    print "A:  ",acceptedPerms[proc][len(acceptedPerms[proc])-1]
    for proc in range(0,len(acceptedPerms)):
        acceptedPerms[proc]=acceptedPerms[proc][len(acceptedPerms[proc])-1]
    print "B: ",acceptedPerms
    abort()
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
