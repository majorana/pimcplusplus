import numpy
import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
import stats


def Process(infiles,summaryDoc,detailedDoc,StartCut):
    PermuteTable=[]
    totalAcceptRatioVec=infiles.ReadVar("AcceptRatio")
    for proc in range(0,len(totalAcceptRatioVec)):
        totalAcceptRatioVec[proc]=sum(totalAcceptRatioVec[proc])/len(totalAcceptRatioVec[proc])
    totalAcceptRatio=sum(totalAcceptRatioVec)/len(totalAcceptRatioVec)
    numStages=infiles.CountSections()
    AcceptTable=BuildTable(2,numStages+1)
    AcceptTable.body[0][0:1]=["Move","Total Acceptance Ratio"]
    AcceptTable.body[0][2:]=["Stage "+str(i) for i in range(0,numStages)]
    AcceptTable.body[1][0:1]=["Bisection Block",totalAcceptRatio]
    for stage in range(0,numStages):
        infiles.OpenSection(stage)
        acceptedPerms=infiles.ReadVar("AcceptRatio")
        if acceptedPerms[0]==None:
            print "Reading permutation data"
            acceptedPerms=infiles.ReadVar("Acceptance Ratio")
            permsTried=infiles.ReadVar("Perms Tried")
            acceptedPerms = [x*y for (x,y) in zip(acceptedPerms,permsTried)]
            acceptedPerms = [sum(x) for x in acceptedPerms] 
            permsTried = [sum(x) for x in permsTried]
            acceptedPerms=sum(acceptedPerms)
            permsTried=sum(permsTried)
            PermuteTable=BuildTable(4,5)
            PermuteTable.body[0][0]=""
            PermuteTable.body[1][0]="Accepted Permutations"
            PermuteTable.body[2][0]="Attempted Permutations"
            PermuteTable.body[3][0]="Percent Accepted"
            for dim in range(1,5):
                PermuteTable.body[0][dim]=str(dim)+" ptcl"
                PermuteTable.body[1][dim]=acceptedPerms[dim-1]
                PermuteTable.body[2][dim]=permsTried[dim-1]
                PermuteTable.body[3][dim]=acceptedPerms[dim-1]/permsTried[dim-1]
        else:
            acceptedPerms=[average(x) for x in acceptedPerms]
            acceptedPerms=sum(acceptedPerms)/len(acceptedPerms)
            AcceptTable.body[1][stage+2]=acceptedPerms
#        AcceptTable.body[0][stage+2]="Stage "+str(stage)
        infiles.CloseSection()
    summaryDoc.append(AcceptTable)
    if PermuteTable!=[]:
        summaryDoc.append(PermuteTable)
