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
      if i != 0 and i != 10:
        varList.append(i)
    numProcs = infiles.len()
    #for i in range(0,N):
    #    print i
    #    data = infiles.ReadVar(i)
    #    print data[0], infiles.GetVarName(i)
    #    if (type(data[0])==numpy.ndarray):
    #      varList.append(i)
    #      numProcs=len(data)
    print varList, numProcs
    scalarVarTable=BuildTable(numProcs+1,len(varList)+1)
    scalarVarTable.body[0][0]="Proc"
    scalarTracePageHTMLList=[]
    summaryTable=BuildTable(len(varList)+1,5)
    summaryTable.body[0]=["Energy","Mean","Error","Variance","Kappa"]

#    for varCounter in varList:
#        data = infiles.ReadVar(varCounter)
#        for proc in range(0,numProcs):
#            nodeblows = len(data[proc])
#            data[proc] = data[proc][data[proc]<1e90]
#            nodeblows -= len(data[proc])
#            #for i in range(0,len(data[proc])):
#            #  if (data[proc][i].ndim == 0 and data[proc][i] > 1e90):
#            #    print data[proc][i]
#            #    print data[proc][i]
#            #    nodeblows += 1
#            if len(data[proc]) > StartCut:
#              (mean,var, error,kappa) = stats.Stats(data[proc][StartCut:-1])
#            else:
#              (mean,var, error,kappa) = stats.Stats(data[proc])
#            Etable[proc,2*col]   = mean
#            Etable[proc,2*col+1] = error
#            if nodeblows > 0:
#              print varCounter, proc, "Node Blow Ups:", nodeblows
#        col = col + 1


    Etable = 1.0*numpy.zeros((numProcs, 2*len(varList)))
    col=0
    row = 0
    AvgEfile = open ('AvgEnergies.dat','w')
    for varCounter in varList:
        row = row + 1
        meanList=[]
        errorList=[]
        variList=[]
        kappaList=[]
        #if (len(data[0]) > 1 and numProcs<50):
        #if (len(data[0]) > 1):
            #scalarTracePageHTMLList.append(BuildScalarTracePage(data,baseName,varName,StartCut))
        #scalarVarTable.body[0][row]=varName
        for proc in range(0,numProcs):
            data = infiles.IOlist[proc].ReadVar(varCounter)
            varName = infiles.IOlist[proc].GetVarName(varCounter)
            baseName=varName+"Energy"
            nodeblows = len(data)
            data = data[data<1e90]
            nodeblows -= len(data)
            #for i in range(0,len(data[proc])):
            #  if (data[proc][i].ndim == 0  and data[proc][i] > 1e90):

            #    numpy.delete(data[proc],i)
            #    nodeblows += 1
            if len(data) > StartCut:
              (mean,var,error,kappa)=stats.Stats(data[StartCut:-1])
            else:
              (mean,var,error,kappa)=stats.Stats(data)
            Etable[proc,2*col]   = mean
            Etable[proc,2*col+1] = error
            meanList.append(mean)
            errorList.append(error)
            variList.append(var)
            kappaList.append(kappa)
            (meanstr, errorstr) = stats.MeanErrorString (mean, error)
            #scalarVarTable.body[proc+1][0]=repr(proc)
            #scalarVarTable.body[proc+1][row]=meanstr + "+/-" + errorstr
            if nodeblows > 0:
              print varCounter, proc, "Node Blow Ups:", nodeblows
        (totalMean,totalError)=stats.UnweightedAvg(meanList,errorList) 
        print varName, totalMean, totalError
        AvgEfile.write (varName+' '+str(totalMean)+' '+str(totalError)+'\n')
        totalVar=sum(variList)/len(variList)
        totalKappa=sum(kappaList)/len(kappaList)
        (totalMeanStr,totalErrorStr)=stats.MeanErrorString(totalMean,totalError)
        # If the variance is essentially zero, the autocorrelation time is undefined
        if (totalVar < 1.0e-20):
            kappaString = "N/A"
        else:
            kappaString = '%1.2f' % totalKappa
        #summaryTable.body[row]=[varName,totalMeanStr,totalErrorStr, '%1.2e' % totalVar, kappaString]
        col = col + 1
    AvgEfile.close()
    Efile = open ('Energies.dat', 'w')
    for row in range(0, numProcs):
        for col in range(0,2*len(varList)):
            Efile.write ('%1.6f ' % Etable[row][col])
        Efile.write ('\n')
    Efile.close()
    summaryDoc.append(summaryTable)
    detailedDoc.append(scalarVarTable)
    for pageName in scalarTracePageHTMLList:
        myFrame=IFrame()
        myFrame.src=pageName
        myFrame.width="100%"
        myFrame.height="375"
        detailedDoc.append(myFrame)
    return infiles,summaryDoc,detailedDoc


