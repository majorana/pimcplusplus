import numarray
import IO
from Tables import *
from HTMLgen import *
from HTMLPlots import *
import stats

def ProcessTimeAnalysis(infiles,summaryDoc,detailedDoc,StartCut):
    MoveTime     = infiles.ReadVar("MoveTime")
    numMoves = MoveTime[0].shape[1]
    numProcs = len(MoveTime)
    print MoveTime[0].shape
    MT = 0.0*arange(0,numMoves)
    for t in MoveTime:
        numT = t.shape[0]
        print t[-1,:].shape
        MT = MT + t[-1,:]
    MT = MT / numProcs
    print MT

    ObserveTime  = infiles.ReadVar("ObservableTime")
    numObserves = ObserveTime[0].shape[1]
    print ObserveTime[0].shape
    OT = 0.0*arange(0,numObserves)
    for t in ObserveTime:
        numT = t.shape[0]
        print t[-1,:].shape
        OT = OT + t[-1,:]
    OT = OT / numProcs
    print OT

    TotalTime    = infiles.ReadVar("TotalTime")
    MoveNames    = infiles.ReadVar("MoveNames")[0]
    ObserveNames = infiles.ReadVar("ObserveNames")[0]
    if (MoveNames == None):
        MoveNames = [ "1" ]
        for i in range(2, numMoves+1):
            MoveNames.append(repr(i))
    print MoveNames
    if (ObserveNames == None):
        ObserveNames = [ "1" ]
        for i in range(2, numObserves+1):
            ObserveNames.append(repr(i))

    MoveTable = BuildTable(len(MT)+1,2)
    MoveTable.width = '20%'
    MoveTable.body[0] = ["Move Name", "% Time"]
    row = 0
    for i in range(0,numMoves):
        MoveTable.body[i+1] = [MoveNames[i], "%5.2f" % (MT[i]*100.0)]
    summaryDoc.append(MoveTable)

    ObserveTable = BuildTable(len(OT)+1,2)
    ObserveTable.width = '20%'
    ObserveTable.body[0] = ["Observe Name", "% Time"]
    row = 0
    for i in range(0,numObserves):
        ObserveTable.body[i+1] = [ObserveNames[i], "%5.2f" % (OT[i]*100.0)]
    summaryDoc.append(ObserveTable)

