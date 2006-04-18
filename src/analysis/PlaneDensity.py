def ProcessPlaneDensity(infiles,summaryDoc,detailedDoc,StartCut):
    yData=infiles.ReadVar("y")
    print "done"
    yData=sum(yData[0])/len(yData)
    for counter in range(0,len(yData)):
        for counter2 in range(0,len(yData[counter])):
            print yData[counter][counter2],
            print " ",
        print
