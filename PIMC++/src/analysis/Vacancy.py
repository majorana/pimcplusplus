def ProcessVacancy(infiles,summaryDoc,detailedDoc,StartCut):
    grData=infiles.ReadVar("g(r)")
#    yData=infiles.ReadVar("y")
    procAverages=[]
#    yprocAverages=[]
    for counter in range(0,len(grData)):
        procAverages.append(sum(grData[counter][StartCut:])/len(grData[counter][StartCut:]))
#        yprocAverages.append(sum(yData[counter][StartCut:])/len(yData[counter][StartCut:]))
#    yfullAverage=sum(yprocAverages)/len(yprocAverages)
    fullAverage=sum(procAverages)/len(procAverages)
    print "Printing average"
    for counter in range(0,len(fullAverage)):
        print fullAverage[counter]    
#    for counter in range(0,len(yfullAverage)):
#        print counter, yfullAverage[counter]
#    print sum(yfullAverage)/len(yfullAverage)
