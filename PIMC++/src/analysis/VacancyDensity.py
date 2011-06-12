def ProcessVacancyDensity(infiles,summaryDoc,detailedDoc,StartCut):
    yData=infiles.ReadVar("y")
    print "done"
    yData=sum(yData[0])/len(yData)
    for z in range(0,len(yData[0][0])):
        for x in range(0,len(yData)):
            for y in range(0,len(yData[x])):
                print x,y,yData[x][y][z]

