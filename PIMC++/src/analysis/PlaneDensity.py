##from pylab import *

def ProcessPlaneDensity(infiles,summaryDoc,detailedDoc,StartCut,box):
    yData=infiles.ReadVar("y")
    print "done",len(yData)
#    yData=sum(yData[0][StartCut:])/len(yData[0][StartCut:])
    yData=sum(yData[0][90:98])/len(yData[0][90:98])
    print shape(yData)
    newyData=zeros((len(yData)*2,len(yData[0])*2))+0.0
    for i in range(0,len(yData)):
        for j in range(0,len(yData[i])):
            newyData[i,j]=yData[i][j]
            newyData[i+len(yData),j]=yData[i][j]
            newyData[i+len(yData),j+len(yData[i])]=yData[i][j]
            newyData[i,j+len(yData[i])]=yData[i][j]
#    newyData=zeros((len(yData)*2,len(yData[0])*2))+0.0

    contour(yData)
    print "hello",newyData
    savefig("contour.png",dpi=300)
    print "goodby",newyData
    for counter in range(0,len(yData)):
        for counter2 in range(0,len(yData[counter])):
            print yData[counter][counter2],
            print " ",
        print
