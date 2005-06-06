#!/usr/bin/python
from IO import *
import sys
import os
import stats
import numarray
from matplotlib.matlab import *
from HTMLgen import *
#import povexport
#from visual import *

 
def IsMonotonic (x):
     isMono = True
#     for i in range(0,x.size()-2):
     for i in range(0,len(x)-2):
          isMono = isMono and (x[i+1] > x[i])
     return isMono
     

#takes a vector and returns the unweighted average
def Avg (x):
     if x[0] == None:
          return None
     else:
          return sum(x)/len(x)

#Takes a vector of means and erros and returns the weighted average and error
def WeightedAvg (means, errors):
     if (errors[0] != 0.0):
          weights = map (lambda x: 1.0/(x*x), errors)
          norm = 1.0/sum(weights)
          weights = map(lambda x: x*norm, weights)
          avg = 0.0
          error2 = 0.0
          for i in range (0,len(means)):
               avg = avg + means[i]*weights[i]
               error2 = error2 + weights[i]**2*errors[i]*errors[i]
          return (avg, math.sqrt(error2))
     else:
          return (Avg(means), 0.0)



#takes a vector of vectors and returns a vector of unweighted averages
def VecAvg (x):
     if x[0] == None:
          return None
     else:
          return map(lambda y:sum(y)/len(y),x)

# Takes a list of 2D arrays and returns the unweighted average of the
# last row.
def AvgLastVec (data):
     if data==None or data[0]==None:
          return 0
     s = data[0][-1]
     for i in range(0,len(data)):
          s = s+data[i][-1]
     return s/len(data)



def ProduceCorrelationPicture(x,y,fileBase,hlabel,vlabel):
     clf()
     if (IsMonotonic(x)):
          plot(x, y)
     else:
          plot(x, y, 'o')
     h1=xlabel(hlabel)
     set(h1,"FontSize",20)
     v1=ylabel(vlabel)
     set(v1,"FontSize",20)
     labels = get(gca(), 'xticklabels')
     set(labels, 'fontsize', 16)
     labels = get(gca(), 'yticklabels')
     set(labels, 'fontsize', 16)
     currAxis=axis()
     currAxis[0]=cutoff
     axis(currAxis)
     savefig(fileBase+".png",dpi=60)
     savefig(fileBase+".ps")
     myImg=Image(fileBase+".png")
     return myImg



def WriteAsciiFile (asciiFileName,x,y):
     asciiFile = open (asciiFileName, "w")
     n = len(x)
     for i in range(0,n):
##          asciiFile.write(repr(x[i]) + ' ' + repr(data[-1,i]) +'\n')
          asciiFile.write('%20.16e %20.16e\n' % (x[i], y[i]))
     asciiFile.close()
     return

def BuildTable():
     myTable=Table()
     myTable.border=0
     myTable.width='40%'
     myTable.column1_align='center'
     myTable.cell_align='center'
     return myTable

def ProcessCorrelationSection(infiles,doc,currNum):
     #acquire data about the correlation section
     sectionName=infiles.GetName()
     hlabel=infiles.ReadVar("xlabel")[0]
     vlabel=infiles.ReadVar("ylabel")[0]
     data=infiles.ReadVar("y")
     if (data==[None]):
          return currNum
     isCumulative=infiles.ReadVar("Cumulative")
     print "My cumulativeness is",isCumulative
     if (isCumulative=="false"):
          y = Avg(VecAvg(data))
     else:
          y = AvgLastVec(data)
     x=infiles.ReadVar("x")[0]
     if (x==None):
          return currNum
     description=infiles.ReadVar("Description")[0]

     currNum=currNum+1
     baseName=sectionName+repr(currNum)

##Produce Image
     myImg=ProduceCorrelationPicture(x, y,baseName,hlabel,vlabel)
##Produce Ascii file
     asciiFileName = baseName + '.dat'
     WriteAsciiFile(asciiFileName,x,y)
     psFileName=baseName+'.ps'

##create table with file names in it
     fileTable=BuildTable()
     fileTable.body = [[Href(psFileName,'PostScript'), Href(asciiFileName,'ASCII data')]]
##Write things to document
     doc.append(Heading(1,sectionName))
     doc.append(Heading(4,description))
     doc.append(myImg)
     doc.append(BR())
     doc.append(fileTable)
     return currNum

def compare(a):
     return a[0]

def ProcessStructureFactor(infiles,doc,currNum):
     #acquire data about the structure factor
     sectionName=infiles.GetName()
     hlabel=infiles.ReadVar("xlabel")[0]
     vlabel=infiles.ReadVar("ylabel")[0]
     data=infiles.ReadVar("y")
     if (data==[None]):
          return currNum
     y = AvgLastVec(data)

     x=infiles.ReadVar("x")[0]
     if (x==None):
          return currNum
     toSort=[]
     for counter in range(0,len(x)):
          toSort.append((x[counter],y[counter]))
     toSort.sort()
     for counter in range(0,len(x)):
          x[counter]=toSort[counter][0]
          y[counter]=toSort[counter][1]

     xNew=[]
     yNew=[]
     counter=1
#     xNew.append(x[0])
#     yNew.append(y[0])
     while (counter<len(x)):
          totalY=y[counter]
          numY=1
          counter=counter+1
          while (counter<len(x) and x[counter]-x[counter-1]<1e-10):
               totalY=totalY+y[counter]
               numY=numY+1
               counter=counter+1
          xNew.append(x[counter-1])
          yNew.append(totalY/(numY+0.0))
     description=infiles.ReadVar("Description")[0]
     
     currNum=currNum+1
     baseName=sectionName+repr(currNum)

##Produce Image
     myImg=ProduceCorrelationPicture(x, y,baseName,hlabel,vlabel)
     myImg_avg=ProduceCorrelationPicture(xNew,yNew,baseName+"_avg",hlabel,vlabel)
##Produce Ascii file
     asciiFileName = baseName + '.dat'
     asciiFileName_avg = baseName + '_avg.dat'
     WriteAsciiFile(asciiFileName,x,y)
     WriteAsciiFile(asciiFileName_avg,xNew,yNew)
     psFileName=baseName+'.ps'
     psFileName_avg=baseName+'_avg.ps'

##create table with file names in it
     fileTable=BuildTable()
     fileTable.body = [[Href(psFileName,'PostScript'), Href(asciiFileName,'ASCII data')]]
     fileTable_avg=BuildTable()
     fileTable_avg.body = [[Href(psFileName_avg,'PostScript'), Href(asciiFileName_avg,'ASCII data')]]
     
##Write things to document
     doc.append(Heading(1,sectionName))
     doc.append(Heading(4,description))
     doc.append(myImg)
     doc.append(BR())
     doc.append(fileTable)
     doc.append(myImg_avg)
     doc.append(BR())
     doc.append(fileTable_avg)
     return currNum


def LongRangeImage(basename,r,long,short,myTitle,labelY):
     clf()
     hold ("off")
     l1 = plot (r, long)
     a = axis()
     hold ("on")
     l2 = plot (r[1:-1], long[1:-1]+short[1:-1], 'r-')
     set (l1, 'linewidth', 2);
     set (l2, 'linewidth', 2);
     h1 = xlabel("r")
     axis(a)
     set (h1, "FontSize", 20)
     h2 = ylabel (labelY)
     set (h2, "FontSize", 20)
     labels = get(gca(), 'xticklabels')
     set(labels, 'fontsize', 16)
     labels = get(gca(), 'yticklabels')
     set(labels, 'fontsize', 16)
     h3 = legend ('Ulong')
     h4 = title (myTitle)
     set (h4, "FontSize", 20)
     savefig(basename+".png",dpi=60)
     return Image(basename+".png")


def ProcessLongRangeAction(infiles):
     doc = SeriesDocument()
     doc.logo=""
     doc.author="Ken and Bryan"
     doc.email="esler@uiuc.edu and bkclark@uiuc.edu"
     doc.banner=("http://esler.physics.uiuc.edu/pimcLogo.png")
     doc.place_nav_buttons=0
     doc.header()

     infiles.OpenSection("U");
     r = infiles.ReadVar ("r")[0]
     numActions = infiles.CountSections2("PairAction")
     for i in range (0,numActions):
          infiles.OpenSection2("PairAction", i)
          ptcl1 = infiles.ReadVar("Particle1")[0]
          ptcl2 = infiles.ReadVar("Particle2")[0]
          print "Found long range breakup for ", ptcl1, " and ", ptcl2
          numLevels = infiles.CountSections2("Level")
          for level in range(0,numLevels):
               basename = "LongRangeU(" + ptcl1 + "," + ptcl2 + ")_" +repr(level)
               infiles.OpenSection2("Level", level)
               Ulong  = infiles.ReadVar("Ulong")[0]
               Ushort = infiles.ReadVar("Ushort")[0]
               myTitle = "Ulong for (" + ptcl1 + "," + ptcl2 +") level " + repr(level)
               doc.append(LongRangeImage(basename, r, Ulong, Ushort, myTitle, 'U(r)'))
               infiles.CloseSection() # "Level"
          infiles.CloseSection()      # "PairAction"
     infiles.CloseSection()           # "U"

     infiles.OpenSection("dU");
     r = infiles.ReadVar ("r")[0]
     numActions = infiles.CountSections2("PairAction")
     for i in range (0,numActions):
          infiles.OpenSection2("PairAction", i)
          ptcl1 = infiles.ReadVar("Particle1")[0]
          ptcl2 = infiles.ReadVar("Particle2")[0]
          print "Found long range breakup for ", ptcl1, " and ", ptcl2
          numLevels = infiles.CountSections2("Level")
          for level in range(0,numLevels):
               basename = "LongRangedU(" + ptcl1 + "," + ptcl2 + ")_" +repr(level)
               infiles.OpenSection2("Level", level)
               dUlong  = infiles.ReadVar("dUlong")[0]
               dUshort = infiles.ReadVar("dUshort")[0]
               myTitle = "dUlong for (" + ptcl1 + "," + ptcl2 +") level " + repr(level)
               doc.append(LongRangeImage(basename, r, dUlong, dUshort, myTitle, 'dU(r)'))
               infiles.CloseSection() # "Level"
          infiles.CloseSection()      # "PairAction"
     infiles.CloseSection()           # "dU"

     
     doc.write("longrange.html")               
     return Href("longrange.html", "Long Range Breakups")

     
##      numVars=infiles.CountVars()
##      for counter in range(0,numVars):
##           data=infiles.ReadVar(counter)
##           if type(data)==numarray.numarraycore.NumArray:

##                varName=infiles.GetVarName(counter)
##                doc.append(Name(sectionName+varName+repr(currNum)))
##                doc.append(Heading(2,varName))
##                myImg=ProduceCorrelationPicture(data[-1],varName+repr(currNum),'r',varName)

def ProduceTracePicture(data,fileBase,hlabel,vlabel,myTitle=''):
#produce scalar trace image and ps with data as a function of index
    clf()
    x=fromfunction(lambda i:i,(len(data),))
    plot(x,data)
    h1=xlabel(hlabel)
    set(h1,"FontSize",20)
    v1=ylabel(vlabel)
    set(v1,"FontSize",20)
    if len(myTitle) != 0:
         t1=title(myTitle)
         set(t1,"FontSize",20)
    labels = get(gca(), 'xticklabels')
    set(labels, 'fontsize', 16)
    labels = get(gca(), 'yticklabels')
    set(labels, 'fontsize', 16)
    savefig(fileBase+".png",dpi=45)
    savefig(fileBase+".ps")
    myImg=Image(fileBase+".png")

#produce asciiFile
    asciiFileName = fileBase + '.dat'
    asciiFile = open (asciiFileName, "w")
    n = len(data)
    for i in range(0,len(data)):
         asciiFile.write('%20.16e\n' % data[i])
    asciiFile.close()

#build table with image, ps, and ascii file in it
    fileTable=BuildTable()
    fileTable.width='100%'
    fileTable.body= [[Href(fileBase+".ps",'Postscript'),Href(asciiFileName,'ASCII data')]]
    myTable=BuildTable()
    myTable.body=[[myImg]]
    myTable.body.append([fileTable])
    return myTable


def MeanErrorString (mean, error):
     if (mean!=0.0):
          meanDigits = math.floor(math.log(abs(mean))/math.log(10))
     else:
          meanDigits=2
     if (error!=0.0):
          rightDigits = -math.floor(math.log(error)/math.log(10))+1
     else:
          rightDigits=2
     if (rightDigits < 0):
          rightDigits = 0
     formatstr = '%1.' + '%d' % rightDigits + 'f'
     meanstr  = formatstr % mean
     errorstr = formatstr % error
     return (meanstr, errorstr)

def BuildScalarTracePage(data,baseName,varName):
     doc=SimpleDocument()
     #wrong
     picTable = Table()
     picTable.width = 350*len(data)
     picTable.border = 0
     picTable.body = [[]]
     row = []
     for i in range(0,len(data)):
          d = data[i];
          row.append(ProduceTracePicture(d, baseName+"_"+repr(i), 'Blocks', varName,"Proc "+repr(i)))
     picTable.body.append(row)
     doc.append(picTable)
     doc.write(baseName+'.html')
     return baseName+'.html'                    
#     map(lambda x: doc.append(ProduceTracePicture(x,baseName,'Blocks',varName)),data)
#     doc.write(baseName+'.html')
#     return baseName+'.html'
          
def ProcessScalarSection(infiles,doc,currNum):
     sectionName=infiles.GetName()
     doc.append(Heading(1,sectionName))
     toAddList=[]
    #put description
     myTable=Table()
     myTable.body=[['','Mean','Error','Variance', 'Kappa']]
     myTable.width='50%'
     numVars=infiles.CountVars()
#     print "Num vars is ",numVars
     for counter in range(0,numVars):
          data = infiles.ReadVar(counter)
##          print "data is ",data
##          print type(data[0])
          data=data[scalarCutoff:len(data)]
          if type(data[0])==numarray.numarraycore.NumArray:
               currNum=currNum+1
               varName=infiles.GetVarName(counter)
               baseName = varName+repr(currNum)
               toAddList.append(Name(sectionName+varName+repr(currNum)))
               toAddList.append(Heading(2,varName))
               pageName=BuildScalarTracePage(data,baseName,varName)
               print "My page name is",pageName
               myFrame=IFrame()
#               myFrame.scrolling='no'
               myFrame.src=pageName
               myFrame.width="100%"
               myFrame.height="375"
#               myImg=ProduceTracePicture(data[0], baseName,'Blocks',varName)
               toAddList.append(myFrame)

##   Write ASCII data to a file
               asciiFileName = baseName + '.dat'
               asciiFile = open (asciiFileName, "w")
               n = len(data[0])
               for i in range(0,n):
                    asciiFile.write('%20.16e\n' % data[0][i])
               asciiFile.close()
               psFileName=baseName+'.ps'
               doc.append(BR())
               fileTable = Table()
               fileTable.body = [[Href(psFileName,'PostScript'), Href(asciiFileName,'ASCII data')]]
               fileTable.border=0
               fileTable.width='40%'
               fileTable.column1_align='center'
               fileTable.cell_align='center'
               toAddList.append(fileTable)
               meanlist = []
               varlist = []
               errorlist = []
               kappalist = []
               for d in data:
                    (mean,var,error,kappa)=stats.Stats(d)
                    meanlist.append(mean)
                    errorlist.append(error)
                    varlist.append(var)
                    kappalist.append(kappa)
               (mean,error) = WeightedAvg(meanlist, errorlist)
               print repr(mean) + "+/-" + repr(error)
                    
               (meanstr, errorstr) = MeanErrorString (mean, error)
               myTable.body.append([Href("#"+sectionName+varName+repr(currNum),varName),\
                                    meanstr,errorstr, '%1.2e' % var ,'%1.2f' % kappa])
     doc.append(myTable)
     for counter in range(0,len(toAddList)):
          doc.append(toAddList[counter])
     return currNum


def ProcessSystemInfo(infiles):
     tau=infiles.ReadVar("tau")[0]
     box=infiles.ReadVar("Box")[0]
     numTimeSlices=infiles.ReadVar("NumTimeSlices")[0]
     beta=tau*numTimeSlices
     temp=1.0/beta
     systemTable=Table("System")
     systemTable.body=[["tau",repr(tau)]]
     systemTable.body.append(["# of Slices",repr(numTimeSlices)])
     systemTable.body.append(["beta", '%1.2f' % beta])
     systemTable.body.append(["temperature", '%1.4e' % temp])
     systemTable.body.append(["Box","[ "+'%1.2f' % box[0]+",  "+ '%1.2f' % box[1]+",  "+'%1.2f' % box[2]+" ]"])
     speciesTable=Table("Species")
     speciesTable.body=[]
     speciesTable.body.append(["Name","NumParticles","lambda","Type"])
     numSections=infiles.CountSections2("Species")
     print 'NumSpecies = ', numSections
     for spec in range(0,numSections):
          infiles.OpenSection2("Species",spec)
          name=infiles.ReadVar("Name")[0]
          numPtcl=infiles.ReadVar("NumParticles")[0]
          lambdam = infiles.ReadVar("lambda")[0]
          type=infiles.ReadVar("ParticleType")[0]
          speciesTable.body.append([name,numPtcl,lambdam,type])
          infiles.CloseSection()
     totalTable=Table()
     totalTable.body=[]
     totalTable.border=0
     totalTable.body.append([systemTable])
     totalTable.body.append([speciesTable])
     return totalTable
     
                       
                       


def ProcessRunInfo(infiles):
     myTable=Table("Run Information")
     myTable.body=[]
     myTable.width='40%'
     numVars=infiles.CountVars()
     for counter in range(0,numVars):
          data=infiles.ReadVar(counter)[0]
          varName=infiles.GetVarName(counter)
          myTable.body.append([varName,data])
     return myTable



def ProcessTopTable(doc,infiles):
     largeTable=Table()
     largeTable.border=0
     infiles.OpenSection("RunInfo")
     runTable=ProcessRunInfo(infiles)
     infiles.CloseSection()
     infiles.OpenSection("System")
     speciesTable=ProcessSystemInfo(infiles)
     infiles.CloseSection()
     
     # Write the input file to the output directory
     InputFile = infiles.ReadVar("InputFile")[0]
     file = open ("pimc.in", "w")
     file.write(InputFile)
     file.close()

     runTable.body.append (["Input File", Href("pimc.in", "pimc.in")])
     largeTable.body.append([runTable,speciesTable])
     doc.append(largeTable)
     doc.append(HR())
     return doc

def ProcessMove(doc,infiles):
     doc.append(Heading(1,"Moves"))
     myTable=Table()
     myTable.body=[['Moves','Acceptance']]
     myTable.width='100%'
     infiles.OpenSection("Moves")
     numMoves=infiles.CountSections()
     for i in range (0, numMoves):
          infiles.OpenSection(i)
          name=infiles.GetName()
##          print "they are ",infiles.ReadVar("AcceptRatio")
          ar = VecAvg(infiles.ReadVar("AcceptRatio"))
          if (ar!=None):
               totAccept=ar[0]
               numAccept=1
               for counter in range(1,len(ar)):
                    totAccept=totAccept+ar[counter]
                    numAccept=numAccept+1
               myTable.body.append([name+" ",totAccept/numAccept])
          else:
               myTable.body.append([name,"No acceptance available"])
          numStages=infiles.CountSections()
          print "The number of stages is",numStages
          if numStages!=0:
               stageTable=Table()
               stageTable.body=[['Stages','Acceptance',"Total Attempts"]]
               stageTable.width='100%'
          for i in range (0, numStages):
               infiles.OpenSection(i)
               name=infiles.GetName()
               ar = VecAvg(infiles.ReadVar("AcceptRatio"))
               if (ar!=None):
                    totAccept=ar[0]
                    numAccept=1
                    for counter in range(1,len(ar)):
                         totAccept=totAccept+ar[counter]
                         numAccept=numAccept+1
                    stageTable.body.append([name+repr(i)+" ",totAccept/numAccept,""])
               else:
                    ar=VecAvg(infiles.ReadVar("Acceptance Ratio"))
                    if (ar!=None):
                         totAccept=ar[0]
                         numAccept=1
                         for counter in range(1,len(ar)):
                              totAccept=totAccept+ar[counter]
                              numAccept=numAccept+1
                         ar=VecAvg(infiles.ReadVar("Perms Tried"))
                         if (ar!=None):
                              totPerm=ar[0]
                              numPerm=1
                              for counter in range(1,len(ar)):
                                   totPerm=totPerm+ar[counter]
                                   numPerm=numPerm+1
                              stageTable.body.append(["PermutationStage"+" ",totAccept/numAccept,(totPerm+0.0)/(numPerm+0.0)])
                    else:
                         myTable.body.append([name,"No acceptance available"])
               infiles.CloseSection()
          if numStages!=0:
               myTable.body.append(["StageInfo",stageTable])
          infiles.CloseSection()
     infiles.CloseSection() # "Moves"
     doc.append(myTable)



basename = sys.argv[1]
infiles = IOSectionClassList()
infiles.OpenFiles(basename);
      
print 'Found ' +repr(infiles.len()) + ' output files.'






dirName=basename 
cutoff=None
StartCut = None
if (os.access(dirName+"/.pref",os.F_OK)):
     print dirName+"/.pref"
     prefFile=IOSectionClass()
     prefFile.OpenFile(dirName+"/.pref")
     prefFile.ReadVar("StartCut")
     print "StartCut = ", StartCut
     print prefFile.ReadVar("cutoff")
     cutoff=prefFile.ReadVar("cutoff")
     print "cutoff = ", cutoff
#     prefFile.CloseFile() 
if cutoff==None:
     cutoff=0
if scalarCutoff==None:
     scalarCutoff=0
if not(os.access(dirName,os.F_OK)):
     os.mkdir(dirName)
os.chdir(dirName)
doc=SeriesDocument()
ProcessTopTable(doc,infiles)

ProcessMove(doc,infiles)

currNum=0
infiles.OpenSection("Observables")
numSections=infiles.CountSections()




#print "The number of sections is ",numSections
for counter in range(0,numSections):
     infiles.OpenSection(counter)
     myName= infiles.GetName()
     myType=infiles.ReadVar("Type")[0]
     if myName=="StructureFactor":
          currNum=ProcessStructureFactor(infiles,doc,currNum)
          doc.append(HR())
     elif myType=="Scalar":
          currNum=ProcessScalarSection(infiles,doc,currNum)
          doc.append(HR())
     elif myType=="CorrelationFunction":
          currNum=ProcessCorrelationSection(infiles,doc,currNum)
          doc.append(HR()) 
     else:
          a=5
     infiles.CloseSection()
infiles.CloseSection() # "Observables"

if (infiles.CountSections2("LongRangeAction") > 0):
     infiles.OpenSection("LongRangeAction")
     LRsection = ProcessLongRangeAction (infiles)
     doc.append(LRsection)
     infiles.CloseSection() # "LongRangeAction"



#myFrame=IFrame("index.html","blah")
#myFrame.src="hi"
#doc.append(myFrame)
doc.logo=""
doc.author="Ken and Bryan"
doc.email="esler@uiuc.edu and bkclark@uiuc.edu"
doc.banner=("http://esler.physics.uiuc.edu/pimcLogo.png")
doc.place_nav_buttons=0
doc.header()

doc.write("index.html")
