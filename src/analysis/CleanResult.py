from IO import *
import sys
import stats
#import povexport
#from visual import *
import numarray
from matplotlib.matlab import *
from HTMLgen import *

## def Gofr(infile,data,doc):
##     tempImg=ProduceCorrelationPicture(infile,data,"gofr")
##     doc.append(Heading(2,"G(r) Plot"))
##     doc.append(tempImg)
##     numVars=infile.CountVars()
##     for counter in range(0,numVars):
##          data=infile.ReadVar(counter)
##          if type(data)==numarray.numarraycore.NumArray:
##               varName=infile.GetVarName(counter)
##               doc.append(Name(sectionName+varName+repr(currNum)))
##               doc.append(Heading(2,varName))
##               myImg=ProduceCorrelationPicture(data[-1],varName+repr(currNum),'r',varName)
              


def ProduceCorrelationPicture(data,fileBase,hlabel,vlabel):
     clf()
     infile.OpenSection("grid")
     points=infile.ReadVar("Points")
     infile.CloseSection()
     plot(points,data)
     h1=xlabel(hlabel)
     set(h1,"FontSize",20)
     v1=ylabel(vlabel)
     set(v1,"FontSize",20)
     labels = get(gca(), 'xticklabels')
     set(labels, 'fontsize', 16)
     labels = get(gca(), 'yticklabels')
     set(labels, 'fontsize', 16)
     savefig(fileBase+".png",dpi=60)
     savefig(fileBase+".ps")
     myImg=Image(fileBase+".png")
     return myImg



def ProcessCorrelationSection(infile,doc,currNum):
     sectionName=infile.GetName()
     doc.append(Heading(1,infile.GetName()))
     numVars=infile.CountVars()
     for counter in range(0,numVars):
          data=infile.ReadVar(counter)
          if type(data)==numarray.numarraycore.NumArray:
               currNum=currNum+1
               varName=infile.GetVarName(counter)
               doc.append(Name(sectionName+varName+repr(currNum)))
               doc.append(Heading(2,varName))
               myImg=ProduceCorrelationPicture(data[-1],varName+repr(currNum),'r',varName)
               doc.append(myImg)
     return currNum

def ProduceTracePicture(data,fileBase,hlabel,vlabel):
    clf()
    x=fromfunction(lambda i:i,(len(data),))
    plot(x,data)
    h1=xlabel(hlabel)
    set(h1,"FontSize",20)
    v1=ylabel(vlabel)
    set(v1,"FontSize",20)
    labels = get(gca(), 'xticklabels')
    set(labels, 'fontsize', 16)
    labels = get(gca(), 'yticklabels')
    set(labels, 'fontsize', 16)
    savefig(fileBase+".png",dpi=60)
    savefig(fileBase+".ps")
    myImg=Image(fileBase+".png")
    return myImg



def ProcessScalarSection(infile,doc,currNum):
     sectionName=infile.GetName()
     doc.append(Heading(1,sectionName))
     toAddList=[]
    #put description
     myTable=Table()
     myTable.body=[['','Mean','Variance','Error','Kappa']]
     numVars=infile.CountVars()
     print "Num vars is ",numVars
     for counter in range(0,numVars):
          data=infile.ReadVar(counter)
          print type(data)
          if type(data)==numarray.numarraycore.NumArray:
               currNum=currNum+1
               varName=infile.GetVarName(counter)
               toAddList.append(Name(sectionName+varName+repr(currNum)))
               toAddList.append(Heading(2,varName))
               myImg=ProduceTracePicture(data,varName+repr(currNum),'Blocks',varName)
               toAddList.append(myImg)
               (mean,var,error,kappa)=stats.Stats(data)
               myTable.body.append([Href("#"+sectionName+varName+repr(currNum),varName),mean,var,error,kappa])
     doc.append(myTable)
     for counter in range(0,len(toAddList)):
          doc.append(toAddList[counter])
     return currNum

## def Energy(infile,doc):

##     myTable=Table("Energy Table")
##     myTable.body=[['','Mean','Variance','Error','Kappa']]
    
    
##     data=infile.ReadVar("TotalEnergy")
##     (mean,var,error,kappa)= stats.Stats(data)
##     doc.append(Heading(2,"Total Energy"))
##     myImg=ProduceTracePicture(data,"TotalEnergy.png")
##     doc.append(myImg)
    
##     myTable.body.append([Href('#totE','Total Energy'),mean,var,error,kappa])
    

##     data=infile.ReadVar("SpringEnergy")
##     (mean,var,error,kappa)= stats.Stats(data)
##     doc.append(Heading(2,"Spring Energy"))
##     myImg=ProduceTracePicture(data,"SpringEnergy.png")
##     doc.append(myImg)
##     myTable.body.append(['Spring Energy',mean,var,error,kappa])

##     data=infile.ReadVar("PotentialEnergy")
##     (mean,var,error,kappa)= stats.Stats(data)
##     doc.append(Name('potE'))
##     doc.append(Heading(2,"Potential Energy"))
##     myImg=ProduceTracePicture(data,"PotentialEnergy.png")
##     doc.append(myImg)
##     myTable.body.append([Href('#potE','Potential Energy'),mean,var,error,kappa])

               
##     data=infile.ReadVar("DBetaEnergy")
##     (mean,var,error,kappa)= stats.Stats(data)
##     doc.append(Heading(2,"du/dBeta Energy"))
##     myImg=ProduceTracePicture(data,"dBetaEnergy.png")
##     doc.append(myImg)
##     myTable.body.append(['du/dBeta Energy',mean,var,error,kappa])

##     doc.prepend(myTable)
    
    

def ProcessRunInfo(doc,infile):
     doc.append(Heading(2,"Run information"))
     myTable=Table()
     myTable.body=[]
     myTable.width=50
     numVars=infile.CountVars()
     for counter in range(0,numVars):
          data=infile.ReadVar(counter)
          varName=infile.GetVarName(counter)
          myTable.body.append([varName,data])
     doc.append(myTable)
     doc.append(HR())

infile=IOSectionClass()
infile.OpenFile(sys.argv[1])
doc=SeriesDocument()
infile.OpenSection("RunInfo")
ProcessRunInfo(doc,infile)
BuildDate=infile.ReadVar("BuildDate")
BuildTime=infile.ReadVar("BuildTime")
HostName=infile.ReadVar("HostName")
ProgramName=infile.ReadVar("ProgramName")
RunTime=infile.ReadVar("RunTime")
UserName=infile.ReadVar("UserName")
Version=infile.ReadVar("Version")
infile.CloseSection()


currNum=0
numSections=infile.CountSections()
print "The number of sections is ",numSections
for counter in range(0,numSections):
     infile.OpenSection(counter)
     myType=infile.ReadVar("Type")
     if myType=="Scalar":
          print "I'm processing"
          currNum=ProcessScalarSection(infile,doc,currNum)
     elif myType=="CorrelationFunction":
          currNum=ProcessCorrelationSection(infile,doc,currNum)
     doc.append(HR())
#     data=infile.ReadVar("gofr")
#     if data!=None:
#         Gofr(infile,data,doc)
#     elif infile.GetName()=="Energies":
#         Energy(infile,doc)
     infile.CloseSection()
doc.logo="beach.jpg"
doc.author="Ken and Bryan"
doc.email="bkclark@uiuc.edu"
doc.banner=("PICT0067.JPG",640,300)
doc.place_nav_buttons=0
doc.header()

doc.write("out.html")
