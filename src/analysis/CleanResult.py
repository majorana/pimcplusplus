from IO import *
import sys
import stats
#import povexport
#from visual import *
import numarray
from matplotlib.matlab import *
from HTMLgen import *


def ProduceCorrelationPicture(infile,data,varName):
    clf()
    infile.OpenSection("grid")
    points=infile.ReadVar("Points")
    infile.CloseSection()
    plot(points,data[-1])
    fileName=varName+".png"
    savefig(fileName)
    myImg=Image(fileName)
    return myImg



def Gofr(infile,data,doc):
    tempImg=ProduceCorrelationPicture(infile,data,"gofr")
    doc.append(Heading(2,"G(r) Plot"))
    doc.append(tempImg)


def ProduceTracePicture(data,fileName):
    clf()
    x=fromfunction(lambda i:i,(len(data),))
    plot(x,data)
    savefig(fileName)
    myImg=Image(fileName)
    return myImg




def Energy(infile,doc):

    myTable=Table("Energy Table")
    myTable.body=[['','Mean','Variance','Error','Kappa']]
    
    
    data=infile.ReadVar("TotalEnergy")
    (mean,var,error,kappa)= stats.Stats(data)
    doc.append(Heading(2,"Total Energy"))
    myImg=ProduceTracePicture(data,"TotalEnergy.png")
    doc.append(myImg)
    
    myTable.body.append([Href('#totE','Total Energy'),mean,var,error,kappa])
    

    data=infile.ReadVar("SpringEnergy")
    (mean,var,error,kappa)= stats.Stats(data)
    doc.append(Heading(2,"Spring Energy"))
    myImg=ProduceTracePicture(data,"SpringEnergy.png")
    doc.append(myImg)
    myTable.body.append(['Spring Energy',mean,var,error,kappa])

    data=infile.ReadVar("PotentialEnergy")
    (mean,var,error,kappa)= stats.Stats(data)
    doc.append(Name('potE'))
    doc.append(Heading(2,"Potential Energy"))
    myImg=ProduceTracePicture(data,"PotentialEnergy.png")
    doc.append(myImg)
    myTable.body.append([Href('#potE','Potential Energy'),mean,var,error,kappa])

               
    data=infile.ReadVar("DBetaEnergy")
    (mean,var,error,kappa)= stats.Stats(data)
    doc.append(Heading(2,"du/dBeta Energy"))
    myImg=ProduceTracePicture(data,"dBetaEnergy.png")
    doc.append(myImg)
    myTable.body.append(['du/dBeta Energy',mean,var,error,kappa])

    doc.prepend(myTable)
    
    
    
    


infile=IOSectionClass()
infile.OpenFile(sys.argv[1])
doc=SeriesDocument()

numSections=infile.CountSections()
print "The number of sections is ",numSections
for counter in range(0,numSections):
     infile.OpenSection(counter)
     getType=infile.ReadVar("Type")
     print "My type is ",getType
     data=infile.ReadVar("gofr")
     if data!=None:
         Gofr(infile,data,doc)
     elif infile.GetName()=="Energies":
         Energy(infile,doc)
     infile.CloseSection()
doc.logo="beach.jpg"
doc.author="Ken and Bryan"
doc.email="bkclark@uiuc.edu"
doc.banner=("PICT0067.JPG",640,300)
doc.place_nav_buttons=0
doc.header()

doc.write("out.html")
