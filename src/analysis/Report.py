#!/usr/bin/env python
from IO import *
import sys
import os
import math
import stats
import numarray
from Energy import *
from nofr import *
from PairCorrelation import *
from HTMLgen import *
from SystemInfo import *
from Josephson import *
from Coupling import *
from Vacancy import *
from PlaneDensity import *
from Pressure import *
from Langevin import *
from TimeAnalysis import *

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
     StartCut = prefFile.ReadVar("StartCut")
     print "StartCut = ", StartCut
     print prefFile.ReadVar("cutoff")
     cutoff=prefFile.ReadVar("cutoff")
     print "cutoff = ", cutoff
#     prefFile.CloseFile() 
if cutoff==None:
     cutoff=0
if not(os.access(dirName,os.F_OK)):
     os.mkdir(dirName)
os.chdir(dirName)
detailedDoc=SeriesDocument()
summaryDoc=SeriesDocument()
summaryDoc.append(Href("detail.html","Detailed HTML Page"))
summaryDoc.append(HR())


#ProcessMove(doc,infiles)


(_,tau,numTimeSlices)=ProcessTopTable(summaryDoc,infiles)
beta = tau*numTimeSlices

###########
## Moves ##
###########
infiles.OpenSection("Moves")
numSections=infiles.CountSections()

for secNum in range(0, numSections):
     infiles.OpenSection(secNum)
     moveName = infiles.GetName()
     if moveName == "Langevin":
          print "Processing Langevin move."
          ProcessLangevin(infiles, summaryDoc, detailedDoc, StartCut, beta)
          print "Done Langevin."
     infiles.CloseSection() # Current move section
infiles.CloseSection() # "Moves"


#################
## Observables ##
#################
currNum=0
infiles.OpenSection("Observables")
numSections=infiles.CountSections()

for counter in range(0,numSections):
     infiles.OpenSection(counter)
     myName= infiles.GetName()
     myType=infiles.ReadVar("Type")[0]
     print "Currently processing ",myName
     if myName=="PairCorrelation":
         ProcessPairCorrelation(infiles,summaryDoc,detailedDoc,StartCut)
     if myName=="Vacancy":
         ProcessVacancy(infiles,summaryDoc,detailedDoc,StartCut)
         summaryDoc.append(HR())
         detailedDoc.append(HR())
     elif myName=="Energy":
         ProcessEnergy(infiles,summaryDoc,detailedDoc,StartCut)
         summaryDoc.append(HR())
         detailedDoc.append(HR())
     elif myName=="nofr":
         Processnofr(infiles,summaryDoc,detailedDoc,StartCut)
         summaryDoc.append(HR())
         detailedDoc.append(HR())
     elif myName=="Pressure":
         ProcessPressure(infiles,summaryDoc,detailedDoc,StartCut)
         summaryDoc.append(HR())
         detailedDoc.append(HR())
     elif myName=="PlaneDensity":
          ProcessPlaneDensity(infiles,summaryDoc,detailedDoc,StartCut)
          summaryDoc.append(HR())
          detailedDoc.append(HR())
     elif myName=="TimeAnalysis":
          ProcessTimeAnalysis(infiles,summaryDoc,detailedDoc,StartCut)
          summaryDoc.append(HR())
          detailedDoc.append(HR())
     elif myName=="PhiK":
         ProcessJosephson(infiles,summaryDoc,detailedDoc,\
                          StartCut,tau,numTimeSlices)
#     elif myName=="Coupling":
#          ProcessCoupling(infiles,summaryDoc,detailedDoc,StartCut)
     else:
         a=5
     infiles.CloseSection()
infiles.CloseSection() # "Observables"
print "Made it here"

summaryDoc.logo=""
summaryDoc.author="Ken and Bryan"
summaryDoc.email="esler@uiuc.edu and bkclark@uiuc.edu"
summaryDoc.banner=("http://esler.physics.uiuc.edu/pimcLogo.png")
summaryDoc.place_nav_buttons=0
summaryDoc.header()

summaryDoc.write("index.html")

detailedDoc.logo=""
detailedDoc.author="Ken and Bryan"
detailedDoc.email="esler@uiuc.edu and bkclark@uiuc.edu"
detailedDoc.banner=("http://esler.physics.uiuc.edu/pimcLogo.png")
detailedDoc.place_nav_buttons=0
detailedDoc.header()

detailedDoc.write("detail.html")

        
