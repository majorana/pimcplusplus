#!/usr/bin/python
from IO import *
import sys
import os
import math
import stats
import numarray
from Energy import *
from PairCorrelation import *
from HTMLgen import *




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

#ProcessMove(doc,infiles)

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
     elif myName=="Energy":
         ProcessEnergy(infiles,summaryDoc,detailedDoc,StartCut)
         summaryDoc.append(HR())
         detailedDoc.append(HR())
     else:
         a=5
     infiles.CloseSection()
infiles.CloseSection() # "Observables"


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

        
