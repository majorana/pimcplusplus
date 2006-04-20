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
     return (totalTable,tau,numTimeSlices)
     

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
     (speciesTable,tau,numTimeSlices)=ProcessSystemInfo(infiles)
     #speciesTable="broken"
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
     return (doc,tau,numTimeSlices)
