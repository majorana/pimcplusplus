from IO import *
import sys
from numpy import *

class Squarer2HDFParser:
  def __init__(self,filename):
    self.f = open(filename,'r')
    self.basename = filename[:-2]
    outfilename = self.basename + 'h5'
    self.a=IOSectionClass()
    self.a.NewFile(outfilename)
    self.numFits = 0
    self.SpitPotential = True
    self.Ucounter = 0 
    self.dUcounter = 0

  def Done(self):
    print "Done."
    self.f.close()
    self.a.FlushFile()
    self.a.CloseFile()

  def next(self):
    w = ''
    c = self.f.read(1)
    # check for EOF
    if(c == ''):
      print "ENCOUNTERED EOF"
      sys.exit()
      return w
    empty = True
    while(empty):
      while(c!=' ' and c!='\t' and c!='\n' and c!=''):
        w += c
        c = self.f.read(1)
        empty = False
      if(empty):
        c = self.f.read(1)
    #print "parsed",w
    return w

  def find(self, target):
    s = self.next()
    while(s != target):
      s = self.next()

  def ProcessSquarerInfo(self):
    print "ProcessSquarerInfo"
    ### collect squarer run info 
    self.a.NewSection("Squarer")
    self.a.NewSection("Units")
    g = self.next()
    check(g,'UNITS')
    self.a.WriteVar("Temp",self.next())
    self.a.WriteVar("Length",self.next())
    self.a.CloseSection()
    g = self.next() 
    check(g,'TYPE')
    self.a.NewSection("Type")
    self.a.WriteVar("Species",self.next())
    self.a.WriteVar("Lambda",float(self.next()))
    self.a.CloseSection()
    g = self.next() 
    check(g,'TYPE')
    self.a.NewSection("Type")
    self.a.WriteVar("Species",self.next())
    self.a.WriteVar("Lambda",float(self.next()))
    self.a.CloseSection()
    
    ### get important stats for remainder of read
    self.find("SQUARER")
    self.next()
    self.next()
    self.next()
    self.numFits = int(self.next())
    self.a.WriteVar("NumFits",self.numFits)
    self.a.CloseSection()
    #self.a.WriteVar("NumFits",self.numFits)
    self.a.FlushFile()
  ### end SquarerInfo ###

  ### get potential
  def ProcessPotential(self):
    print "ProcessPotential"
    self.find("RANK")
    self.next()
    size = int(self.next())
    self.find("BEGIN")
    self.next()
    self.next()
    u = zeros(size) + 0.
    for i in range(0,size):
      u[i] = float(self.next())
    self.a.NewSection("Potential")
    self.a.WriteVar("Data",u)
    self.a.CloseSection()
    self.a.FlushFile()
    # dump potential to ASCII
    if(self.SpitPotential):
      pout = open(self.basename + 'potential.dat', 'w')
      for i in range(0,len(u)):
        pout.write(str(u[i]) + '\n')
      pout.close()
  ### end ProcessPotential ###

  def ProcessU(self):
    print "ProcessU"
    self.find("RANK")
    rnk = int(self.next())
    check(rnk,3)
    numPts = int(self.next())
    numUkj = int(self.next())
    numTau = int(self.next())
  
    self.find("GRID")
    self.next()
    gridType = self.next()
    start = float(self.next())
    end = float(self.next())
    self.find("GRID")
    self.next()
    checkgrid = self.next()
    check(checkgrid, "LOG")
    loTau = float(self.next())
    hiTau = float(self.next())
   
    self.find("BEGIN")
    self.next()
    self.next()
    self.next()
    self.next()
    NMax = int(self.next())
    self.next()
    self.next()
    derv = int(self.next())
    #check(derv,order)
  
    # load array from .dm file
    Ukj = zeros([numPts, numUkj, numTau]) + 0.
    for cT in range(0,numTau):
      for cU in range(0,numUkj):
        for cG in range(0,numPts):
          Ukj[cG, cU, cT] = float(self.next())
  
    Taus = zeros(numTau) + 0.
    tau0 = loTau/2
    for t in range(0,numTau):
      tau0 *= 2.
      Taus[t] = tau0
   
    print NMax,derv
    SectionTitle = 'Ukj' + str(self.Ucounter)
    self.Ucounter += 1
    #SectionTitle = ''
    #print "Sectiontitle reset"
    #if(order == 1):
    #  self.Ucounter += 1
    #  print "setting sectiontitle",self.SectionTitle,"to"
    #  self.SectionTitle = 'Ukj' + str(self.Ucounter)
    #  print order,self.SectionTitle
    #elif(order == 2):
    #  self.dUcounter += 1
    #  print "setting sectiontitle",self.SectionTitle,"to"
    #  self.SectionTitle = 'dUkj_dBeta' + str(self.dUcounter)
    #  print order,self.SectionTitle
    #else:
    #  print "What's going on? order is",order
    #  sys.exit()
    #print order,self.SectionTitle
    print "Going to write section",SectionTitle
    self.a.NewSection(SectionTitle)
    self.a.NewSection("Grid")
    self.a.WriteVar("NumGridPoints",numPts)
    self.a.WriteVar("Type",gridType)
    self.a.WriteVar("Start",start)
    self.a.WriteVar("End",end)
    self.a.CloseSection()
    self.a.WriteVar("NumUkj",numUkj)
    self.a.WriteVar("NumTau",numTau)
    self.a.WriteVar("NMax",NMax)
    self.a.WriteVar("Derv",derv)
    self.a.WriteVar("Rank",rnk)
    self.a.WriteVar("Taus",Taus)
    self.a.WriteVar("Data",Ukj)
    self.a.CloseSection()
    self.a.FlushFile()
  ### end ProcessU ###

  def ProcessdU_dBeta(self):
    print "ProcessdU_dBeta"
    self.find("RANK")
    rnk = int(self.next())
    check(rnk,3)
    numPts = int(self.next())
    numUkj = int(self.next())
    numTau = int(self.next())
  
    self.find("GRID")
    self.next()
    gridType = self.next()
    start = float(self.next())
    end = float(self.next())
    self.find("GRID")
    self.next()
    checkgrid = self.next()
    check(checkgrid, "LOG")
    loTau = float(self.next())
    hiTau = float(self.next())
   
    self.find("BEGIN")
    self.next()
    self.next()
    self.next()
    self.next()
    NMax = int(self.next())
    self.next()
    self.next()
    derv = int(self.next())
    #check(derv,order)
  
    # load array from .dm file
    Ukj = zeros([numPts, numUkj, numTau]) + 0.
    for cT in range(0,numTau):
      for cU in range(0,numUkj):
        for cG in range(0,numPts):
          Ukj[cG, cU, cT] = float(self.next())
  
    Taus = zeros(numTau) + 0.
    tau0 = loTau/2
    for t in range(0,numTau):
      tau0 *= 2.
      Taus[t] = tau0
   
    print NMax,derv
    SectionTitle = 'dUkjdBeta' + str(self.dUcounter)
    self.dUcounter += 1
    #SectionTitle = ''
    #print "Sectiontitle reset"
    #if(order == 1):
    #  self.Ucounter += 1
    #  print "setting sectiontitle",self.SectionTitle,"to"
    #  self.SectionTitle = 'Ukj' + str(self.Ucounter)
    #  print order,self.SectionTitle
    #elif(order == 2):
    #  self.dUcounter += 1
    #  print "setting sectiontitle",self.SectionTitle,"to"
    #  self.SectionTitle = 'dUkj_dBeta' + str(self.dUcounter)
    #  print order,self.SectionTitle
    #else:
    #  print "What's going on? order is",order
    #  sys.exit()
    #print order,self.SectionTitle
    print "Going to write section",SectionTitle
    self.a.NewSection(SectionTitle)
    self.a.NewSection("Grid")
    self.a.WriteVar("NumGridPoints",numPts)
    self.a.WriteVar("Type",gridType)
    self.a.WriteVar("Start",start)
    self.a.WriteVar("End",end)
    self.a.CloseSection()
    self.a.WriteVar("NumUkj",numUkj)
    self.a.WriteVar("NumTau",numTau)
    self.a.WriteVar("NMax",NMax)
    self.a.WriteVar("Derv",derv)
    self.a.WriteVar("Rank",rnk)
    self.a.WriteVar("Taus",Taus)
    self.a.WriteVar("Data",Ukj)
    self.a.CloseSection()
    self.a.FlushFile()
  ### end ProcessdU_dBeta ###
  
  def ProcessSampling(self):
    print "ProcessSampling",
    self.find("RANK")
    rnk = int(self.next())
    check(rnk,3)
    numPts = int(self.next())
    numUkj = int(self.next())
    numTau = int(self.next())
  
    self.find("GRID")
    self.next()
    checkgrid = self.next()
    check(checkgrid, "LOG")
    loTau = float(self.next())
    hiTau = float(self.next())
    
    self.find("BEGIN")
    self.next()
    self.next()
    derv = int(self.next())
    print derv
  
    # load array from .dm file
    Ukj = zeros([numPts, numUkj, numTau]) + 0.
    for cT in range(0,numTau):
      for cU in range(0,numUkj):
        for cG in range(0,numPts):
          #print derv,cT,"/",numTau,cU,"/",numUkj,cG,"/",numPts,
          Ukj[cG, cU, cT] = float(self.next())
  
    Taus = zeros(numTau) + 0.
    tau0 = loTau/2
    for t in range(0,numTau):
      tau0 *= 2.
      Taus[t] = tau0
    
    self.a.NewSection("Sampling")
    self.a.WriteVar("Taus",Taus)
    self.a.WriteVar("NumUkj",numUkj)
    self.a.WriteVar("NumTau",numTau)
    self.a.WriteVar("Derv",derv)
    self.a.WriteVar("Data",Ukj)
    self.a.CloseSection()
    self.a.FlushFile()
  ### end ProcessSampling ###

### end class Squarer2HDFParser ###

def check(s, sCheck):
  if(s != sCheck):
    print "MISMATCH: expected",sCheck,"got",s
    return(False)
  return(True);

### main ###
if(len(sys.argv)!=2):
  print "Usage: python dmparse.py squarer_output_file.dm"
  sys.exit()

sq = Squarer2HDFParser(sys.argv[1])
sq.ProcessSquarerInfo()
sq.ProcessPotential()
for sec in range(0,sq.numFits + 1):
  sq.ProcessU()
for sec in range(0,sq.numFits + 1):
  sq.ProcessdU_dBeta()
sq.ProcessSampling()
sq.ProcessSampling()
sq.Done()
