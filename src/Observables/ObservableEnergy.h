#ifndef OBSERVABLE_ENERGY_H
#define OBSERVABLE_ENERGY_H

#include "ObservableBase.h"

class EnergyClass : public ObservableClass
{

private:
  double TotalSum, KineticSum, dUShortSum, dULongSum, NodeSum, 
    VShortSum, VLongSum; //, TIP5PSum;

  ObservableDouble TotalVar, KineticVar, dUShortVar, dULongVar, NodeVar,
    VShortVar, VLongVar; //, TIP5PVar;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);
  EnergyClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      TotalVar  ("Total",  IOSection,myPathData.Path.Communicator),
      KineticVar("Kinetic",IOSection,myPathData.Path.Communicator), 
      dUShortVar("dUShort",IOSection,myPathData.Path.Communicator), 
      dULongVar ("dULong", IOSection,myPathData.Path.Communicator), 
      NodeVar   ("Node",   IOSection,myPathData.Path.Communicator), 
      VShortVar ("VShort",IOSection,myPathData.Path.Communicator), 
      VLongVar  ("VLong",IOSection,myPathData.Path.Communicator)
    //      TIP5PVar  ("TIP5P",IOSection,myPathData.Path.Communicator)
  {
    TotalSum   = 0.0;
    KineticSum = 0.0;
    dUShortSum = 0.0;
    dULongSum  = 0.0;
    NodeSum    = 0.0;
    VShortSum  = 0.0;
    VLongSum   = 0.0;
    NumSamples = 0;
    TimesCalled=0;
    //    TIP5PSum = 0;
  }
};

class EnergySignClass : public ObservableClass
{

private:
  double TotalSum, KineticSum, dUShortSum, dULongSum, NodeSum, 
    VShortSum, VLongSum;

  ObservableDouble TotalVar, KineticVar, dUShortVar, dULongVar, NodeVar,
    VShortVar, VLongVar;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);
  EnergySignClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      TotalVar  ("Total",  IOSection,myPathData.Path.Communicator),
      KineticVar("Kinetic",IOSection,myPathData.Path.Communicator), 
      dUShortVar("dUShort",IOSection,myPathData.Path.Communicator), 
      dULongVar ("dULong", IOSection,myPathData.Path.Communicator), 
      NodeVar   ("Node",   IOSection,myPathData.Path.Communicator), 
      VShortVar ("VShort",IOSection,myPathData.Path.Communicator), 
      VLongVar  ("VLong",IOSection,myPathData.Path.Communicator)
  {
    TotalSum   = 0.0;
    KineticSum = 0.0;
    dUShortSum = 0.0;
    dULongSum  = 0.0;
    NodeSum    = 0.0;
    VShortSum  = 0.0;
    VLongSum   = 0.0;
    NumSamples = 0;
    TimesCalled=0;
  }
};


#endif 
