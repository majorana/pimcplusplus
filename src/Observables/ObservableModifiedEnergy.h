#ifndef OBSERVABLE_MODIFIED_ENERGY_H
#define OBSERVABLE_MODIFIED_ENERGY_H

#include "ObservableBase.h"

class ModifiedEnergyClass : public ObservableClass
{

private:
  double TotalSum, KineticSum, dUShortSum, dULongSum, NodeSum, VShortSum, VLongSum, TIP5PSum, ST2Sum, RotKinSum, P2RotKinSum;

  ObservableDouble TotalVar, KineticVar, dUShortVar, dULongVar, NodeVar,
    VShortVar, VLongVar, TIP5PVar, ST2Var, RotKinVar, P2RotKinVar;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);
  ModifiedEnergyClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      TotalVar  ("Total",  IOSection,myPathData.Path.Communicator),
      KineticVar("Kinetic",IOSection,myPathData.Path.Communicator), 
      dUShortVar("dUShort",IOSection,myPathData.Path.Communicator), 
      dULongVar ("dULong", IOSection,myPathData.Path.Communicator), 
      NodeVar   ("Node",   IOSection,myPathData.Path.Communicator), 
      VShortVar ("VShort",IOSection,myPathData.Path.Communicator), 
      VLongVar  ("VLong",IOSection,myPathData.Path.Communicator),
      TIP5PVar  ("TIP5P",IOSection,myPathData.Path.Communicator),
      ST2Var  ("ST2",IOSection,myPathData.Path.Communicator),
      RotKinVar  ("RotKin",IOSection,myPathData.Path.Communicator),
      P2RotKinVar  ("P2RotKin",IOSection,myPathData.Path.Communicator)
  {
    TotalSum   = 0.0;
    KineticSum = 0.0;
    dUShortSum = 0.0;
    dULongSum  = 0.0;
    NodeSum    = 0.0;
    VShortSum  = 0.0;
    VLongSum   = 0.0;
    TIP5PSum = 0.0;
    ST2Sum = 0.0;
    RotKinSum = 0.0;
    P2RotKinSum = 0.0;
    NumSamples = 0;
    TimesCalled=0;
  }
};

/*
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
*/

#endif 
