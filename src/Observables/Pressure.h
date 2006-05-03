#ifndef PRESSURE_H
#define PRESSURE_H

#include "ObservableBase.h"

class PressureClass : public ObservableClass
{
protected:
  ObservableDouble PressureVar, KineticVar, ShortRangeVar, LongRangeVar,
    NodeVar;
  ObservableVecDouble1 ShortRangeVecVar;
  ObservableVecDouble1 KineticVecVar;
  double Psum;
  double KineticSum, ShortRangeSum, LongRangeSum, NodeSum;
  int NumSamples;
  double PartitionPressure();
  double KineticPressure();
  double ShortRangePressure();
  double LongRangePressure();
  double NodePressure();
  Array<double,1> ShortRangePressureVec;
  Array<double,1> KineticPressureVec;
public:
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass &in);
  PressureClass(PathDataClass &pathData, IOSectionClass &ioSection)
    : ObservableClass (pathData, ioSection),
      PressureVar  ("Total",      IOSection, pathData.Path.Communicator),
      KineticVar   ("Kinetic",    IOSection, pathData.Path.Communicator),
      ShortRangeVar("ShortRange", IOSection, pathData.Path.Communicator),
      ShortRangeVecVar("ShortRangeComponents",IOSection,pathData.Path.Communicator),
      KineticVecVar("KineticComponents",IOSection,pathData.Path.Communicator),
      LongRangeVar ("LongRange",  IOSection, pathData.Path.Communicator),
      NodeVar      ("Node",       IOSection, pathData.Path.Communicator),
      NumSamples(0), Psum(0.0), KineticSum(0.0), ShortRangeSum(0.0),
      LongRangeSum(0.0), NodeSum(0.0)
  {

  }
};

#endif
