#ifndef PRESSURE_H
#define PRESSURE_H

#include "ObservableBase.h"

class PressureClass : public ObservableClass
{
protected:
  ObservableDouble PressureVar, KineticVar, ShortRangeVar, LongRangeVar;
  double Psum;
  double KineticSum, ShortRangeSum, LongRangeSum;
  int NumSamples;
  double PartitionPressure();
  double KineticPressure();
  double ShortRangePressure();
  double LongRangePressure();
public:
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass &in);
  PressureClass(PathDataClass &pathData, IOSectionClass &ioSection)
    : ObservableClass (pathData, ioSection),
      PressureVar  ("Total",      IOSection, pathData.Path.Communicator),
      KineticVar   ("Kinetic",    IOSection, pathData.Path.Communicator),
      ShortRangeVar("ShortRange", IOSection, pathData.Path.Communicator),
      LongRangeVar ("LongRange",  IOSection, pathData.Path.Communicator),
      NumSamples(0), Psum(0.0), KineticSum(0.0), ShortRangeSum(0.0),
      LongRangeSum(0.0)
  {

  }
};

#endif
