#ifndef WEIGHT_H
#define WEIGHT_H



#include "ObservableBase.h"




class WeightClass : public ObservableClass
{

private:
  Array<double,1> Weight;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
  ObservableVecDouble1 Tot;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);
  WeightClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection),
      Tot("Total",IOSection,PathData.Path.Communicator)
  {
    Weight.resize(4);
    Weight = 0.0;
    NumSamples = 0;
    TimesCalled=0;
  }
};


#endif 
