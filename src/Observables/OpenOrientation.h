#ifndef OPEN_ORIENTATION_H
#define OPEN_ORIENTATION_H

#include "ObservableBase.h"

class OpenOrientationClass : public ObservableClass
{

private:
  //Must be initialized
  double R2;
  double Z2;
  double R2OverZ2;

  ObservableDouble R2Var;
  ObservableDouble Z2Var;
  ObservableDouble R2OverZ2Var;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  void WriteInfo();
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& in);
  OpenOrientationClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
    R2Var("R2",IOSection,myPathData.Path.Communicator),
    Z2Var("Z2",IOSection,myPathData.Path.Communicator),
    R2OverZ2Var("R2OverZ2",IOSection,myPathData.Path.Communicator)
    {
      NumSamples = 0;
      TimesCalled=0;
    }
};


#endif 
