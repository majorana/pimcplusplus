#ifndef ANGULAR_H
#define ANGULAR_H

#include "ObservableBase.h"

class AngularClass : public ObservableClass
{
 private:
  int Freq,DumpFreq;
  int TotalCounts;
  int TimesCalled;
  ObservableVecDouble2 CorVar;
  
public:
  ///l x t
  Array<double,2> Correlation;
  int Species;
  void Accumulate();
  void Read(IOSectionClass& in);
  void WriteBlock();
  AngularClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    CorVar("y", IOSection, myPathData.Path.Communicator)
  {
    TimesCalled=0;
  }

};


#endif
