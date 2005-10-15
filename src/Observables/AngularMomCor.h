#ifndef ANGULARMOMCOR_H
#define ANGULARMOMCOR_H

#include "ObservableBase.h"

class AngularMomCor : public ObservableClass
{
 private:
  int Freq,DumpFreq;
  int TotalCounts;
  int TimesCalled;
  ObservableVecDouble1 CorVar;
public:
  Array<double,1> Correlation;
  int Species;
  void Accumulate();
  void Read(IOSectionClass& in);
  void WriteBlock();
  AngularMomCor(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    CorVar("y", IOSection, myPathData.Path.Communicator)
  {
    TimesCalled=0;
  }

};


#endif
