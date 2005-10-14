#ifndef PERMUTATIONCOUNT_H
#define PERMUTATIONCOUNT_H

#include "ObservableBase.h"

class PermutationCountClass : public ObservableClass
{
 private:
  int Freq,DumpFreq;
  int TotalCounts;
  int TimesCalled;
  ObservableVecDouble1 CycleCountVar;
  
public:
  Array<double,1> CycleCount;
  int Species;
  void Accumulate();
  void Read(IOSectionClass& in);
  void WriteBlock();
  PermutationCountClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    CycleCountVar("y", IOSection, myPathData.Path.Communicator)
  {
    TimesCalled=0;
  }

};


#endif
