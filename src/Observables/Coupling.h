#ifndef COUPLING_H
#define COUPLING_H

#include "ObservableBase.h"

/// A pair correlation function observable.
class CouplingClass : public ObservableClass
{
private:
  ObservableVecDouble1 couplingVar;
  /// Stores number of counts in each bin
  Array<int,1> Coupling;
  /// Stores the total number of counts
  int TotalCounts;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  /// My specialization of the virtual function.
  void Accumulate();
  /// My specialization of the virtual function.
  void Initialize();
  /// My specialization of the virtual function.
  void WriteBlock();
  void WriteInfo();
  void Read(IOSectionClass& IO);
  CouplingClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    couplingVar("y", IOSection, myPathData.Path.Communicator)
  {
    Initialize();
  }


};


#endif
