#ifndef STRUCTURE_FACTOR_H
#define STRUCTURE_FACTOR_H

#include "ObservableBase.h"

/// A pair correlation function observable.
class StructureFactorClass : public ObservableClass
{
  /// Stores the total number of counts

  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  int TotalCounts;
  /// Stores number of counts in each bin
  Array<double,1> Sk;
  /// The species between which I am calculating the pair correlation
  /// function.
  int Species1, Species2;
  /// My specialization of the virtual function.
  void Accumulate();
  /// My specialization of the virtual function.
  void Initialize();
  /// My specialization of the virtual function.
  void WriteBlock();
  void WriteInfo();
  void Read(IOSectionClass& IO);
  void Calculate();
  void Clear();
  StructureFactorClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection){
    TimesCalled=0;
  }
  StructureFactorClass(PathDataClass &myPathData, IOSectionClass &ioSection,
		       int species1, int species2) : 
    ObservableClass(myPathData, ioSection), 
    Species1(species1), Species2(species2) 
  { Initialize(); }
};



#endif
