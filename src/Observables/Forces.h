#ifndef FORCES_H
#define FORCES_H

#include "ObservableBase.h"

/// A pair correlation function observable.
class ForcesClass : public ObservableClass
{
private:
  ObservableVecDouble2 ForcesVar;
  /// Stores number of counts in each bin
  Array<dVec,1> Forces;
  Array<double,2> ForcesArray;
  int SpeciesNum;
  Array<int,1> Ptcls;
  /// Stores the total number of counts
  int TimesCalled, Counts;
  int Freq;
  int DumpFreq;
public:
  /// My specialization of the virtual function.
  void Accumulate();
  /// My specialization of the virtual function.
  void Initialize();
  /// My specialization of the virtual function.
  void Print();
  void WriteBlock();
  void WriteInfo();
  void Read(IOSectionClass& IO);
  void SetSpecies (int speciesNum);
  ForcesClass(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), Counts(0),
    ForcesVar("F", IOSection, myPathData.Path.Communicator)
  {
    TimesCalled=0;
  }
  ForcesClass(PathDataClass &myPathData, IOSectionClass &ioSection,
	      int speciesNum) : 
    ObservableClass(myPathData, ioSection), Counts(0),
    ForcesVar("y", IOSection, myPathData.Path.Communicator)
  { 
    SetSpecies (speciesNum);
  }
};

#endif
