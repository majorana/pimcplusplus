#ifndef PATHDUMP__H
#define PATHDUMP_H

#include "ObservableBase.h"



class WindingNumberClass : public ObservableClass
{
 protected:
  ObservableVecDouble1 WNVar;
  /// This stores a list of integers corresponding to the species that
  /// are included in the winding number calculation.
  Array<int,1> SpeciesList;
  /// This stores the block sum of the winding numbers
  dVec W2Sum;
  /// Vector which stores all of the local winding numbers since the
  /// last writeblock.
  std::vector<dVec> WNVec;
  Array<double,1> WN2Array;
  int SamplesInBlock;
  void CalcWN2();
 public:
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& IO);
  WindingNumberClass(PathDataClass &myPathData,IOSectionClass &ioSection) :
    ObservableClass(myPathData,ioSection),
    WNVar("W2", IOSection, myPathData.Path.Communicator)
  {
    W2Sum = 0.0;
    SamplesInBlock = 0;
  }

};



#endif
