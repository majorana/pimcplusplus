#ifndef HBOND_H
#define HBOND_H

#include "ObservableBase.h"

class HbondClass : public ObservableClass
{

private:
  double TotalSum;

  ObservableDouble TotalVar;

  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  bool CheckPair(int slice, int obond, int ohome, int p);
  bool IsHBond(int slice, int OA, int OB);
  void Accumulate();
  void WriteBlock();
  void Initialize();
  void Read(IOSectionClass& in);
  HbondClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      TotalVar  ("Total",  IOSection,myPathData.Path.Communicator)
  {
    Initialize();
  }
};
#endif 
