#ifndef OBSERVABLE_ENERGY_H
#define OBSERVABLE_ENERGY_H



#include "ObservableClass.h"




class TotalEnergyClass : public ObservableClass
{

private:
  double ESum;
  double VSum;
  double SSum;
  double FSum;
  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
  VarClass *IOVVar;
  VarClass *IOSVar;
  VarClass *IOUVar;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in) {  
    ObservableClass::Read(in);
    assert(in.ReadVar("freq",Freq));
    assert(in.ReadVar("dumpFreq",DumpFreq));
  }
  TotalEnergyClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) 
  {
    ESum = 0.0;
    VSum = 0.0;
    SSum = 0.0;
    FSum = 0.0;
    NumSamples = 0;
    TimesCalled=0;
  }
};


#endif 
