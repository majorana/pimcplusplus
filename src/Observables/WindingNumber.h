#ifndef PATHDUMP__H
#define PATHDUMP_H

#include "ObservableBase.h"



class WindingNumberClass : public ObservableClass
{
 private:
  int Freq;
  int DumpFreq;
  dVec TotalW2;
  Array<dVec,1> TotalDisp;
  Array<dVec,1> TempDisp;
 public:
  int TimesCalled;
  int NumSamples;
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& IO);
  WindingNumberClass(PathDataClass &myPathData,IOSectionClass &ioSection):
    ObservableClass(myPathData,ioSection) {
    TimesCalled=0;
    NumSamples=0;
    TotalDisp.resize(PathData.Path.NumParticles());
    TempDisp.resize(PathData.Path.NumParticles());
    TotalW2=0.0;
  }

};



#endif
