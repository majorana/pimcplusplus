#ifndef VACANCY_LOC_H
#define VACANCY_LOC_H

#include "ObservableBase.h"

class VacancyLocClass : public ObservableClass
{

private:
  //Must be initialized
  Array<int,1> Loc;
  ///This is the set of locations you should compare against to decide
  ///the location of the head and the tail
  Array<dVec,1> FixedLoc;
  ObservableVecDouble1 VacancyLocVar;
  ObservableDouble R2Var;
  double R2Dist;
  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  void WriteInfo();
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& in);
  VacancyLocClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
    VacancyLocVar("y",IOSection,myPathData.Path.Communicator),
    R2Var("R2",IOSection,myPathData.Path.Communicator)
  {
    NumSamples = 0;
    TimesCalled=0;
    R2Dist=0.0;
  }
};


#endif 
