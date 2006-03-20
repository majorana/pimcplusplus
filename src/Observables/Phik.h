#ifndef PHI_K_H
#define PHI_K_H

#include "ObservableBase.h"

class PhiKClass : public ObservableClass
{

private:
  Array<double,1> PhiK;
  ObservableVecDouble1 PhiKVar;
  
  

  int NumSamples;
  int TimesCalled;
public:
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& in);
  PhiKClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      PhiKVar("PhiK",IOSection,myPathData.Path.Communicator)

  {
    PhiK.resize(PathData.Path.NumTimeSlices());
    PhiK=0;
    TimesCalled=0;
    
  }
};


#endif 
