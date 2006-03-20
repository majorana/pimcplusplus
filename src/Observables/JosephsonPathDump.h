#ifndef JOSEPHSON_PATH_DUMP_H
#define JOSEPHSON_PATH_DUMP_H

#include "ObservableBase.h"

class JosephsonPathDumpClass : public ObservableClass
{

private:
  Array<double,1> Phase;
  ObservableVecDouble1 PhaseVar;
  
  int NumSamples;
  int TimesCalled;
public:
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& in);
  JosephsonPathDumpClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      PhaseVar("Phase",IOSection,myPathData.Path.Communicator)

  {
    Phase.resize(PathData.Path.NumTimeSlices());
    TimesCalled=0;
    
  }
};


#endif 
