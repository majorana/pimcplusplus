#ifndef PATHDUMP__H
#define PATHDUMP_H

#include "ObservableBase.h"


class PathDumpClass : public ObservableClass
{
private:
  ObservableVecDouble3 PathVar;
  ObservableVecInt1 PermVar;
public:
  int TimesCalled;
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& IO);
  int DumpFreq;
  PathDumpClass(PathDataClass &myPathData, IOSectionClass &ioSection) :
    ObservableClass(myPathData, ioSection),
    PathVar ("Path", IOSection, myPathData.Path.Communicator),
    PermVar ("Permutation", IOSection, myPathData.Path.Communicator)
  { 
    Name="PathDump";
    TimesCalled=0;
  }
};



#endif
