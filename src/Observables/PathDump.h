#ifndef PATHDUMP__H
#define PATHDUMP_H

#include "ObservableBase.h"


class PathDumpClass : public ObservableClass
{
private:
public:
  int TimesCalled;
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& IO);
  PathDumpClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection)  { 
    Name="PathDump";
    TimesCalled=0;
  }
};



#endif
