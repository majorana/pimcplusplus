#ifndef PATHDUMP__H
#define PATHDUMP_H

#include "ObservableBase.h"


class PathDumpClass : public ObservableClass
{
private:
  ObservableVecDouble3 PathVar;
  ObservableVecInt1 PermVar;
  ObservableInt OpenLinkVar;
  ObservableInt OpenLinkPtclVar;
  ObservableInt RefLinkVar;
  ObservableVecDouble1 TailLocVar;
  // Variables for node dumps.
  ObservableVecDouble3 NodeVar;
  bool DumpNodes;
  int  NodePtcl;
  int  NodeSlice;
  LinearGrid Xgrid, Ygrid, Zgrid;
  void NodeDump();
public:
  int TimesCalled;
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& IO);
  int DumpFreq;
  PathDumpClass(PathDataClass &myPathData, IOSectionClass &ioSection) :
    ObservableClass(myPathData, ioSection),
    PathVar ("Path", IOSection, myPathData.Path.Communicator),
    PermVar ("Permutation", IOSection, myPathData.Path.Communicator),
    OpenLinkVar("OpenLinkSlice",IOSection,myPathData.Path.Communicator),
    TailLocVar("TailLocation",IOSection,myPathData.Path.Communicator),
    OpenLinkPtclVar("OpenPtcl",IOSection,myPathData.Path.Communicator),
    RefLinkVar("RefLink",IOSection,myPathData.Path.Communicator),
    NodeVar("Nodes",IOSection,myPathData.Path.Communicator)
  
  { 
    Name="PathDump";
    TimesCalled=0;
  }
};



#endif
