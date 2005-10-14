#ifndef SUPERFLUIDROP_H
#define SUPERFLUIDROP_H

#include "ObservableBase.h" // what are we inheriting?

class SuperfluiDrop : public ObservableClass
{
 private:
  int Freq,DumpFreq;
  int TotalCounts;
  int TimesCalled;
  //double superfluidity;
  ObservableDouble areaSquared,momInertia;
public:
  double area,anorm;
  double momi,mominorm;
  int Species;
  void Accumulate();
  void Read(IOSectionClass& in);
  void WriteBlock();
  SuperfluiDrop(PathDataClass &myPathData, IOSectionClass &ioSection) : 
    ObservableClass(myPathData,ioSection), 
    areaSquared("area", IOSection, myPathData.Path.Communicator),
    momInertia("mominertia", IOSection, myPathData.Path.Communicator)
  {
    TimesCalled=0;
  }
};


#endif
