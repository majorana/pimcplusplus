#ifndef PLANE_DENSITY_H
#define PLANE_DENSITY_H

#include "ObservableBase.h"
#include <list>
#include <vector>

class PlaneDensityClass : public ObservableClass
{

private:
  Array<double,2> Grid;
  ObservableVecDouble2 GridVar;
  int IntoGrid(double num);
  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  void WriteInfo();
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& in);
  PlaneDensityClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      GridVar("y",IOSection,myPathData.Path.Communicator)
  {
    NumSamples = 0; 
    TimesCalled=0;
    int numPoints=50;
    Grid.resize(numPoints,numPoints);
    Grid=0.0;
  }
};


#endif 
