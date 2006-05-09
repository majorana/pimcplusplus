#ifndef VARIATIONALPI_ENERGY_H
#define VARIATIONALPI_ENERGY_H



#include "ObservableBase.h"




class VariationalPIEnergyClass : public ObservableClass
{

private:
  double Energy;
  Array<double,2> DetMatrix;
  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
  ObservableDouble EnergyVar;
public:
  double rho(int i,int j);
  double DRho(int i, int j);
  int NumImages;
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& in);
  VariationalPIEnergyClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection),
      EnergyVar("Energy",IOSection,PathData.Path.Communicator)
  {
    Energy=0.0;
    NumSamples = 0;
    TimesCalled=0;
    NumImages=0;
  }
};


#endif 
