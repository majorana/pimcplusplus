#ifndef OBSERVABLE_ENERGY_H
#define OBSERVABLE_ENERGY_H



#include "ObservableBase.h"




class EnergyClass : public ObservableClass
{

private:
  double ESum;
  double VSum;
  double SSum;
  double FSum;
  double NodeSum;

  ObservableDouble TotAvg;
  ObservableDouble VAvg;
  ObservableDouble SAvg;
  ObservableDouble FAvg;
  ObservableDouble NAvg;



  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
  VarClass *IOVVar;
  VarClass *IOSVar;
  VarClass *IOUVar;
  VarClass *IONodeVar;
public:
  void Accumulate();
  void WriteBlock();
  void ShiftData(int numTimeSlices);
  void Read(IOSectionClass& in);
  EnergyClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      TotAvg("Total",  IOSection,myPathData.Path.Communicator),
      VAvg("Potential",IOSection,myPathData.Path.Communicator), 
      SAvg("Spring",   IOSection,myPathData.Path.Communicator) , 
      FAvg("dUdBeta",  IOSection,myPathData.Path.Communicator),
      NAvg("Node",     IOSection,myPathData.Path.Communicator)
  {
    ESum = 0.0;
    VSum = 0.0;
    SSum = 0.0;
    FSum = 0.0;
    NodeSum = 0.0;
    NumSamples = 0;
    TimesCalled=0;
  }
};


#endif 
