#ifndef PARTICLE_AVERAGE_LOC_H
#define PARTICLE_AVERAGE_LOC_H

#include "ObservableBase.h"
#include <list>
#include <vector>

class ParticleAverageLocClass : public ObservableClass
{

private:
  int Species;
  Array<double,2> ParticleCenterOfMass;
  ObservableVecDouble2 ParticleAverageLocVar;
  int NumSamples;
  int TimesCalled;
  int Freq;
  int DumpFreq;
public:
  void WriteInfo();
  void Accumulate();
  void WriteBlock();
  void Read(IOSectionClass& in);
  ParticleAverageLocClass(PathDataClass &myPathData, IOSectionClass &ioSection)
    : ObservableClass(myPathData, ioSection) , 
      ParticleAverageLocVar("y",IOSection,myPathData.Path.Communicator)
  {
    NumSamples = 0; 
    TimesCalled=0;
  }
};


#endif 
