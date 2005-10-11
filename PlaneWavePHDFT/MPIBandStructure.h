#ifndef MPI_BAND_STRUCTURE_H
#define MPI_BAND_STRUCTURE_H

#include "PlaneWavesMPI.h"
#include "../MPI/Communication.h"

class MPIBandStructureClass
{
private:
  MPISystemClass *System;
  Array<Vec3,1> Rions;
  Potential *PH;
  double kCut;
  Array<Vec3,1> kPoints;
  int NumBands;
  Vec3 Box;
  string OutFilename;
  int InterpPoints;
  CommunicatorClass Communicator;
public:
  void Read(IOSectionClass &in);
  void CalcBands();
};


#endif
