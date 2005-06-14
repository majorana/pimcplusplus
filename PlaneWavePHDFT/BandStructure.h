#ifndef BAND_STRUCTURE_H
#define BAND_STRUCTURE_H

#include "PlaneWaves.h"

class BandStructureClass
{
private:
  SystemClass *System;
  Array<Vec3,1> Rions;
  Potential *PH;
  double kCut;
  Array<Vec3,1> kPoints;
  int NumBands;
  Vec3 Box;
  string OutFilename;
  int InterpPoints;
public:
  void Read(IOSectionClass &in);
  void CalcBands();

};


#endif
