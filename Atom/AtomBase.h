#ifndef ATOM_BASE_H
#define ATOM_BASE_H

#include "RadialWF.h"

class Atom
{
protected:
  Grid *grid;
public:
  Array<RadialWF,1> RadialWFs;
  virtual void UpdateVHXC() = 0;
  virtual void CalcEnergies (double &kinetic, double &potential,
			     double &hartree, double &HXC) = 0;
  virtual void Solve() = 0;
  virtual void Write(IOSectionClass &out) = 0;
  virtual void Read(IOSectionClass &in) = 0;
  virtual void SetGrid(Grid *newGrid) = 0;
  virtual void SetBarePot (Potential *pot) = 0;
};



#endif
