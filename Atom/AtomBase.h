#ifndef ATOM_BASE_H
#define ATOM_BASE_H

#include "RadialWF.h"

typedef enum {DFTType, HFType, OEPType} AtomType;

class Atom
{
protected:
  Grid *grid;
public:
  Array<RadialWF,1> RadialWFs;
  virtual AtomType Type() = 0;
  virtual void UpdateVHXC() = 0;
  virtual void CalcEnergies (double &kinetic, double &potential,
			     double &hartree, double &HXC) = 0;
  virtual void Solve() = 0;
  virtual void Write(IOSectionClass &out) = 0;
  virtual void Read(IOSectionClass &in) = 0;
  virtual void SetGrid(Grid *newGrid) = 0;
  inline Grid* GetGrid(){ return grid; }
  inline double NumElecs();
  virtual void SetBarePot (Potential *pot) = 0;
};

inline double
Atom::NumElecs()
{
  double num = 0.0;
  for (int i=0; i<RadialWFs.size(); i++) 
    num += RadialWFs(i).Occupancy;
  return num;
}


#endif
