#ifndef FORCES_H
#define FORCES_H

#include "../Common.h"

class PathDataClass;

class ForcesClass
{
protected:
  PathDataClass &PathData;
  /// This specifies which species I'm computing the forces on;
  int SpeciesNum;
  /// This is the variable in which we accumlate forces.
  Array<dVec,1> Fsum;
  int Counts;
  /// This stores the current forces experienced by the ions at each
  /// link.  Indexed by (link,ptcl)
  /// Array<dVec,2> CurrentF;
public:
  void ShiftData (int slices);
  void SetSpecies (int speciesNum);
  /// This computes the current force at all time steps
  void Accumulate();
  /// Resets the internal average force variable to zero.
  void Reset();
  void GetForces(Array<dVec,1> &F);
};

#endif
