#ifndef ACTION_BASE_H
#define ACTION_BASE_H

#include "../Common.h"
#include "../Common/IO/InputOutput.h"

class PathDataClass;
class PathClass;
class ActionBaseClass
{
protected:
  PathDataClass &PathData;
  PathClass &Path;
public:
  virtual double Action(int slice1, int slice2,
			const Array<int,1> &activeParticles,
			int level) = 0;
  virtual double d_dBeta (int slice1, int slice2,
			  int level) = 0;
  // This returns the sum over all time slices, using MPI
  // to sum over processors if necessary.
  virtual void Read (IOSectionClass &in);
  ActionBaseClass(PathDataClass &pathData);				   
};

class PotentialBaseClass
{
protected:
  PathDataClass &PathData;
  PathClass &Path;
public:
  virtual double V (int slice) = 0;
  PotentialBaseClass (PathDataClass &pathData);
};

#endif
