#ifndef ACTION_BASE_H
#define ACTION_BASE_H

#include "../Common.h"

class PathDataClass;

class ActionBaseClass
{
protected:
  PathDataClass &PathData;
public:
  virtual double Evaluate(int slice1, int slice2,
			  const Array<int,1> &activeParticles,
			  int level) = 0;
  ActionBaseClass(PathDataClass &pathData) : PathData(pathData)
  {
    /* Do nothing */
  }
};

#endif
