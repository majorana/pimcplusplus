#ifndef ACTION_BASE_H
#define ACTION_BASE_H

#include "../PathDataClass.h"

class ActionBase
{
protected:
  PathDataClass &PathData;
public:
  virtual double Evaluate(int slice1, int slice2,
			  const Array<int,1> &activeParticles,
			  int level) = 0;
  ActionBase(PathDataClass &pathData) : PathData(pathData)
  {
    /* Do nothing */
  }
};

#endif
