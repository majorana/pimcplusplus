#ifndef ACTION_BASE_H
#define ACTION_BASE_H

#include "../Common.h"

class PathDataClass;
class PathClass;
class ActionBaseClass
{
protected:
  PathDataClass &PathData;
  PathClass &Path;
public:
  virtual double Evaluate(int slice1, int slice2,
			  const Array<int,1> &activeParticles,
			  int level) = 0;
  ActionBaseClass(PathDataClass &pathData);
					   
};

#endif
