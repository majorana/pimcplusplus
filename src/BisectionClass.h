#ifndef BISECTION_CLASS_H
#define BISECTION_CLASS_H

#include "PathDataClass.h"


class BisectionClass 
{

 public:
  PathDataClass  &PathData;
  bool Bisect(int startTimeSlice,int numLevels, Array<int,1> activeParticles);
  double SamplePaths(int startSlice, int endSlice, Array<int,1> particles, int level);
  double LogSampleProb(int startSlice, int endSlice, 
		       Array<int,1> particles, 
		       int level);
  BisectionClass(PathDataClass &pathData) : PathData(pathData)
    { /* Do nothing for now */}
};

#endif
