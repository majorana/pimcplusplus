#ifndef BISECTION_CLASS_H
#define BISECTION_CLASS_H

#include "../PathDataClass.h"


class BisectionClass 
{

 public:
  PathDataClass  &PathData;
  bool Bisect(int startTimeSlice,int numLevels, Array<int,1> activeParticles);
  bool Bisect(int startTimeSlice,int numLevels, Array<int,1> activeParticles,
	      double permActionChange);

  ///This picks a new location in space for the particles in the
  ///particles Array at all of the time slices between startSlice and
  ///endSlice (at the appropriate skip for the level)
  double SamplePaths(int startSlice, int endSlice, Array<int,1> particles, int level);

  /// This calculates the sample probability for going from the state
  /// that is currently in the newMode of MirroredArrayClass to the
  /// state that is currently in oldMode of MirroredArrayClass 
  double LogSampleProb(int startSlice, int endSlice, 
		       Array<int,1> particles, 
		       int level);
  BisectionClass(PathDataClass &pathData) : PathData(pathData)
    { /* Do nothing for now */}
};

#endif
