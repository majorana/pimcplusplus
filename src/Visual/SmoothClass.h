#ifndef SMOOTH_CLASS_H
#define SMOOTH_CLASS_H

#include "OnePath.h"

class SmoothClass
{
private:
  int PointsPerPath;
  double SmoothLevel;
public:
  // Set smoothness level from 0 to 1
  void SetLevel (double level);

  // Perform fourier smoothing of the paths
  void SmoothPaths (vector<OnePath*> &inList);
  SmoothClass() : PointsPerPath(300), SmoothLevel (0.5)
  {

  }
};


#endif
