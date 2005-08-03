#ifndef ISOSURFACE_H
#define ISOSURFACE_H

#include <Common/Splines/MyTricubicSpline.h>
#include <GLObject.h>

/// This class uses to Marching Cubes algorithm to construct an
/// isosurface from data tabulated on a 3D grid.  It inherits from
/// both MyTricubicSpline and GLObject.  You simply set up the
/// tricubic spline as usual and then call SetIsoval to compute set up
/// the isosurface with the Marching Cubes algorithm.
class Isosurface : public MyTricubicSpline, public GLObject
{
private:
  static int EdgeData[256][13];
  static int EdgeTable[12][7];
  Vec3 FindEdge(int ix, int iy, int iz, int edgeNum);
  bool UseCubicInterp;
  bool UseNormals;
  double Isoval;
  Vec3 Color;
  void Set();
public:
  inline void SetIsoval (double val) { Isoval = val; Set(); }
  void DrawPOV(FILE* out, string rotMatrix);
  Isosurface() : UseCubicInterp(false), UseNormals(true)
  {
    Color[0] = 0.0; Color[1] = 0.8; Color[2] = 0.0;
  }
};
#endif
