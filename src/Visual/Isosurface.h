#ifndef ISOSURFACE_H
#define ISOSURFACE_H

#include <Common/Splines/MyTricubicSpline.h>
#include <GLObject.h>
#include <vector>

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
  /// Returns the number of real solutions
  int CubicFormula (double a, double b, double c, double d,
		     double &x1, double &x2, double &x3);
  Vec3 FindEdge(int ix, int iy, int iz, int edgeNum);
  bool UseCubicInterp;
  bool UseNormals;
  vector<double> Isovals;
  double Isoval;
  TinyVector<double,4> Color;
  vector<TinyVector<double,4> > Colors;
  void Set();
public:
  int NumTriangles();
  inline void SetIsoval (double val) 
  { Isovals.resize(1);  Isovals[0] = val; Set();}
  inline void SetIsoval (vector<double> vals) 
  { Isovals.resize(vals.size()); Isovals = vals; Set(); }
  void DrawPOV(FILE* out, string rotMatrix);
  inline void SetColor (Vec3 color) 
  { Color[0]=color[0]; Color[1]=color[1]; Color[2]=color[2]; }
  inline void SetAlpha (double alpha) { Color[3] = alpha; }
  inline void SetColor (double r, double g, double b) 
  { SetColor (TinyVector<double,3>(r,g,b)); }
  Isosurface() : UseCubicInterp(true), UseNormals(true)
  {
    Color[0] = 0.0; Color[1] = 0.8; Color[2] = 0.0; Color[3] = 0.5;
  }
};
#endif
