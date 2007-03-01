#ifndef PLANE_OBJECT_H
#define PLANE_OBJECT_H

#include "GLObject.h"
#include <GL/glut.h>
#include <Common/Blitz.h>
#include "ColorMap.h"
#include <Common/Splines/MyTricubicSpline.h>

class PlaneObject : public GLObject
{
private:
  MyTricubicSpline &Spline;
  int Direction;
  double Position;
  void Set();
  ColorMap CMap;
  ColorMapType MapType;
  double MinVal, MaxVal;
  GLuint TextureNum;
  bool HaveTexture, BuiltTexture;
  bool IsInitialized;
  Mat3 Lattice;
public:
  void Init();
  void DrawPOV (FILE *out, string rotMatrix);
  void SetPosition (int dir, double pos);
  void SetLattice(Mat3 lattice);
  void SetColorMap (ColorMapType map);
  PlaneObject (MyTricubicSpline &spline) : 
    Spline (spline), Direction(0), Position(0),
    HaveTexture(false), BuiltTexture(false),
    IsInitialized(false), MapType (BLUE_WHITE_RED)
  {
    Lattice = 1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 1.0;
  }
  ~PlaneObject();
};

#endif
