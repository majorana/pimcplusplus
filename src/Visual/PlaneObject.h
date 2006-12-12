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
  double MinVal, MaxVal;
  GLuint TextureNum;
  bool HaveTexture, BuiltTexture;
public:
  void Init();
  void DrawPOV (FILE *out, string rotMatrix);
  void SetPosition (int dir, double pos);
  PlaneObject (MyTricubicSpline &spline) : 
    Spline (spline), Direction(0), Position(0),
    HaveTexture(false), BuiltTexture(false)
  {
  }
  ~PlaneObject();
};

#endif
