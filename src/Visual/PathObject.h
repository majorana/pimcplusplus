#ifndef PATH_OBJECT_H
#define PATH_OBJECT_H

#include "GLObject.h"
#include "../Common/Blitz.h"
#include <vector>

using namespace std;

class PathObject : public GLObject
{
protected:
  Vec3 Color;
  void Cylinder(const Vec3 &r1, const Vec3 &r2);
  bool Closed;
  double Radius;
public:
  void LinesSet (vector<Vec3> &path);
  void TubesSet (vector<Vec3> &path);
  void SetColor (double red, double green, double blue);
  inline void SetRadius (double radius) 
  { Radius = radius; }
  PathObject() : Closed(true), Radius (0.1)
  {
    Color = Vec3(0.0, 0.0, 1.0);
  }
};



#endif

