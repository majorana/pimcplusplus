#ifndef PATH_OBJECT_H
#define PATH_OBJECT_H

#include "GLObject.h"
#include "../Common/Blitz.h"


class PathObject : public GLObject
{
protected:
  Vec3 Color;
  void Cylinder(const Vec3 &r1, const Vec3 &r2);
  bool Closed;
public:
  void Set (Array<Vec3, 1> &path);
  void SetColor (double red, double green, double blue);
  PathObject() : Closed(true)
  {
    Color = Vec3(0.0, 0.0, 1.0);
  }
};



#endif

