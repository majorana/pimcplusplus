#ifndef SPHERE_OBJECT_H
#define SPHERE_OBJECT_H

#include "GLObject.h"
#include "../Common/Blitz.h"

class SphereObject : public GLObject
{
protected:
  Vec3 Color, Pos;
  double Radius;
  void Set();
public:
  void SetPos (Vec3 pos);
  void SetRadius (double radius);
  void SetColor (Vec3 color);
  SphereObject() : Radius (1.0)
  {
    Color = Vec3 (1.0, 0.0, 0.0);
    Pos   = Vec3 (0.0, 0.0, 0.0);
  }
};


#endif
