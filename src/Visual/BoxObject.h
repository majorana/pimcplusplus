#ifndef BOX_OBJECT_H
#define BOX_OBJECT_H

#include "GLObject.h"
#include "../Common/Blitz.h"

class BoxObject : public GLObject
{
  Vec3 Color;
public:
  void Set (Vec3 box);
  void Set (double lx, double ly, double lz);
  void SetColor (double red, double blue, double green);
  BoxObject()
  {
    Color = Vec3(0.2, 0.2, 0.2);
  }
};

#endif
