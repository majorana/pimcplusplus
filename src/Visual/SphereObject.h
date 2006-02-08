#ifndef SPHERE_OBJECT_H
#define SPHERE_OBJECT_H

#include "GLObject.h"
#include <Common/Blitz.h>

class SphereObject : public GLObject
{
protected:
  Vec3 Color, Pos;
  double Radius;
  static int SphereListNum, OffScreenListNum;
  static bool SphereListCreated, OffScreenListCreated;
  bool OffScreen;
  void Set();
public:
  void SetPos (Vec3 pos);
  void SetRadius (double radius);
  void SetColor (Vec3 color);
  void DrawPOV (FILE *fout, string rotString);
  SphereObject(bool offScreen=false);
};


#endif
