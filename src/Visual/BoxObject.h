#ifndef BOX_OBJECT_H
#define BOX_OBJECT_H

#include "GLObject.h"
#include "../Common/Blitz.h"

class BoxObject : public GLObject
{
  Vec3 Color;
  double Lx, Ly, Lz;
  void POVLine (FILE *fout, 
		double x1, double y1, double z1,
		double x2, double y2, double z2, double radius);
public:
  void Set (Vec3 box);
  void Set (double lx, double ly, double lz);
  void SetColor (double red, double blue, double green);
  void DrawPOV (FILE *fout);
  BoxObject()
  {
    Color = Vec3(0.2, 0.2, 0.2);
  }
};

#endif
