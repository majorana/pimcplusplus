#include "BoxObject.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/gle.h>
#include <GL/glut.h>

void 
BoxObject::Set(Vec3 box)
{
  Set (box[0], box[1], box[2]);
}

void 
BoxObject::Set (double lx, double ly, double lz)
{
  Start();
  glPushMatrix();
  glLineWidth (2.0);
  float fcolor[4];
  fcolor[0] = Color[0]; fcolor[1] = Color[1]; fcolor[2] = Color[2];
  fcolor[3] = 1.0;
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
  glColor3d (Color[0], Color[1], Color[2]);
  glScaled (lx, ly, lz);
  glutWireCube (1.0);
  glTranslated (-0.5, -0.5, -0.5);
  glPopMatrix();
  End();
}

void BoxObject::SetColor (double red, double green, double blue)
{
  Color[0] = red;
  Color[1] = green;
  Color[2] = blue;
}
