#include "SphereObject.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

void 
SphereObject::Set()
{
  GLUquadricObj* qobj = gluNewQuadric();
  gluQuadricDrawStyle(qobj, GLU_FILL);
  Start();
  // Set color
  float fcolor[4];
  fcolor[0] = Color[0]; fcolor[1] = Color[1]; fcolor[2] = Color[2];
  fcolor[3] = 1.0;
  glColor3f (Color[0], Color[1], Color[2]);
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
  float spec[4] = { 1.0, 1.0, 1.0, 1.0 };
  glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, spec);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30.0);

  // Create sphere
  glPushMatrix();
  glTranslated (Pos[0], Pos[1], Pos[2]);
  gluSphere(qobj, Radius, 20, 20);
  glPopMatrix();
  End();
}

void 
SphereObject::SetPos (Vec3 pos)
{
  Pos = pos;
  Set();
}

void 
SphereObject::SetRadius (double radius)
{
  Radius = radius;
  Set();
}

void 
SphereObject::SetColor (Vec3 color)
{
  Color = color;
  Set();
}


void
SphereObject::DrawPOV (FILE *fout, string rotString)
{
  fprintf (fout, "sphere {\n");
  fprintf (fout, "  <%10.8f, %10.8f, %10.8f>, %10.8f\n",
	   Pos[0], Pos[1], Pos[2], Radius);
  fprintf (fout, "%s", rotString.c_str());
  fprintf (fout, "  pigment { color rgb <%1.5f %1.5f %1.5f> }\n", 
	   Color[0], Color[1], Color[2]);
  fprintf (fout, "  finish {\n");
  fprintf (fout, "    ambient 0.1\n  diffuse 0.8\n");
  //  fprintf (fout, "    reflection 0.5\n");
  fprintf (fout, "    specular 0.6\n");
  fprintf (fout, "    roughness 0.025 \n");
  fprintf (fout, "  }\n");
  fprintf (fout, "}\n");

}
