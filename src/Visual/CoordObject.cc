#include "CoordObject.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/glut.h>

void
CoordObject::Set(Vec3 box)
{
  Set (box[0], box[1], box[2]);
}

void
CoordObject::Set (double lx, double ly, double lz)
{
  Lx=1.1*lx; Ly=1.1*ly; Lz=1.1*lz;
  Start();
  
  // Make sure the coordinates don't get clipped
  glDisable(GL_CLIP_PLANE0);
  glDisable(GL_CLIP_PLANE1);
  glDisable(GL_CLIP_PLANE2);
  glDisable(GL_CLIP_PLANE3);
  glDisable(GL_CLIP_PLANE4);
  glDisable(GL_CLIP_PLANE5);

  glPushMatrix();
  
  glTranslated (-0.5*Lx, -0.5*Ly, -0.5*Lz);
  float red[4]   = {1.0, 0.0, 0.0, 1.0 };
  float green[4] = {0.0, 1.0, 0.0, 1.0 };
  float blue[4]  = {0.0, 0.0, 1.0, 1.0 };
  
 
  glPushMatrix();
  glRotated (90.0, 0.0, 1.0, 0.0);
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, red);
  glColor3d (red[0], red[1], red[2]);
  GLUquadricObj* qobj = gluNewQuadric();
  gluCylinder(qobj, 0.015*Lx, 0.015*Lx, 0.25*Lx, 20, 1);
  gluDeleteQuadric(qobj);
  glTranslated(0.0, 0.0, 0.25*Lx);
  glutSolidCone(0.025*Lx, 0.04*Lx, 20, 1);
  glTranslated(0.0, -0.02*Lx, 0.06*Lx);
  glRotated(-90.0, 0.0, 1.0, 0.0);
  glScaled (0.06*Lx/119.05, 0.06*Lx/119.05, 0.06*Lx/119.05);
  glutStrokeCharacter (GLUT_STROKE_ROMAN, (int)'x');
  glPopMatrix();

  glPushMatrix();
  glRotated (-90.0, 1.0, 0.0, 0.0);
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, green);
  glColor3d (green[0], green[1], green[2]);
  qobj = gluNewQuadric();
  gluCylinder(qobj, 0.015*Lx, 0.015*Lx, 0.25*Lx, 20, 1);
  gluDeleteQuadric(qobj);
  glTranslated(0.0, 0.0, 0.25*Lx);
  glutSolidCone(0.025*Lx, 0.04*Lx, 20, 1);
  glTranslated(0.0, -0.02*Lx, 0.06*Lx);
  glRotated(-90.0, 0.0, 1.0, 0.0);
  glScaled (0.06*Lx/119.05, 0.06*Lx/119.05, 0.06*Lx/119.05);
  glutStrokeCharacter (GLUT_STROKE_ROMAN, (int)'y');
  glPopMatrix();

  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, blue);
  qobj = gluNewQuadric();
  glColor3d (blue[0], blue[1], blue[2]);
  gluCylinder(qobj, 0.015*Lx, 0.015*Lx, 0.25*Lx, 20, 1);
  glTranslated(0.0, 0.0, 0.25*Lx);
  glutSolidCone(0.025*Lx, 0.04*Lx, 20, 1);
  glTranslated(0.0, -0.02*Lx, 0.06*Lx);
  glRotated(-90.0, 0.0, 1.0, 0.0);
  glScaled (0.06*Lx/119.05, 0.06*Lx/119.05, 0.06*Lx/119.05);
  glutStrokeCharacter (GLUT_STROKE_ROMAN, (int)'z');
  gluDeleteQuadric(qobj);

  
  glPopMatrix();
  
  End();
}


void
CoordObject::DrawPOV(FILE *fout, string rotString)
{



}
