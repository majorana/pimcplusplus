#include "BoxObject.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/gle.h>
#include <GL/glut.h>

void 
BoxObject::Set(Vec3 box, bool useClip)
{
  Set (box[0], box[1], box[2], useClip);
}

void 
BoxObject::Set (double lx, double ly, double lz, bool useClip)
{

  Lx=lx; Ly=ly; Lz=lz;
  Start();
  if (useClip) {
    GLdouble eqn0[4] = {  1.0,  0.0,  0.0,  0.5001*lx };
    GLdouble eqn1[4] = { -1.0,  0.0,  0.0,  0.5001*lx };
    GLdouble eqn2[4] = {  0.0,  1.0,  0.0,  0.5001*ly };
    GLdouble eqn3[4] = {  0.0, -1.0,  0.0,  0.5001*ly };
    GLdouble eqn4[4] = {  0.0,  0.0,  1.0,  0.5001*lz };
    GLdouble eqn5[4] = {  0.0,  0.0, -1.0,  0.5001*lz };
    glClipPlane(GL_CLIP_PLANE0, eqn0);
    glClipPlane(GL_CLIP_PLANE1, eqn1);
    glClipPlane(GL_CLIP_PLANE2, eqn2);
    glClipPlane(GL_CLIP_PLANE3, eqn3);
    glClipPlane(GL_CLIP_PLANE4, eqn4);
    glClipPlane(GL_CLIP_PLANE5, eqn5);
    glEnable(GL_CLIP_PLANE0);
    glEnable(GL_CLIP_PLANE1);
    glEnable(GL_CLIP_PLANE2);
    glEnable(GL_CLIP_PLANE3);
    glEnable(GL_CLIP_PLANE4);
    glEnable(GL_CLIP_PLANE5);
  }
  else {
    glDisable(GL_CLIP_PLANE0);
    glDisable(GL_CLIP_PLANE1);
    glDisable(GL_CLIP_PLANE2);
    glDisable(GL_CLIP_PLANE3);
    glDisable(GL_CLIP_PLANE4);
    glDisable(GL_CLIP_PLANE5);
  }
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


void 
BoxObject::POVLine (FILE *fout, 
		    double x1, double y1, double z1,
		    double x2, double y2, double z2,
		    double radius, string rotString)
{
  fprintf (fout, "cylinder {\n");
  fprintf (fout, "  <%12.8f, %12.8f, %12.8f>,\n",
	   x1, y1, z1);
  fprintf (fout, "  <%12.8f, %12.8f, %12.8f>,\n",
	   x2, y2, z2);
  fprintf (fout, "  %10.8f\n", radius);
  fprintf (fout, "  pigment { color rgb <%1.5f %1.5f %1.5f> }\n", 
	   0.0, 0.0, 0.0);
  fprintf (fout, "%s", rotString.c_str());
  fprintf (fout, "}\n\n");
}

void
BoxObject::DrawPOV (FILE *fout, string rotString)
{
  double minDim = min(min(Lx,Ly),Lz);
  double radius = minDim/450.0;
  POVLine (fout, -0.5*Lx, -0.5*Ly, -0.5*Lz,
	          0.5*Lx ,-0.5*Ly, -0.5*Lz, radius, rotString);
  POVLine (fout, -0.5*Lx,  0.5*Ly, -0.5*Lz,
	          0.5*Lx , 0.5*Ly, -0.5*Lz, radius, rotString);
  POVLine (fout, -0.5*Lx, -0.5*Ly, -0.5*Lz,
	         -0.5*Lx , 0.5*Ly, -0.5*Lz, radius, rotString);
  POVLine (fout,  0.5*Lx, -0.5*Ly, -0.5*Lz,
	          0.5*Lx , 0.5*Ly, -0.5*Lz, radius, rotString);

  POVLine (fout, -0.5*Lx, -0.5*Ly,  0.5*Lz,
	          0.5*Lx ,-0.5*Ly,  0.5*Lz, radius, rotString);
  POVLine (fout, -0.5*Lx,  0.5*Ly,  0.5*Lz,
	          0.5*Lx , 0.5*Ly,  0.5*Lz, radius, rotString);
  POVLine (fout, -0.5*Lx, -0.5*Ly,  0.5*Lz,
	         -0.5*Lx , 0.5*Ly,  0.5*Lz, radius, rotString);
  POVLine (fout,  0.5*Lx, -0.5*Ly,  0.5*Lz,
	          0.5*Lx , 0.5*Ly,  0.5*Lz, radius, rotString);

  POVLine (fout, -0.5*Lx, -0.5*Ly, -0.5*Lz,
	         -0.5*Lx ,-0.5*Ly,  0.5*Lz, radius, rotString);
  POVLine (fout, -0.5*Lx,  0.5*Ly, -0.5*Lz,
	         -0.5*Lx , 0.5*Ly,  0.5*Lz, radius, rotString);
  POVLine (fout,  0.5*Lx, -0.5*Ly, -0.5*Lz,
	          0.5*Lx ,-0.5*Ly,  0.5*Lz, radius, rotString);
  POVLine (fout,  0.5*Lx,  0.5*Ly, -0.5*Lz,
	          0.5*Lx , 0.5*Ly,  0.5*Lz, radius, rotString);

}
