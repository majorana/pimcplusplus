#include "PathObject.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/gle.h>

void PathObject::Set(Array<Vec3,1> &path)
{
  Start();
  glColor3d (Color[0], Color[1], Color[2]);
  float fcolor[4];
  fcolor[0] = Color[0]; fcolor[1] = Color[1]; fcolor[2] = Color[2];
  fcolor[3] = 1.0;
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
  float spec[4] = { 1.0, 1.0, 1.0, 1.0 };
  glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, spec);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30.0);

  gleSetJoinStyle (TUBE_JN_ROUND /* | TUBE_JN_CAP */ | TUBE_CONTOUR_CLOSED);
//   glBegin(GL_LINE_STRIP);
//   for (int i=0; i<path.size(); i++)
//     glVertex3dv ((double*)&path(i));
//   glEnd();
  
  int N = path.size();
  gleDouble pointArray[N+3][3];
  if (Closed) {
    for (int i=0; i<3; i++) {
      pointArray[0][i] = path(N-2)[i];
      pointArray[N+1][i] = path(1)[i];
      pointArray[N+2][i] = path(2)[i];
    }
  }
  else {  // Fix this
    for (int i=0; i<3; i++) {
      pointArray[0][i] = path(N-2)[i];
      pointArray[N+1][i] = path(1)[i];
      pointArray[N+2][i] = path(2)[i];
    }
  }

  for (int i=0; i<N; i++)
    for (int j=0; j<3; j++)
      pointArray[i+1][j] = path(i)[j];

  float colors[N+3][3];
  for (int i=0; i<N+3; i++) {
    colors[i][0] = Color[0];    
    colors[i][1] = Color[1];
    colors[i][2] = Color[2];
  }
  glePolyCylinder (N+3, pointArray,
		   colors, Radius);

  End();
}

void PathObject::SetColor(double red, double green, double blue)
{
  Color = Vec3 (red, green, blue);
}


void PathObject::Cylinder(const Vec3 &r1, const Vec3 &r2)
{
  Vec3 delta = r2-r1;
  Vec3 deltaNorm = (1.0/sqrt(dot(delta,delta)))*delta;

  Vec3 unit(1.0, 0.0, 0.0);
  if (dot(unit, deltaNorm) > 0.99) 
    unit = Vec3(0.0, 1.0, 0.0);

  Vec3 normal = cross (deltaNorm, unit);
  

}
