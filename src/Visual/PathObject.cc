#include "PathObject.h"
#include <GL/gl.h>
#include <GL/glu.h>
#include <GL/gle.h>

void PathObject::LinesSet(vector<Vec3> &path)
{
  Start();
  glColor3d (Color[0], Color[1], Color[2]);
  float fcolor[4];
  fcolor[0] = Color[0]; fcolor[1] = Color[1]; fcolor[2] = Color[2];
  fcolor[3] = 1.0;
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
  glLineWidth(3.0);
  glBegin(GL_LINE_STRIP);
  for (int i=0; i<path.size(); i++) {
//     fcolor[0] = Color[0]*(double)i/(double)(path.size()-1);
//     fcolor[1] = Color[1]*(double)i/(double)(path.size()-1);
//     fcolor[2] = Color[2]*(double)i/(double)(path.size()-1);
//     glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
    glVertex3dv ((double*)&path[i]);
  }
  glEnd();
  End();
}

void PathObject::TubesSet(vector<Vec3> &path)
{
  Start();
  double alpha;
  int N = path.size();
  Vec3 centroid(0.0, 0.0, 0.0);

  for (int i=0; i<N; i++) 
    centroid += path[i];  

  centroid *= 1.0/N;

  float fcolor[4];
  if (Closed && (fabs(centroid[0]) > 2.0)) {
    alpha = 0.2;
    glBlendFunc(GL_SRC_ALPHA, GL_ONE);
    glBlendEquation(GL_FUNC_REVERSE_SUBTRACT);
    fcolor[0]=1.0-Color[0]; fcolor[1]=1.0-Color[1]; fcolor[2]=1.0-Color[2]; 
  }
  else {
    alpha = 1.0;
    glBlendFunc(GL_SRC_ALPHA, GL_ONE_MINUS_SRC_ALPHA);
    glBlendEquation(GL_FUNC_ADD);
    fcolor[0]=Color[0]; fcolor[1]=Color[1]; fcolor[2]=Color[2];
    //    glBlendFunc(GL_ONE, GL_ZERO);
  }

  glColor3d (Color[0], Color[1], Color[2]);

  fcolor[3] = alpha;
  glMaterialfv (GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE, fcolor);
  float spec[4] = { 1.0, 1.0, 1.0, alpha};
  glMaterialfv (GL_FRONT_AND_BACK, GL_SPECULAR, spec);
  glMaterialf(GL_FRONT_AND_BACK, GL_SHININESS, 30.0);

  glColorMaterial(GL_FRONT_AND_BACK, GL_AMBIENT_AND_DIFFUSE);

  if (Closed)
    gleSetJoinStyle (TUBE_JN_ROUND /* | TUBE_JN_CAP */ | TUBE_CONTOUR_CLOSED);
  else
    gleSetJoinStyle (TUBE_JN_ROUND  | TUBE_JN_CAP);


  gleDouble pointArray[N+3][3];
  if (Closed) {
    for (int i=0; i<3; i++) {
      pointArray[0][i] = path[N-2][i];
      pointArray[N+1][i] = path[1][i];
      pointArray[N+2][i] = path[2][i];
    }
  }
  else {  // Fix this
    for (int i=0; i<3; i++) {
      pointArray[0][i] = path[0][i];//2*path[0][i]-path[1][i];
      pointArray[N+1][i] = path[N-1][i];
      pointArray[N+2][i] = path[N-1][i];//2*path[N][i]-path[N-1][i];
    }
  }


  for (int i=0; i<N; i++) {
    for (int j=0; j<3; j++)
      pointArray[i+1][j] = path[i][j];
  }
  centroid *= 1.0/N;
  
  float colors[N+3][4];
  //  float colors[N+3][e];
  for (int i=0; i<N+3; i++) {
    double a = (double)i/(double)(path.size()-1);
    double b = 1.0-a;
    colors[i][0] = fcolor[0];//*a + (1.0-Color[0])*b;
    colors[i][1] = fcolor[1];//*a + (1.0-Color[1])*b;
    colors[i][2] = fcolor[2];//*a + (1.0-Color[2])*b;
    colors[i][3] = alpha;
//     colors[i][1] = Color[1]*(double)i/(double)(path.size()-1);
//     colors[i][2] = Color[2]*(double)i/(double)(path.size()-1);
  }
  //  glePolyCylinder (N+3, pointArray, colors, Radius);
  glePolyCylinder_c4f (N+3, pointArray, colors, Radius);

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
