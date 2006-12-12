#include "PlaneObject.h"

void
PlaneObject::SetPosition(int dir, double pos)
{
  Direction = dir;
  Position = pos;
  Set();
}


void
PlaneObject::Init()
{
  MinVal = Spline(0,0,0);
  MaxVal = Spline(0,0,0);
//   double xi = Spline.Xgrid->Start; double xf = Spline.Xgrid->End;
//   double nxInv = 1.0/(double)Spline.Nx;
//   double yi = Spline.Ygrid->Start; double yf = Spline.Ygrid->End;
//   double nyInv = 1.0/(double)Spline.Ny;
//   double zi = Spline.Zgrid->Start; double zf = Spline.Zgrid->End;
//   double nzInv = 1.0/(double)Spline.Nz;
  cerr << "Nx = " << Spline.Nx << endl;
  cerr << "Ny = " << Spline.Ny << endl;
  cerr << "Nz = " << Spline.Nz << endl;
  for (int ix=0; ix<Spline.Nx; ix++) 
    for (int iy=0; iy<Spline.Ny; iy++)
      for (int iz=0; iz<Spline.Nz; iz++) {
	double val = Spline (ix,iy,iz);
	MinVal = (val < MinVal) ? val : MinVal;
	MaxVal = (val > MaxVal) ? val : MaxVal;
      }
  cerr << "MinVal = " << MinVal << endl;
  cerr << "MaxVal = " << MaxVal << endl;
  CMap. Init (MinVal, MaxVal);
  // Allocate a texture
  glGenTextures (1, &TextureNum);
  glBindTexture (GL_TEXTURE_2D, TextureNum);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
  glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
}

PlaneObject::~PlaneObject()
{
  glDeleteTextures(1, &TextureNum);
}

void
PlaneObject::Set()
{
  const int N = 256;
  Array<TinyVector<GLubyte,4>,2> texData(N,N);
  
  Vec3 r0, sVec, tVec;
  Vec3 dr[3];
  dr[0] = (Spline.Xgrid->End - Spline.Xgrid->Start) * Vec3(1.0, 0.0, 0.0);
  dr[1] = (Spline.Ygrid->End - Spline.Ygrid->Start) * Vec3(0.0, 1.0, 0.0);
  dr[2] = (Spline.Zgrid->End - Spline.Zgrid->Start) * Vec3(0.0, 0.0, 1.0);
  r0 = Vec3 (Spline.Xgrid->Start, Spline.Ygrid->Start, Spline.Zgrid->Start);
  r0 += Position*dr[Direction];
  if (Direction == 0) {
    sVec = dr[1];
    tVec = dr[2];
  }
  else if (Direction == 1) {
    sVec = dr[2];
    tVec = dr[0];
  }
  else if (Direction == 2) {
    sVec = dr[0];
    tVec = dr[1];
  }

  // Create the data
  double nInv = 1.0/(double)(N-1);
  for (int is=0; is<N; is++) {
    double s = nInv * (double)is;
    for (int it=0; it<N; it++) {
      double t = nInv * (double)it;
      Vec3 r = r0 + s*sVec + t*tVec;
      double val = Spline(r[0], r[1], r[2]);
      TinyVector<double,4> color;
      CMap (val, color);
      //      cerr << "val = " << val << "  color = " << color << endl;
      texData(is, it)[0] = (GLubyte) floor (256.0 * color[0]);
      texData(is, it)[1] = (GLubyte) floor (256.0 * color[1]);
      texData(is, it)[2] = (GLubyte) floor (256.0 * color[2]);
      texData(is, it)[3] = (GLubyte) floor (256.0 * color[3]);
    }
  }
  
  glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, N, N, 0, GL_RGBA,
	       GL_UNSIGNED_BYTE, texData.data());

  // Create the polygons
  Start();
  glEnable(GL_TEXTURE_2D);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  glBindTexture(GL_TEXTURE_2D, TextureNum);
  glBegin (GL_QUADS);
  Vec3 r1 = r0 + sVec;
  Vec3 r2 = r0 + sVec + tVec;
  Vec3 r3 = r0 + tVec;
  glTexCoord2f(0.0, 0.0); glVertex3f(r0[0], r0[1], r0[2]);
  glTexCoord2f(0.0, 1.0); glVertex3f(r1[0], r1[1], r1[2]);
  glTexCoord2f(1.0, 1.0); glVertex3f(r2[0], r2[1], r2[2]);
  glTexCoord2f(1.0, 0.0); glVertex3f(r3[0], r3[1], r3[2]);
  glEnd();
  glDisable(GL_TEXTURE_2D);
  End();
}
	

void
PlaneObject::DrawPOV(FILE *out, string rotMatrix)
{


}
