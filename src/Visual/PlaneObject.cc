#include "PlaneObject.h"
#include <sstream>
#include <gtkmm.h>

void
PlaneObject::SetPosition(int dir, double pos)
{
  Direction = dir;
  Position = pos;
  Set();
}


void
PlaneObject::SetColorMap(ColorMapType map)
{
  MapType = map;
  CMap.Init (MinVal, MaxVal, map);
}

void
PlaneObject::Init()
{
  MinVal = Spline(0,0,0);
  MaxVal = Spline(0,0,0);
  for (int ix=0; ix<Spline.Nx; ix++) 
    for (int iy=0; iy<Spline.Ny; iy++)
      for (int iz=0; iz<Spline.Nz; iz++) {
	double val = Spline (ix,iy,iz);
	MinVal = (val < MinVal) ? val : MinVal;
	MaxVal = (val > MaxVal) ? val : MaxVal;
      }
  CMap. Init (MinVal, MaxVal, MapType);
  IsInitialized = true;
}

PlaneObject::~PlaneObject()
{
  if (HaveTexture) 
    glDeleteTextures(1, &TextureNum);
}

void
PlaneObject::Set()
{
  if (!IsInitialized)
    Init();
  if (!HaveTexture) {
    // Allocate a texture
    glGenTextures (1, &TextureNum);
    glBindTexture (GL_TEXTURE_2D, TextureNum);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_S, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_WRAP_T, GL_REPEAT);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri(GL_TEXTURE_2D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    HaveTexture = true;
  }


  const int N = 256;
  Array<TinyVector<GLubyte,4>,2> texData(N,N);
  
  Vec3 u0, sVec, tVec;
  Vec3 dr[3];
  dr[0] = (Spline.Xgrid->End - Spline.Xgrid->Start) * Vec3(1.0, 0.0, 0.0);
  dr[1] = (Spline.Ygrid->End - Spline.Ygrid->Start) * Vec3(0.0, 1.0, 0.0);
  dr[2] = (Spline.Zgrid->End - Spline.Zgrid->Start) * Vec3(0.0, 0.0, 1.0);
  u0 = Vec3 (Spline.Xgrid->Start, Spline.Ygrid->Start, Spline.Zgrid->Start);
  u0 += Position*dr[Direction];
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
      Vec3 r = u0 + s*sVec + t*tVec;
      double val = Spline(r[0], r[1], r[2]);
      TinyVector<double,4> color;
      CMap (val, color);
      texData(is, it)[0] = (GLubyte) floor (256.0 * color[0]);
      texData(is, it)[1] = (GLubyte) floor (256.0 * color[1]);
      texData(is, it)[2] = (GLubyte) floor (256.0 * color[2]);
      texData(is, it)[3] = (GLubyte) floor (256.0 * color[3]);
    }
  }
  
  glBindTexture(GL_TEXTURE_2D, TextureNum);
  if (!BuiltTexture) {
    glTexImage2D(GL_TEXTURE_2D, 0, GL_RGBA, N, N, 0, GL_RGBA,
		 GL_UNSIGNED_BYTE, texData.data());
    BuiltTexture = true;
  }
  else
    glTexSubImage2D(GL_TEXTURE_2D, 0, 0, 0, N, N, 
		    GL_RGBA, GL_UNSIGNED_BYTE, texData.data());
  

  // Create the polygons
  Start();
  glEnable(GL_TEXTURE_2D);
  glTexEnvf(GL_TEXTURE_ENV, GL_TEXTURE_ENV_MODE, GL_REPLACE);
  glBindTexture(GL_TEXTURE_2D, TextureNum);
  glBegin (GL_QUADS);
  Vec3 u1 = u0 + sVec;
  Vec3 u2 = u0 + sVec + tVec;
  Vec3 u3 = u0 + tVec;
  Vec3 r0 = u0*Lattice;
  Vec3 r1 = u1*Lattice;
  Vec3 r2 = u2*Lattice;
  Vec3 r3 = u3*Lattice;
  // When commensurate with the face of the cell, place slightly
  // outside the cell.
  if ((fabs(Position)<1.0e-6) || (fabs(1.0-Position)<1.0e-6)) {
    r0 = 1.0001*r0;  r1 = 1.0001*r1;  r2 = 1.0001*r2;  r3 = 1.0001*r3;
  }
  glTexCoord2f(0.0, 0.0); glVertex3f(r0[0], r0[1], r0[2]);
  glTexCoord2f(0.0, 1.0); glVertex3f(r1[0], r1[1], r1[2]);
  glTexCoord2f(1.0, 1.0); glVertex3f(r2[0], r2[1], r2[2]);
  glTexCoord2f(1.0, 0.0); glVertex3f(r3[0], r3[1], r3[2]);
  glEnd();
  glDisable(GL_TEXTURE_2D);
  End();
}
	

void
PlaneObject::DrawPOV(FILE *fout, string rotMatrix)
{
  // First, we have to create an image file for the texture map
  const int N = 1024;
  Array<TinyVector<guint8,4>,2> texData(N,N);
  
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
      texData(is, it)[0] = (guint8) floor (256.0 * color[0]);
      texData(is, it)[1] = (guint8) floor (256.0 * color[1]);
      texData(is, it)[2] = (guint8) floor (256.0 * color[2]);
      texData(is, it)[3] = (guint8) floor (256.0 * color[3]);
    }
  }
  
  Glib::RefPtr<Gdk::Pixbuf> pixbuf = 
    Gdk::Pixbuf::create_from_data((guint8*)texData.data(), 
				  Gdk::COLORSPACE_RGB,
				  true, 8, N, N, 4*N);
  stringstream fname;
  fname << "ColorPlane" << Direction << ".png";
  pixbuf->save (fname.str(), (Glib::ustring)"png"); 

  // Now, we actually create a box geometry onto which to apply this
  // texture.  
//   fprintf (fout, "intersection {\n");
//   fprintf (fout, "  box {\n");
//   fprintf (fout, "    <%10.8f, %10.8f, %10.8f>,\n",
// 	   Spline.Xgrid->Start, Spline.Ygrid->Start, Spline.Zgrid->Start);
//   fprintf (fout, "    <%10.8f, %10.8f, %10.8f>\n",
//  	   Spline.Xgrid->End, Spline.Ygrid->End, Spline.Zgrid->End);
//   fprintf (fout, "%s", rotMatrix.c_str());
//   fprintf (fout, "  }\n");
  fprintf (fout, "box {\n");
  fprintf (fout, "  <0.0000,0.0000,0.0000>, \n");
  fprintf (fout, "  <1.0000,1.0000,0.0001>  \n");
  fprintf (fout, "  pigment {\n");
  fprintf (fout, "    image_map {\"%s\"}\n", fname.str().c_str());
  fprintf (fout, "  }\n"); // pigment
  fprintf (fout, "  translate <-0.5, -0.5, 0.0>\n");
  if (Direction == 0)
    fprintf (fout, "  scale <-1, -1, 1>\n");
  else if (Direction == 1)
    fprintf (fout, "  scale <1, -1, 1>\n");
  else
    fprintf (fout, "  scale <1, 1, 1>\n");
  // Now rotate
  if (Direction == 0)
    fprintf (fout, "  rotate <0, 90, 0>\n");
  else if (Direction == 1)
    fprintf (fout, "  rotate <90, 0, 0>\n");
  else if (Direction == 2)
    fprintf (fout, "  rotate <0, 0, 90>\n");

  // First scale
  fprintf (fout, "  scale <%1.6f, %1.6f, %1.6f>\n",
	   (Spline.Xgrid->End-Spline.Xgrid->Start),
	   (Spline.Ygrid->End-Spline.Ygrid->Start),
	   (Spline.Zgrid->End-Spline.Zgrid->Start));
  fprintf (fout, "  matrix <%12.8f, %12.8f, %12.8f, \n",
	   Lattice(0,0), Lattice(0,1), Lattice(0,2));
  fprintf (fout, "          %12.8f, %12.8f, %12.8f, \n",
	   Lattice(1,0), Lattice(1,1), Lattice(1,2));
  fprintf (fout, "          %12.8f, %12.8f, %12.8f, \n",
	   Lattice(2,0), Lattice(2,1), Lattice(2,2));
  fprintf (fout, "          %12.8f, %12.8f, %12.8f> \n",
	   0.0, 0.0, 0.0);

  fprintf (fout, "  translate <%1.6f, %1.6f, %1.6f>\n",
	   (Position-0.5) * Lattice(Direction,0), 
	   (Position-0.5) * Lattice(Direction,1), 
	   (Position-0.5) * Lattice(Direction,2));
  

  
//   if (Direction == 0) 
//     fprintf (fout, "  translate <%1.6f, %1.6f, %1.6f>\n",
//  	     (1.0-Position)*Spline.Xgrid->Start
//  	     +Position*Spline.Xgrid->End, 0.0, 0.0);
//   else if (Direction == 1)
//      fprintf (fout, "  translate <%1.6f, %1.6f, %1.6f>\n",
//  	     0.0, 
//  	     (1.0-Position)*Spline.Ygrid->Start
//  	     +Position*Spline.Ygrid->End, 0.0);
//   else if (Direction == 2)
//     fprintf (fout, "  translate <%1.6f, %1.6f, %1.6f>\n",
//  	     0.0, 0.0, (1.0-Position)*Spline.Zgrid->Start
//  	     + Position*Spline.Zgrid->End);
  


  // Now translate
  //   fprintf (fout, "  translate <%1.6f, %1.6f, %1.6f>\n",
  // 	   Spline.Xgrid->Start, Spline.Ygrid->Start, Spline.Zgrid->Start);

  fprintf (fout, "%s", rotMatrix.c_str());
  fprintf (fout, "  }\n"); // box
  //  fprintf (fout, "}\n");   // intersection
}


void
PlaneObject::SetLattice (Mat3 lattice)
{
  Lattice = lattice;
}
