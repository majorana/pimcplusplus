#include "GVecs.h"

void GVecsClass::Set (Vec3 box, double kcut)
{
  Box = box;
  kCut = kcut;
  kBox[0]=2.0*M_PI/box[0]; kBox[1]=2.0*M_PI/box[1]; kBox[2]=2.0*M_PI/box[2];

  int maxX, maxY, maxZ;
  maxX = (int) ceil(kcut/kBox[0]);
  maxY = (int) ceil(kcut/kBox[1]);
  maxZ = (int) ceil(kcut/kBox[2]);

  // The FFT box must be twice the size as the maximum G in each direction.
  Nx = 4*maxX+1;
  Ny = 4*maxY+1;
  Nz = 4*maxZ+1;

  ////////////////////////////////////////////
  // First, set up G-vectors and difference //
  ////////////////////////////////////////////
  int numG = 0;
  Vec3 G;
  /// First, count k-vectors
  for (int ix=-maxX; ix<=maxX; ix++) {
    G[0] = ix*kBox[0];
    for (int iy=-maxY; iy<=maxY; iy++) {
      G[1] = iy*kBox[1];
      for (int iz=-maxZ; iz<=maxZ; iz++) {
	G[2] = iz*kBox[2];
	if (dot(G,G) < (kcut*kcut))
	  numG++;
      }
    }
  }
  GVecs.resize(numG);
  Indices.resize(numG);
  numG = 0;
  /// Now assign
  for (int ix=-maxX; ix<=maxX; ix++) {
    G[0] = ix*kBox[0];
    for (int iy=-maxY; iy<=maxY; iy++) {
      G[1] = iy*kBox[1];
      for (int iz=-maxZ; iz<=maxZ; iz++) {
	G[2] = iz*kBox[2];
	if (dot(G,G) < (kcut*kcut)) {
	  GVecs(numG) = G;
	  Indices(numG) = Int3(ix,iy,iz);
	  numG++;
	}
      }
    }
  }
  cerr << "Using " << numG << " G-vectors.\n";
  ////////////////////////////////////////////
  // Now, set up G-vector differences.      //
  ////////////////////////////////////////////
  numG = 0;
  /// First, count k-vectors
  for (int ix=-2*maxX; ix<=2*maxX; ix++) {
    G[0] = ix*kBox[0];
    for (int iy=-2*maxY; iy<=2*maxY; iy++) {
      G[1] = iy*kBox[1];
      for (int iz=-2*maxZ; iz<=2*maxZ; iz++) {
	G[2] = iz*kBox[2];
	if (dot(G,G) < (4.0*kcut*kcut))
	  numG++;
      }
    }
  }
  GDiff.resize(numG);
  IDiff.resize(numG);
  numG = 0;
  /// Now assign
  for (int ix=-2*maxX; ix<=2*maxX; ix++) {
    G[0] = ix*kBox[0];
    for (int iy=-2*maxY; iy<=2*maxY; iy++) {
      G[1] = iy*kBox[1];
      for (int iz=-2*maxZ; iz<=2*maxZ; iz++) {
	G[2] = iz*kBox[2];
	if (dot(G,G) < (4.0*kcut*kcut)) {
	  GDiff(numG) = G;
	  IDiff(numG) = Int3(ix,iy,iz);
	  numG++;
	}
      }
    }
  }
}


void GVecsClass::GetFFTBoxSize(int &nx, int &ny, int &nz)
{
  nx=Nx; ny=Ny; nz=Nz;
}
