#include "GVecs.h"
#include <vector>
#include <algorithm>

class GVec
{
public:
  Vec3 G;
  Int3 I;
};

inline bool operator<(const GVec &g1, const GVec &g2)
{
  return (dot(g1.G, g1.G) < dot(g2.G, g2.G));
}

inline bool operator==(const GVec &g1, const GVec &g2)
{
  return (fabs(dot(g1.G,g1.G)-dot(g2.G,g2.G)) < 1.0e-12);
}

void GVecsClass::Set (Vec3 box, Vec3 kVec, double kcut)
{
  k = kVec;
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
  vector<GVec> vecs;
  GVec vec;
  int actxmin, actxmax, actymin, actymax, actzmin, actzmax;
  actxmin = actxmax = actymin = actymax = actzmin = actzmax = 0;
  TinyVector<double,3> g;
  for (int ix=-2*maxX; ix<=2*maxX; ix++) {
    g[0] = ix*kBox[0];
    for (int iy=-2*maxY; iy<=2*maxY; iy++) {
      g[1] = iy*kBox[1];
      for (int iz=-2*maxZ; iz<=2*maxZ; iz++) {
	g[2] = iz*kBox[2];
	if (dot(k+g,k+g) < (4.0*kcut*kcut)) {
	  actxmin = min (ix, actxmin);
	  actxmax = max (ix, actxmax);
	  actymin = min (iy, actymin);
	  actymax = max (iy, actymax);
	  actzmin = min (iz, actzmin);
	  actzmax = max (iz, actzmax);
	}
      }
    }
  }
  Nx = actxmax-actxmin+2;
  Ny = actymax-actymin+2;
  Nz = actzmax-actzmin+2;

  /// First, count k-vectors
  for (int ix=-maxX; ix<=maxX; ix++) {
    vec.G[0] = ix*kBox[0];
    vec.I[0] = (ix+Nx)%Nx;
    for (int iy=-maxY; iy<=maxY; iy++) {
      vec.G[1] = iy*kBox[1];
      vec.I[1] = (iy+Ny)%Ny;
      for (int iz=-maxZ; iz<=maxZ; iz++) {
	vec.G[2] = iz*kBox[2];
	vec.I[2] = (iz+Nz)%Nz;
	if (dot(k+vec.G,k+vec.G) < (kcut*kcut)) 
	  vecs.push_back(vec);
      }
    }
  }
  sort (vecs.begin(), vecs.end());
  GVecs.resize(vecs.size());
  Indices.resize(vecs.size());
  int numUnique = 1;
  for (int i=0; i<vecs.size(); i++) {
    GVecs(i) = vecs[i].G;
    Indices(i) = vecs[i].I;
    if (i>0)
      if (!(vecs[i] == vecs[i-1]))
	numUnique++;
  }
  cerr << "Using " << vecs.size() << " G-vectors, of which " 
       << numUnique << " have unique magnitudes.\n";
  ////////////////////////////////////////////
  // Now, set up G-vector differences.      //
  ////////////////////////////////////////////
  vecs.resize(0);
  /// First, count k-vectors
  for (int ix=-2*maxX; ix<=2*maxX; ix++) {
    vec.G[0] = ix*kBox[0];
    vec.I[0] = (ix+Nx)%Nx;
    for (int iy=-2*maxY; iy<=2*maxY; iy++) {
      vec.G[1] = iy*kBox[1];
      vec.I[1] = (iy+Ny)%Ny;
      for (int iz=-2*maxZ; iz<=2*maxZ; iz++) {
	vec.G[2] = iz*kBox[2];
	vec.I[2] = (iz+Nz)%Nz;
	if (dot(k+vec.G,k+vec.G) < (4.0*kcut*kcut))
	  vecs.push_back(vec);
      }
    }
  }
  sort (vecs.begin(), vecs.end());
  GDiff.resize(vecs.size());
  IDiff.resize(vecs.size());
  numUnique = 1;
  for (int i=0; i<vecs.size(); i++) {
    GDiff(i) = vecs[i].G;
    IDiff(i) = vecs[i].I;
    if (i>0)
      if (!(vecs[i] == vecs[i-1]))
	numUnique++;
  }
  cerr << "Using " << vecs.size() << " G-difference vectors, of which "
       << numUnique << " have unique magnitudes.\n";
}


void GVecsClass::GetFFTBoxSize(int &nx, int &ny, int &nz)
{
  nx=Nx; ny=Ny; nz=Nz;
}
