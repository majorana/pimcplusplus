/////////////////////////////////////////////////////////////
// Copyright (C) 2003-2006 Bryan Clark and Kenneth Esler   //
//                                                         //
// This program is free software; you can redistribute it  //
// and/or modify it under the terms of the GNU General     //
// Public License as published by the Free Software        //
// Foundation; either version 2 of the License, or         //
// (at your option) any later version.  This program is    //
// distributed in the hope that it will be useful, but     //
// WITHOUT ANY WARRANTY; without even the implied warranty //
// of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. //  
// See the GNU General Public License for more details.    //
// For more information, please see the PIMC++ Home Page:  //
//           http://pathintegrals.info                     //
/////////////////////////////////////////////////////////////

#include "GVecs.h"
#include <vector>
#include <algorithm>

class GVec
{
public:
  Vec3 G;
  Int3 I;
  double G2;
};

inline bool operator<(const GVec &g1, const GVec &g2)
{
  return (g1.G2 < g2.G2);
}

inline bool operator==(const GVec &g1, const GVec &g2)
{
  return (fabs(dot(g1.G,g1.G)-dot(g2.G,g2.G)) < 1.0e-12);
}


Int3
GVecsClass::GetFFTBoxSize(Vec3 box, Vec3 kvec, double kcut)
{
  int maxX, maxY, maxZ;
  Vec3 kbox = Vec3(2.0*M_PI/box[0], 2.0*M_PI/box[1], 2.0*M_PI/box[2]);
  maxX = (int) ceil(kcut/kbox[0])+1;
  maxY = (int) ceil(kcut/kbox[1])+1;
  maxZ = (int) ceil(kcut/kbox[2])+1;

  // The FFT box must be twice the size as the maximum G in each direction.
  int nx = 4*maxX+1;
  int ny = 4*maxY+1;
  int nz = 4*maxZ+1; 

  ////////////////////////////////////////////
  // First, set up G-vectors and difference //
  ////////////////////////////////////////////
  int actxmin, actxmax, actymin, actymax, actzmin, actzmax;
  actxmin=0; actxmax=0;
  actymin=0; actymax=0; 
  actzmin=0; actzmax=0;
  TinyVector<double,3> g;
  for (int ix=-2*maxX; ix<=2*maxX; ix++) {
    g[0] = ix*kbox[0];
    for (int iy=-2*maxY; iy<=2*maxY; iy++) {
      g[1] = iy*kbox[1];
      for (int iz=-2*maxZ; iz<=2*maxZ; iz++) {
	g[2] = iz*kbox[2];
	if (dot(kvec+g,kvec+g) < (4.0*kcut*kcut)) {
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
  nx = actxmax-actxmin+1;
  ny = actymax-actymin+1;
  nz = actzmax-actzmin+1;
  if ((nx%2)==1) nx++;
  if ((ny%2)==1) ny++;
  if ((nz%2)==1) nz++;
  return Int3(nx, ny, nz);
}


void GVecsClass::Set (Vec3 box, Vec3 kVec, double kcut)
{
  Int3 boxSize = GetFFTBoxSize(box, kVec, kcut);
  Set (box, kVec, kcut, boxSize);
}
			       

void GVecsClass::Set (Vec3 box, Vec3 kVec, double kcut, Int3 boxSize)
{
  k = kVec;
  Box = box;
  Lattice = 
    box[0], 0.0   , 0.0,
    0.0   , box[1], 0.0,
    0.0   , 0.0   , box[2];

  kCut = kcut;
  kBox[0]=2.0*M_PI/box[0]; kBox[1]=2.0*M_PI/box[1]; kBox[2]=2.0*M_PI/box[2];

  int maxX, maxY, maxZ;
  maxX = (int) ceil(kcut/kBox[0]);
  maxY = (int) ceil(kcut/kBox[1]);
  maxZ = (int) ceil(kcut/kBox[2]);

  Nx = boxSize[0];
  Ny = boxSize[1];
  Nz = boxSize[2];

  ////////////////////////////////////////////
  // First, set up G-vectors and difference //
  ////////////////////////////////////////////
  vector<GVec> vecs;
  GVec vec;

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
	vec.G2 = dot (vec.G,vec.G);
	//	if (vec.G2 < (kcut*kcut)) 
	if (dot(vec.G+kVec,vec.G+kVec) < (kcut*kcut))
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
  vecs.clear();
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
	vec.G2 = dot (vec.G, vec.G);
	if (vec.G2 < (4.0*kcut*kcut))
	  vecs.push_back(vec);
      }
    }
  }
  sort (vecs.begin(), vecs.end());
  GDiff.resize(vecs.size());
  GDiffInv2.resize(vecs.size());
  IDiff.resize(vecs.size());
  numUnique = 1;
  for (int i=0; i<vecs.size(); i++) {
    GDiff(i) = vecs[i].G;
    if (vecs[i].I != TinyVector<int,3>(0,0,0))
      GDiffInv2(i) = 1.0/dot(vecs[i].G, vecs[i].G);
    else
      GDiffInv2(i) = 0.0;
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


GVecsClass& 
GVecsClass::operator=(const GVecsClass &gvecs)
{
  GVecs.resize(gvecs.GVecs.shape());
  GVecs = gvecs.GVecs;
  Indices.resize(gvecs.Indices.shape());
  Indices = gvecs.Indices;
  GDiff.resize(gvecs.GDiff.size());
  GDiff = gvecs.GDiff;
  GDiffInv2.resize(gvecs.GDiffInv2.shape());
  GDiffInv2 = gvecs.GDiffInv2;
  IDiff.resize (gvecs.IDiff.shape());
  IDiff = gvecs.IDiff;
  Box = gvecs.Box;
  kBox = gvecs.kBox;
  Nx = gvecs.Nx;  Ny = gvecs.Ny;  Nz = gvecs.Nz;
  kCut = gvecs.kCut;
  k = gvecs.k;

  return *this;
}

inline Int3 
GetIndex(Vec3 a[3], Vec3 b[3], Vec3 gvec)
{
  Int3 i;
  const double twoPiInv = 1.0/(2.0*M_PI);
  i[0] = (int)round(twoPiInv*dot(a[0], gvec));
  i[1] = (int)round(twoPiInv*dot(a[1], gvec));
  i[2] = (int)round(twoPiInv*dot(a[2], gvec));

  Vec3 gprime = (double)i[0]*b[0] + (double)i[1]*b[1] + (double)i[2]*b[2];
  Vec3 diff = gprime - gvec;
  assert (dot(diff,diff) < 1.0e-14);

  return i;
}
  


void
GVecsClass::Set (Mat3 &lattice, Array<Vec3,1> &gvecs, double fftFactor)
{
  Lattice = lattice;
  Vec3 a[3], b[3];
  a[0] = lattice(0,0), lattice(0,1), lattice(0,2);
  a[1] = lattice(1,0), lattice(1,1), lattice(1,2);
  a[2] = lattice(2,0), lattice(2,1), lattice(2,2);
  double vol = fabs(dot(cross(a[0],a[1]),a[2]));
  double volInv = 1.0/vol;
  b[0] = 2.0*M_PI*volInv * cross(a[1], a[2]);  
  b[1] = 2.0*M_PI*volInv * cross(a[2], a[0]);
  b[2] = 2.0*M_PI*volInv * cross(a[0], a[1]);

  cerr << "b0 = " << b[0] << endl;
  cerr << "b1 = " << b[1] << endl;
  cerr << "b2 = " << b[2] << endl;


  GVecs.resize(gvecs.size());
  Indices.resize(gvecs.size());
  GVecs = gvecs;
  
  for (int i=0; i<GVecs.size(); i++)
    Indices(i) = GetIndex(a, b, GVecs(i));

  /// Find FFT box size;
  Int3 maxIndex(0,0,0), minIndex(0,0,0);
  double maxG2 = 0.0;
  for (int i=0; i<Indices.size(); i++) {
    for (int j=0; j<3; j++) {
      if (Indices(i)[j] < minIndex[j])
	minIndex[j] = Indices(i)[j];
      if (Indices(i)[j] > maxIndex[j])
	maxIndex[j] = Indices(i)[j];
    }
    if (dot (GVecs(i), GVecs(i)) > maxG2)
      maxG2 = dot (GVecs(i), GVecs(i));
  }
  kCut = sqrt (maxG2);
//   cerr << "minIndex = " << minIndex << endl;
//   cerr << "maxIndex = " << maxIndex << endl;
//   cerr << "maxG2 = " << maxG2 << endl;
  Nx = 2*(maxIndex[0] - minIndex[0] + 1);
  Ny = 2*(maxIndex[1] - minIndex[1] + 1);
  Nz = 2*(maxIndex[2] - minIndex[2] + 1);

  // This is to increase the resolution of the FFT box, for spline
  // representation 
  Nx = (int)ceil(fftFactor * Nx);
  Ny = (int)ceil(fftFactor * Ny);
  Nz = (int)ceil(fftFactor * Nz);

  if ((Nx%2)==1) Nx++;
  if ((Ny%2)==1) Ny++;
  if ((Nz%2)==1) Nz++;

  /// Adjust indices for FFT box
  for (int i=0; i<Indices.size(); i++) {
    Indices(i)[0] = (Indices(i)[0]+Nx)%Nx;
    Indices(i)[1] = (Indices(i)[1]+Ny)%Ny;
    Indices(i)[2] = (Indices(i)[2]+Nz)%Nz;
  }
    
  Int3 iDiff;
  Vec3 gDiff;
  int numDiff = 0;
  // Construct difference vectors
  for (int ix=2*minIndex[0]; ix<=2*maxIndex[0]; ix++) {
    iDiff[0] = (ix+Nx)%Nx;
    for (int iy=2*minIndex[1]; iy<=2*maxIndex[1]; iy++) {
      iDiff[1] = (iy+Ny)%Ny;
      for (int iz=2*minIndex[2]; iz<=2*maxIndex[2]; iz++) {
	iDiff[2] = (iz+Nz)%Nz;
	gDiff   = (double)ix*b[0]+(double)iy*b[1]+(double)iz*b[2];
	if (dot (gDiff, gDiff) <= 4.0000001*maxG2)
	  numDiff++;
      }
    }
  }
  GDiff.resize(numDiff);
  GDiffInv2.resize(numDiff);
  IDiff.resize(numDiff);
  cerr << "Found " << GDiff.size() << " unique difference vectors.\n";

  // Now actually store
  numDiff = 0;
  for (int ix=2*minIndex[0]; ix<=2*maxIndex[0]; ix++) {
    iDiff[0] = (ix+Nx)%Nx;
    for (int iy=2*minIndex[1]; iy<=2*maxIndex[1]; iy++) {
      iDiff[1] = (iy+Ny)%Ny;
      for (int iz=2*minIndex[2]; iz<=2*maxIndex[2]; iz++) {
	iDiff[2] = (iz+Nz)%Nz;
	gDiff   = (double)ix*b[0]+(double)iy*b[1]+(double)iz*b[2];
	if (dot (gDiff, gDiff) <= maxG2) {
	  GDiff(numDiff) = gDiff;
	  IDiff(numDiff) = iDiff;
	  GDiffInv2(numDiff) = 1.0/dot(gDiff, gDiff);
	  numDiff++;
	}
      }
    }
  }


//   for (int g1=0; g1<GVecs.size(); g1++) {
//     Int3 i1 = Indices(g1);
//     Indices(g1) = i1;
//     for (int g2=0; g2<GVecs.size(); g2++) {
//       Int3 i2 = Indices(g1);
//       Int3 idiff = i2-i1;
//       bool duplicate = false;
//       for (int j=0; j<IDiff.size(); j++) 
// 	if (idiff == IDiff(j))
// 	  duplicate = true;
//       if (!duplicate) 
// 	diffCount ++;
//     }
//   }
//   IDiff.resize(diffCount);
//   GDiff.resize(diffCount);




//   diffCount = 0;
//   for (int g1=0; g1<GVecs.size(); g1++) {
//     Int3 i1 = Indices(g1);
//     Indices(g1) = i1;
//     for (int g2=0; g2<GVecs.size(); g2++) {
//       Int3 i2 = Indices(g1);
//       Int3 idiff = i2-i1;
//       bool duplicate = false;
//       for (int j=0; j<IDiff.size(); j++) 
// 	if (idiff == IDiff(j))
// 	  duplicate = true;
//       if (!duplicate) {
// 	IDiff(diffCount) = idiff;
// 	Vec3 G = (double)idiff[0]*b[0] + (double)idiff[1]*b[1] + 
// 	  (double)idiff[2]*b[2];
// 	GDiff(diffCount) = G;
//       }
//     }
//   }

  
//   for (int i=0; i<IDiff.size(); i++) {
//     IDiff(i)[0] = (IDiff(i)[0]+Nx)%Nx;
//     IDiff(i)[1] = (IDiff(i)[1]+Ny)%Ny;
//     IDiff(i)[2] = (IDiff(i)[2]+Nz)%Nz;
//   }
    
}
