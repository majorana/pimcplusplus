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
	if (vec.G2 < (kcut*kcut)) 
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
