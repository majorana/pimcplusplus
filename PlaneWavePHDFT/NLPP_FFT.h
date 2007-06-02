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

#ifndef NLPP_FFT_H
#define NLPP_FFT_H

#include "HamiltonianBase.h"
#include "../PH/NLPP.h"
#include "../PH/CoulombPot.h"
#include "../PH/SplinePot.h"

// Implements the King-Smith method for applying nonlocal
// pseudopotentials in real space.  See PRB 44, 13063.

class KingSmithProjector
{
public:
  int l;
  double R0;
  CubicSpline chi_r;
  double El;

  void Setup (int l, CubicSpline &chi, double R0,
	      double El);
};


// This class stores what is necessary to do Kleinman-Bylander
// projection in real space on a single ion for one l channel.
class Ion_l_Projector
{
private:
  FFTBox *FFT;
  complex<double> Ylm (int l, int m, Vec3 r);
  complex<double> Ylm2(int l, int m, Vec3 r);
public:
  Array<Int3, 1> FFTIndices;
  // The first index specifies the point.  The second index specifies
  // the value of m.  Since chi_r is real and 
  // Y_l{-m} = (-1)^m conj(Y_lm), we only need to store for values of
  // m from 0 to l.  We will store for all -l to l, though.
  Array<complex<double>,2> ChiYlm;
  // The volume of a mesh element
  double MeshVol;
  void Setup(KingSmithProjector &projector,
	     Vec3 rion, FFTBox *fft);
  // This returns the application of the projector to the contents of
  // the FFTBox in chi_psi.  Thus, for l=0, chi_psi has one element,
  // for l=1: 3 elements, l=2:  5 elements...
  void Project (Array<complex<double>,1> &chi_psi);
};



class NLPP_FFTClass : public VionBase
{
private:
  // This computes the reciprocal-space version of the local potential 
  kSpacePH kPH;
  zVec Vc;
  FFTBox      &cFFT;
  void SetupkPotentials();
  void SetuprPotentials();
  NLPPClass NLPP;
  SplinePot Vlocal;
  CoulombPot Vouter;
  Array<Ion_l_Projector,1> Ion_l_Projectors;

public:
  void Setup();
  void SetIons (const Array<Vec3, 1>& rions);
  void Vmatrix (Array<complex<double>,2> &vmat);
  void Setk(Vec3 k);

  void Apply (const zVec &c, zVec &Hc);
  void Apply (const zVec &c, zVec &Hc,
	      Array<double,3> &VHXC);
  NLPP_FFTClass (NLPPClass &nlpp, GVecsClass &gvecs, FFTBox &fft) :
    VionBase (gvecs), kPH(Vlocal), cFFT(fft), NLPP(nlpp)
  {

  }
};

#endif
