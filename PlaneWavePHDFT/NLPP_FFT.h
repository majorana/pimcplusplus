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
