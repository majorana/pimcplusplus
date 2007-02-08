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

#ifndef PAIR_FIXED_PHASE_CLASS_H
#define PAIR_FIXED_PHASE_CLASS_H

#include "NodalActionClass.h"
#include <Common/PlaneWavePHDFT/PlaneWavesMPI.h>
#include <Common/Splines/ComplexMultiTricubicSpline.h>
#include "../MirroredClass.h"

class PairFixedPhaseClass  : public ActionBaseClass
{
public:
  Array<complex<double>,2> DetMatrix;
  void BuildRefDet();
  int Species1Num;
  int Species2Num;
  int NumberAssymetry;
  PathDataClass &PathData;
#ifdef BUILD_DEV
  PathClassDev &Path;
#else
  PathClass &Path;
#endif
  double SingleAction (int slice1, int slice2,
		 const Array<int,1> &activeParticles, 
		 int level);
  PairFixedPhaseClass (PathDataClass &pathData);
  double d_dBeta(int slice1, int slice2, int level);
  string GetName();
  complex<double> Det (int slice,dVec kVec);
  Array<complex<double>,1> DegenerateRefSliceDeterminates;
  complex<double>  CalcDegenerateDet(int slice);

  void Read (IOSectionClass &in);
};

#endif
