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

#ifndef KINETIC_SPHERE_CLASS_H
#define KINETIC_SPHERE_CLASS_H

#include "ActionBase.h"

/// The KineticSphereClass calculates the kinetic part of the action.  This
/// is the "spring term", of the form
/// \f$ K \equiv \left(4\pi\lambda\tau\right)^{-\frac{ND}{2}} 
/// \exp\left[-\frac{(R-R')^2}{4\lambda\tau}\right] \f$
class KineticSphereClass : public ActionBaseClass
{
  int NumImages;
  double SphereRadius;
public:
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		 const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  double d_dBeta_old (int slice1, int slice2, int level);
  inline void SetNumImages (int num) { NumImages = num; }
  KineticSphereClass (PathDataClass &pathData);
  double K(int slice,int nextSlice,int ptcl,int level,double lambda);
  double K2(int slice,int nextSlice,int ptcl,int level,double lambda);

};


#endif
