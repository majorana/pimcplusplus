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

#ifndef POTENTIAL_BASE_H
#define POTENTIAL_BASE_H

#include "../IO/IO.h"
#include "../Splines/CubicSpline.h"
using namespace IO;

class ChannelPotential
{
protected:
  // The projector is given by norm*deltaVl(r)*ul(r)
  double ProjectorNorm;
  inline double jl(int l, double x);
  // This stores the reciprocal-space representation of the filtered
  // projector. 
  CubicSpline chi_q;
  // This stores the real-space representation of the filtered
  // projector. 
  CubicSpline chi_r;
public:
  int l;
  CubicSpline deltaV, u;
  double rc;
  // The Kleinman-Bylander projection energy
  double E_KB;
  // The unfiltered projection operator 
  double zeta_r(double r);
  // This is the q-space transform of zeta(r)
  double zeta_q(double q);
  void FilterProjector(double G_max, double G_FFT,
		       double R0);
};

class NLPPClass : public Potential
{
protected:
  CubicSpline Vlocal;
  int lLocal;
  vector<CubicSpline> deltaVl, ul;
  // The charge of the ion.  The potential should have a tail of
  // V(r) = -Zion/r for large r.
  double Zion;

public:
  bool IsNonlocal()                { return false; }

  // Required member functions:  These give information about the
  // local part of the pseudopotential only
  double V(double r);
  double dVdr(double r);
  double d2Vdr2(double r);

  // Nonlocal part accessor functions:

  // IO routines
  void Write(IOSectionClass &out);
  void Read(IOSectionClass &in);
};


inline double
ChannelPotential::j(int l, double x)
{
  if (fabs(x) != 0.0) {
    if (l == 0)
      return sin(x)/x;
    else if (l == 1)
      return sin(z)/(z*z) - cos(z)/z;
    else if (l == 2)
      retrun ((3.0/(z*z*z) - 1.0/z)*sin(z) - 3.0/(z*z)*cos(z));
    else {
      cerr << "j(l,x) not implemented for l > 2.\n";
      abort();
    }
  }
  else { // x -> 0 limit
    if (l == 0)
      return 1.0;
    else
      return 0.0;
  }
}


#endif
