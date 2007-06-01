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

#ifndef NLPP_CLASS_H
#define NLPP_CLASS_H

#include "PotentialBase.h"
#include "../IO/IO.h"
#include "../Splines/CubicSpline.h"
#include <vector>
using namespace IO;
using namespace blitz;

class ChannelPotential
{
protected:
  // The projector is given by norm*deltaVl(r)*ul(r)
  double ProjectorNorm;
  inline double jl(int l, double x);
  LinearGrid qGrid;
  double qCurrent, rCurrent;
  typedef enum { NORM, EKB, ZETA_Q, CHI_R, CHECK_CHI_R} IntegrandType;
  IntegrandType Job;
  double A(double q, double qp);
public:
  int l;
  // V stores the potential
  // DeltaV stores the potential less to the local potential
  CubicSpline V, DeltaV, u;
  double rc, R0;
  // The Kleinman-Bylander projection energy
  double E_KB;
  // The unfiltered projection operators in real-space and reciprocal
  // space.  
  CubicSpline zeta_r, zeta_q;
  // These store the real and reciprocal-space representation of the
  // filtered projector
  CubicSpline chi_r, chi_q;

  // This is used as the integrand for the various integrals that are required.
  inline double operator()(double x);

  void SetupProjector(double G_max, double G_FFT);
};

class NLPPClass : public Potential
{
protected:
  int lLocal;
  vector<ChannelPotential> Vl;  
  // The charge of the ion.  The potential should have a tail of
  // V(r) = -Zion/r for large r.
  double Zion;
  int AtomicNumber;
public:
  bool IsNonlocal();
  inline int LocalChannel() { return lLocal; }
  inline int NumChannels()  { return Vl.size(); }

  // Required member functions:  These give information about the
  // local part of the pseudopotential only
  double V(double r);
  double dVdr(double r);
  double d2Vdr2(double r);

  // Nonlocal part accessor functions:

  // IO routines
  void Write(IOSectionClass &out);
  void Read(IOSectionClass &in);
  void SetupProjectors(double G_max, double G_FFT);
};


inline double
ChannelPotential::jl(int l, double x)
{
  if (fabs(x) != 0.0) {
    if (l == 0)
      return sin(x)/x;
    else if (l == 1)
      return sin(x)/(x*x) - cos(x)/x;
    else if (l == 2)
      return ((3.0/(x*x*x) - 1.0/x)*sin(x) - 3.0/(x*x)*cos(x));
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


inline double
ChannelPotential::operator()(double x)
{
  switch (Job) {
  case NORM:
    return u(x)*DeltaV(x)*DeltaV(x)*u(x);
  case EKB:
    return u(x)*DeltaV(x)*u(x);
  case ZETA_Q:
    return jl(l,qCurrent*x)*x*x*zeta_r(x);
  case CHI_R:
    return 2.0/M_PI * x*x*chi_q(x)*jl(l,x*rCurrent);
  case CHECK_CHI_R:
    return chi_r(x)*chi_r(x)*x*x;
  default:
    return 0.0;
  }
}



#endif
