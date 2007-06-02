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

#include "NLPP_FFT.h"


////////////////////////////////////////////////////////////
//                  IonProjector stuff                    //
////////////////////////////////////////////////////////////

complex<double> 
Ion_l_Projector::Ylm(int l, int m, Vec3 r)
{
  if (l == 0)
    return complex<double>(1.0/sqrt(4.0*M_PI), 0.0);

  if (m < 0) {
    double sign = 1.0;
    for (int i=0; i<m; i++)
      sign *= -1.0;
    return sign * conj(Ylm(l, -m, r));
  }

  Vec3 omega = 1.0/ sqrt(dot(r, r)) * r;
  double costheta = omega[2];
  double sintheta = sqrt(1.0-costheta*costheta);
  double cosphi, sinphi;
      
  if (sintheta != 0.0) {
    double sinthetaInv = 1.0/sintheta;
    cosphi = sinthetaInv*omega[0];
    sinphi = sinthetaInv*omega[1];
  }
  else { // We're at the poles, phi is undefined
    cosphi = 1.0;
    sinphi = 0.0;
  }
  if (l == 1) {
    if (m==0)
      return -sqrt(3.0/(4.0*M_PI)) * costheta *	complex<double>(1.0, 0.0);
    else if (m==1)
      return -sqrt(3.0/(8.0*M_PI))*sintheta * 
	complex<double> (cosphi, sinphi);
    else {
      cerr << "Invalid l,m\n";
      abort();
    }
  }
  else if (l==2) {
    if (m==0)
      return sqrt(5.0/(4.0*M_PI)) * (1.5*costheta *costheta - 0.5) *
	complex<double>(1.0, 0.0);
    else if (m == 1)
      return -sqrt(15.0/(8.0*M_PI))*sintheta*costheta*
	complex<double>(cosphi, sinphi);
    else if (m==2){
      double cos2phi = cosphi*cosphi - sinphi*sinphi;
      double sin2phi = 2.0*sinphi*cosphi;
      return 0.25*sqrt(15.0/(2.0*M_PI))*sintheta*sintheta *
	complex<double>(cos2phi, sin2phi);
    }
  }
  else {
    cerr << "Ylm not implemented for l > 2.\n";
    abort();
  }
  return 0.0;
}


complex<double> 
Ion_l_Projector::Ylm2(int l, int m, Vec3 r)
{
  if (l == 0)
    return 0.5*sqrt(1.0/M_PI);

  if (m < 0) {
    double sign = ((-m)&1) ? -1.0 : 1.0;
    return sign * conj(Ylm(l, -m, r));
  }

  double nrm = sqrt (1.0/dot(r, r));
  double x = nrm*r[0];
  double y = nrm*r[1];
  double z = nrm*r[2];
  if (l==1) {
    if (m==0) 
      return 0.5*sqrt(3.0/M_PI)*z;
    else if (m==1)
      return -0.5*sqrt(3.0/(2.0*M_PI))*complex<double>(x,y);
  }
  else if (l==2) {
    if (m==0)
      return 0.25*sqrt(5.0/M_PI)*(3.0*z*z - 1.0);
    else if (m==1)
      return -0.5*sqrt(15.0/(2.0*M_PI))*z*complex<double>(x,y);
    else if (m==2) {
      complex<double> xy(x,y);
      return 0.25*(15.0/(2.0*M_PI))*xy*xy;
    }
  }
  else {
    cerr << "Ylm not implemented for l > 2.\n";
    abort();
  }
  return 0.0;
}


void
Ion_l_Projector::Setup(KingSmithProjector &projector,
		       Vec3 rion, FFTBox *fft)
{
  FFT = fft;
  double R0 = projector.R0;
  int nx, ny, nz;
  fft->GetDims (nx, ny, nz);
  Vec3 box = fft->GVecs.GetBox();
  Vec3 boxInv(1.0/box[0], 1.0/box[1], 1.0/box[2]);
  double dx = box[0]/(double)nx;
  double dy = box[0]/(double)ny;
  double dz = box[0]/(double)nz;
  Vec3 r;
  vector<Int3> indices;
  vector<Vec3> disps;
  for (int ix=0; ix<nx; ix++) {
    r[0] = dx * ix;
    for (int iy=0; iy<ny; iy++) {
      r[1] = dy * iy;
      for (int iz=0; iz<nz; iz++) {
	r[2] = dz * iz;
	Vec3 disp = r - rion;
	disp[0] -= round(disp[0]*boxInv[0])*box[0];
	disp[1] -= round(disp[1]*boxInv[1])*box[1];
	disp[2] -= round(disp[2]*boxInv[2])*box[2];
	if (dot(disp, disp) < R0*R0) {
	  // This point is inside the projection core
	  indices.push_back(Int3(ix,iy,iz));
	  disps.push_back(disp);
	}
      }
    }
  }
  cerr << "Found " << indices.size() << " projection points.\n";
  // Now, compute projector inside the core
  int l = projector.l;
  ChiYlm.resize(indices.size(), 2*l+1);
  FFTIndices.resize(indices.size());
  for (int i=0; i<indices.size(); i++) {
    FFTIndices(i) = indices[i];
    double r = sqrt(dot(disps[i], disps[i]));
    for (int m=-l; m<=l; m++)
      ChiYlm(i,l+m) = Ylm(l,m,disps[i]) * projector.chi_r(r);
  }
  MeshVol = dx*dy*dz;
}



// Note:  this presently only includes the local potential
void
NLPP_FFTClass::Vmatrix (Array<complex<double>,2> &vmat)
{
  if (!IsSetup)
    Setup();
  double volInv = 1.0/GVecs.GetBoxVol();
  for (int i=0; i<vmat.rows(); i++) 
    for (int j=0; j<=i; j++) {
      Vec3 diff = GVecs(i) - GVecs(j);
      complex<double> s(0.0,0.0);
      for (int zi=0; zi<Rions.size(); zi++) {
	double cosVal, sinVal, phase;
	phase = dot (diff, Rions(zi));
	sincos(phase, &sinVal, &cosVal);
	s += complex<double> (cosVal,sinVal);
      }
      vmat (i,j) = s*kPH.V(kPoint, GVecs(i), GVecs(j))*volInv;
      vmat (j,i) = conj(vmat(i,j));
    }
}


void 
NLPP_FFTClass::SetupkPotentials()
{
  // First, create Vlocal from the local potential of the NLPP
  Vlocal.Spline = NLPP.GetLocalSpline();
  Vouter.Z1Z2   = -NLPP.GetValenceCharge();
  Vlocal.Vouter = &Vouter;

  // Now, compute the tail coefficients of the local potential
  kPH.CalcTailCoefs (30.0, 60.0);

  // Compute local potential in reciprocal space
  double volInv = 1.0/GVecs.GetBoxVol();
  // Setup V and F tensors in k-space
  //  double gMag, lastMag2;
  double gMag, lastMag2;
  lastMag2 = -1.0;
  double a, bPerp, bPar, V;
  int numCalls = 0;
  for (int i=0; i<GVecs.DeltaSize(); i++) {
    double gMag2 = dot(GVecs.DeltaG(i), GVecs.DeltaG(i));
    if (fabs(gMag2-lastMag2) > 1.0e-12) {
      lastMag2 = gMag2;
      gMag  = sqrt(lastMag2);
      kPH.GetVals(gMag, a, bPerp, bPar, V);
      numCalls++;
    }
    if (gMag2 > 1.0e-10)
      VG(i) = volInv * V;
    else {
      VG(i) = 0.0;
      VG0 = volInv * V;
    }
  }
}

void
NLPP_FFTClass::SetuprPotentials()
{
  cFFT.kBox = complex<FFT_FLOAT>(0.0, 0.0);
  for (int i=0; i<GVecs.DeltaSize(); i++) {
    Int3 I = GVecs.DeltaI(i);
    cFFT.kBox(I)   = StructureFactor(i)*VG(i);
  }
  cFFT.k2r();
  Vr = cFFT.rBox;
}

void
NLPP_FFTClass::SetIons(const Array<Vec3,1> &rions)
{
  // Calculate the structure factor
  VionBase::SetIons(rions);
  if (IsSetup)
    SetuprPotentials();
}



void
NLPP_FFTClass::Setup()
{
  int nx, ny, nz;
  cFFT.GetDims(nx,ny,nz);
  Vr.resize(nx,ny,nz);
  Vc.resize(GVecs.size());
  VG.resize(GVecs.DeltaSize());
  SetupkPotentials();
  SetuprPotentials();

  // Compute the Kleinmain-Bylander projectors:
  double kc = cFFT.GVecs.GetkCut();
  NLPP.SetupProjectors(kc, 4.0*kc);

  IsSetup = true;
}

void
NLPP_FFTClass::Setk (Vec3 k)
{
  kPoint = k;
  
  int nx, ny, nz;
  cFFT.GetDims(nx,ny,nz);
  Vr.resize(nx,ny,nz);
  Vc.resize(GVecs.size());
  VG.resize(GVecs.DeltaSize());

  SetupkPotentials();
  SetIons(Rions);
}


void 
NLPP_FFTClass::Apply (const zVec &c, zVec &Hc)
{
  if (!IsSetup)
    Setup();
  int nx, ny, nz;
  cFFT.GetDims(nx, ny, nz);

  ////////////////////
  // Potential part //
  ////////////////////
  // Transform c into real space
  cFFT.PutkVec (c);
  cFFT.k2r();
  // Multiply by V
  cFFT.rBox *= Vr;
  // Transform back
  cFFT.r2k();

  // Get vector
  cFFT.GetkVec (Vc);
  for (int i=0; i<GVecs.size(); i++)
    Hc(i) += Vc(i);
}

void 
NLPP_FFTClass::Apply (const zVec &c, zVec &Hc,
		      Array<double,3> &VHXC)
{
  if (!IsSetup)
    Setup();
  int nx, ny, nz;
  cFFT.GetDims(nx, ny, nz);

  ////////////////////
  // Potential part //
  ////////////////////
  // Transform c into real space
  cFFT.PutkVec (c);
  cFFT.k2r();

  // Apply local potential and VHXC
  //  cFFT.rBox *= (Vr+VHXC);
  for (int ix=0; ix<nx; ix++)
    for (int iy=0; iy<ny; iy++)
      for (int iz=0; iz<nz; iz++)
	cFFT.rBox(ix,iy,iz) *= (Vr(ix,iy,iz)+VHXC(ix,iy,iz));

  // Transform back
  cFFT.r2k();

  // Get vector
  cFFT.GetkVec (Vc);
  for (int i=0; i<GVecs.size(); i++)
    Hc(i) += Vc(i);
}
