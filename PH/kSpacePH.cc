#include "kSpacePH.h"
#include "../Integration/GKIntegration.h"
#include "../Fitting/Fitting.h"
#include <gsl/gsl_sf_expint.h>


void kSpacePH::CalcTailCoefs (double r1, double r2)
{
  R1 = r1; R2 = r2;
  const int numPoints=1000;
  double nInv = 1.0/(double)(numPoints-1);

  Array<double,1> coefs, error;
  Array<double,1> y(numPoints), sigma(numPoints);
  Array<double,2> F(numPoints,3);
  
  // Setup basis functions:
  for (int i=0; i<numPoints; i++) {
    double r = r1 + (r2-r1)*(double)i*nInv;
    F(i,0) = 1.0/r;
    F(i,1) = 1.0/(r*r);
    F(i,2) = 1.0/(r*r*r);
    sigma(i) = 1.0;
    y(i) = PH.V(r);
    sigma(i) = 1.0;
  }
  LinFitSVD (y, sigma, F, coefs, error, 1.0e-15);
  Ctail1 = coefs(0);
  Ctail2 = coefs(1);
  Ctail3 = coefs(2);
  HaveTailCoefs = true;
}


class aIntegrand
{
protected:
  Potential &PH;
  double k, kinv;
public:
  inline void Setk (double newk) 
  {    k = newk;    kinv = 1.0/newk; }
   
  inline double operator()(double r)
  { 
    if (k < 1.0e-12)
      return 4.0*M_PI*(PH.A(r)-1.0)*r*r;
    else
      return 4.0*M_PI*(PH.A(r)-1.0)*r*sin(k*r)*kinv; 
  }

  aIntegrand (Potential &ph) : PH(ph) 
  { /* do nothing for now */ }
};

double kSpacePH::a(double k)
{
  aIntegrand integrand(PH);
  integrand.Setk(k);
  GKIntegration<aIntegrand,GK31> integrator(integrand);
  
  return (integrator.Integrate(0.0, PH.GetCoreRadius(), 1.0e-10));
}


class bPerpIntegrand
{
protected:
  Potential &PH;
  double k, kinv;
  const double Third, TwoThirds;
public:
  inline void Setk (double newk) 
  {    k = newk;    kinv = 1.0/newk; }
   
  inline double operator()(double r)
  { 
    double kr = k*r;
    double krInv = 1.0/kr;
    double b = PH.B(r) - 1.0;
    double j0, j2;
    if (k < 1.0e-12) {
      j0 = 1.0;
      j2 = 0.0;
    }
    else {
      j0 = sin(kr)/kr;
      j2 = (3.0*krInv*krInv*krInv-krInv)*sin(kr) - 
      3.0*krInv*krInv*cos(kr);
    }
    return (4.0*M_PI*r*r*b*(TwoThirds*j0-Third*j2));
  }
  
  bPerpIntegrand (Potential &ph) : PH(ph), Third(1.0/3.0), TwoThirds(2.0/3.0)
  { /* do nothing for now */ }
};


double kSpacePH::bPerp(double k)
{
  bPerpIntegrand integrand(PH);
  integrand.Setk(k);
  GKIntegration<bPerpIntegrand,GK31> integrator(integrand);
  
  return (integrator.Integrate(0.0, PH.GetCoreRadius(), 1.0e-10));
}



class bParIntegrand
{
protected:
  Potential &PH;
  double k, kinv;
  const double Third, TwoThirds;
public:
  inline void Setk (double newk) 
  {    k = newk;    kinv = 1.0/newk; }
   
  inline double operator()(double r)
  { 
    double kr = k*r;
    double krInv = 1.0/kr;
    double b = PH.B(r) - 1.0;
    double j0, j2;
    if (k < 1.0e-12) {
      j0 = 1.0;
      j2 = 0.0;
    }
    else {
      j0 = sin(kr)/kr;
      j2 = (3.0*krInv*krInv*krInv-krInv)*sin(kr) - 
      3.0*krInv*krInv*cos(kr);
    }
    return (4.0*M_PI*r*r*b*(TwoThirds*j0+TwoThirds*j2));
  }
  
  bParIntegrand (Potential &ph) : PH(ph), Third(1.0/3.0), TwoThirds(2.0/3.0)
  { /* do nothing for now */ }
};


double kSpacePH::bPar(double k)
{
  bParIntegrand integrand(PH);
  integrand.Setk(k);
  GKIntegration<bParIntegrand,GK31> integrator(integrand);
  
  return (integrator.Integrate(0.0, PH.GetCoreRadius(), 1.0e-10));
}


class VIntegrand
{
protected:
  Potential &PH;
  double k, kinv;
public:
  inline void Setk (double newk) 
  {    k = newk;    kinv = 1.0/newk; }
   
  inline double operator()(double r)
  { return 4.0*M_PI*PH.V(r)*r*sin(k*r)*kinv; }

  VIntegrand (Potential &ph) : PH(ph) 
  { /* do nothing for now */ }
};


double kSpacePH::Vk (double k)
{
  if (k == 0.0)
    return 0.0;
  // First, do the part of the integral up to R1 numerically
  VIntegrand integrand(PH);
  integrand.Setk(k);
  GKIntegration<VIntegrand,GK31> integrator(integrand);  
  double result = integrator.Integrate(0.0, R1, 1.0e-10);

  // Now, do the remaining part up to infinity using analytic
  // integratio of our fitted form.
  assert (HaveTailCoefs);
  
  result += 4.0*M_PI*Ctail1/(k*k) * cos(k*R1);
  result -= 4.0*M_PI*Ctail2/k*(gsl_sf_Si(k*R1) - 0.5*M_PI);
  result -= 4.0*M_PI*Ctail3/k*(k*gsl_sf_Ci(k*R1) - sin (k*R1)/R1);
  return result;
}


TinyMatrix<double,3,3> kSpacePH::Ftensor (Vec3 deltaG)
{
  double Gmag = sqrt(dot(deltaG, deltaG));
  Vec3 g; 
  if (Gmag == 0.0)
    g = Vec3 (1.0, 0.0, 0.0);
  else
    g = deltaG / Gmag;

  double aval, bPerpval, bParval, Vval;
  if (UseCache)
    Cache.GetVals (Gmag, aval, bPerpval, bParval, Vval);
  else {
    aval = a(Gmag);
    bPerpval = bPerp(Gmag);
    bParval = bPar(Gmag);
  }

  TinyMatrix<double,3,3> F, G;
  F(0,0)=0.0;   F(0,1)=0.0; F(0,2)=0.0;
  F(1,0)=0.0;   F(1,1)=0.0; F(1,2)=0.0;  
  F(2,0)=0.0;   F(2,1)=0.0; F(2,2)=0.0;
  for (int i=0; i<3; i++)
    F(i,i) += aval + bPerpval;
  for (int i=0; i<3; i++)
    for (int j=0; j<3; j++)
      G(i,j) = g[i]*g[j];
  F = F + (bParval - bPerpval)*G;
  return F;
}

double kSpacePH::V (Vec3 k, Vec3 G, Vec3 Gp)
{
  Vec3 deltaG = G-Gp;
  double Gmag = sqrt(dot(deltaG, deltaG));

  double aval, bPerpval, bParval, Vval;
  if (UseCache) 
    Cache.GetVals (Gmag, aval, bPerpval, bParval, Vval);
  else
    Vval = Vk(Gmag);

  TinyMatrix<double,3,3> F = Ftensor (deltaG);
//   if (G == Gp)
//     cerr << "F = " << F << endl;
  
  Vec3 Gk = G + k;
  Vec3 Gpk = Gp + k;
  Vec3 FGpk (0.0, 0.0, 0.0);
  FGpk[0] = F(0,0)*Gpk[0] + F(0,1)*Gpk[1] + F(0,2)*Gpk[2];
  FGpk[1] = F(1,0)*Gpk[0] + F(1,1)*Gpk[1] + F(1,2)*Gpk[2];
  FGpk[2] = F(2,0)*Gpk[0] + F(2,1)*Gpk[1] + F(2,2)*Gpk[2];
  
  return 0.5*dot(Gk, FGpk) + Vval;
}


void kCache::GetVals(double k, double &a, double &bPerp, double &bPar,
		     double &V)
{
  int i=0;
  while ((i<Cache.size()) && (fabs(Cache[i].k-k) > 1.0e-12))
    i++;
  if (i<Cache.size()) { // It's already in the cache
    a     = Cache[i].a;
    bPerp = Cache[i].bPerp;
    bPar  = Cache[i].bPar;
    V     = Cache[i].V;
  }
  else {
    a     = kPH.a(k);
    bPerp = kPH.bPerp(k);
    bPar  = kPH.bPar (k);
    V     = kPH.Vk    (k);
    kCachePoint newPoint;
    newPoint.a = a;
    newPoint.bPerp = bPerp;
    newPoint.bPar  = bPar;
    newPoint.V     = V;
    newPoint.k = k;
    Cache.push_back(newPoint);
  }
}
