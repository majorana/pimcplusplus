#include "kSpacePH.h"
#include "../Integration/GKIntegration.h"

class aIntegrand
{
protected:
  Potential &PH;
  double k, kinv;
public:
  inline void Setk (double newk) 
  {    k = newk;    kinv = 1.0/newk; }
   
  inline double operator()(double r)
  { return 4.0*M_PI*(PH.A(r)-1.0)*r*sin(k*r)*kinv; }

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
    double j0 = sin(kr)/kr;
    double j2 = (3.0*krInv*krInv*krInv-krInv)*sin(kr) - 
      3.0*krInv*krInv*cos(kr);
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
    double j0 = sin(kr)/kr;
    double j2 = (3.0*krInv*krInv*krInv-krInv)*sin(kr) - 
      3.0*krInv*krInv*cos(kr);
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

