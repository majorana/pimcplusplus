#include "PAFit.h"
#include "../Splines/BicubicSpline.h"

/// The following routines are used only if we are creating fits, not
/// using them.
#ifdef MAKE_FIT
void PAclassicalFitClass::ReadParams(IOSectionClass &inSection)
{
  UsePBC = inSection.ReadVar ("Box", Box);
}

void PAclassicalFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{ }


void PAclassicalFitClass::AddFit (Rho &rho)
{
}


void PAclassicalFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
{
  Uerror = 0.0;
  dUerror = 0.0;
}


void PAclassicalFitClass::WriteFits (IOSectionClass &outSection)
{
}
#endif


double PAclassicalFitClass::U(double q, double z, double s2, int level)
{
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;

  double r = q+0.5*z;
  double rp = q-0.5*z;
  return (0.5*beta*(Pot->V(r)+Pot->V(rp)));
}

double PAclassicalFitClass::dU(double q, double z, double s2, int level)
{
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;

  double r = q+0.5*z;
  double rp = q-0.5*z;
  return (0.5*(Pot->V(r)+Pot->V(rp)));
}



bool PAclassicalFitClass::Read (IOSectionClass &in,
				double smallestBeta, int NumBetas)
{
  SmallestBeta = smallestBeta;
  // Read Particles;
  assert(in.OpenSection("Fits"));
  assert(in.OpenSection("Particle1"));
  Particle1.Read(in);
  in.CloseSection();
  assert(in.OpenSection("Particle2"));
  Particle2.Read(in);
  in.CloseSection();
  lambda = Particle1.lambda + Particle2.lambda;
  assert (lambda == 0.0);

  // Read Potential;
  assert(in.OpenSection("Potential"));
  Pot = ReadPotential(in);
  in.CloseSection();
  in.CloseSection();
  return true;
}



/////////////////////////
/// Long-ranged stuff ///
/////////////////////////
bool PAclassicalFitClass::IsLongRange()
{
  // This needs to be fixed.  We need to add this kind of function to
  // the potential base class. 
  return true;
}

double PAclassicalFitClass::Vlong(double q, int level)
{
  if (q <= 0.0)
    return 2.0/sqrt(M_PI)*Z1Z2*alpha;
  else 
    return Z1Z2/q*erf(alpha*q);
}

double PAclassicalFitClass::Vlong_k(double boxVol, double k, int level)
{
  if (k <= 0.0)
    k = 1.0e-30;
  return 4.0*M_PI/(boxVol*k*k)*exp(-k*k/(4.0*alpha*alpha));
}


void PAclassicalFitClass::DoBreakup(const dVec& box,const Array<dVec,1> &kVecs)
{
  


}
