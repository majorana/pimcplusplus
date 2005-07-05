#include "PAMonopoleFit.h"

/// The following routines are used only if we are creating fits, not
/// using them.

void PACoulombFitClass::ReadParams(IOSectionClass &inSection)
{
}

void PACoulombFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
{ }


void PACoulombFitClass::DoFit (Rho &rho)
{
}


void PACoulombFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
{
  Uerror = 0.0;
  dUerror = 0.0;
}


void PACoulombFitClass::WriteFit (IOSectionClass &outSection)
{
}



double PACoulombFitClass::U(double q, double z, double s2, int level)
{
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  double r = q + 0.5*z;
  double rp = q - 0.5*z;
  double V = 0.5*Z1Z2*(1.0/r + 1.0/rp);
  return (beta*V);
}

double PACoulombFitClass::dU(double q, double z, double s2, int level)
{
  double r = q + 0.5*z;
  double rp = q - 0.5*z;
  double V = 0.5*Z1Z2*(1.0/r + 1.0/rp);
  return (V);
}


bool PACoulombFitClass::Read (IOSectionClass &in,
			      double smallestBeta, int numBetas)
{
  SmallestBeta = smallestBeta;
  NumBetas=numBetas;
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
  in.OpenSection("Potential");
  if (!in.ReadVar ("Z1Z2", Z1Z2))
    Z1Z2 = 0.0;
  in.CloseSection();
  in.CloseSection();
  return true;
}



/////////////////////////
/// Long-ranged stuff ///
/////////////////////////
bool PACoulombFitClass::IsLongRange()
{
  return (Z1Z2 != 0.0);
}


/// The diagonal action only -- used for long-range breakup
double PACoulombFitClass::Udiag(double q, int level)
{  
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return beta*Z1Z2/q;
}

/// The q-derivative of the above
double PACoulombFitClass::Udiag_p(double q, int level) 
{  
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return -beta*Z1Z2/(q*q);
}

/// The q-derivative of the above
double PACoulombFitClass::Udiag_pp(double q, int level) 
{  
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return 2.0*beta*Z1Z2/(q*q*q);
}

/// The beta-derivative of the diagonal action
double PACoulombFitClass::dUdiag    (double q, int level) 
{
  return Z1Z2/q;
}

/// The q-derivative of the above
double PACoulombFitClass::dUdiag_p  (double q, int level) 
{
  return -Z1Z2/(q*q);
}

/// The q-derivative of the above
double PACoulombFitClass::dUdiag_pp (double q, int level) 
{
  return 2.0*Z1Z2/(q*q*q);
}

/// The potential to which this action corresponds.
double PACoulombFitClass::V  (double r) 
{
  return Z1Z2/r;
}

/// The q-derivative of the above
double PACoulombFitClass::Vp (double r)
{
  return -Z1Z2/(r*r);
}

/// The q-derivative of the above
double PACoulombFitClass::Vpp(double r) 
{
  return 2.0*Z1Z2/(r*r*r);
}

void PACoulombFitClass::Setrc(double rc)
{
  rcut = rc;
}

double PACoulombFitClass::Xk_U(double k, int level)
{
  double beta = SmallestBeta;
  for (int i=0; i<level; i++)
    beta *= 2.0;
  return -beta*4.0*M_PI*Z1Z2/(k*k)*cos(k*rcut);
}

double PACoulombFitClass::Xk_dU(double k, int level)
{
  return -4.0*M_PI*Z1Z2/(k*k)*cos(k*rcut);
}

double PACoulombFitClass::Xk_V(double k)
{
  return -4.0*M_PI*Z1Z2/(k*k)*cos(k*rcut);
}

double PACoulombFitClass::Vk(double k)
{
  return 4.0*M_PI*Z1Z2/(k*k);
}
