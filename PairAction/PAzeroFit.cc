#include "PAFit.h"


// /// The following routines are used only if we are creating fits, not
// /// using them.
// #ifdef MAKE_FIT
// void PAclassicalFitClass::ReadParams(IOSectionClass &inSection)
// {
//     UsePBC = inSection.ReadVar ("Box", Box);
// }

// void PAclassicalFitClass::WriteBetaIndependentInfo (IOSectionClass &outSection)
// { }


// void PAclassicalFitClass::AddFit (Rho &rho)
// {
// }


// void PAclassicalFitClass::Error(Rho &rho, double &Uerror, double &dUerror)
// {
//   Uerror = 0.0;
//   dUerror = 0.0;
// }


// void PAclassicalFitClass::WriteFits (IOSectionClass &outSection)
// {
// }
// #endif


double PAzeroFitClass::U(double q, double z, double s2, int level)
{
  //  double beta = SmallestBeta;
  //  for (int i=0; i<level; i++)
  //    beta *= 2.0;

  //  double r = q+0.5*z;
  //  double rp = q-0.5*z;
  //  return (0.5*beta*(Potential->V(r)+Potential->V(rp)));
  return 0.0;
}

double PAzeroFitClass::dU(double q, double z, double s2, int level)
{
  //  double beta = SmallestBeta;
  //  for (int i=0; i<level; i++)
  //    beta *= 2.0;

  //  double r = q+0.5*z;
  //  double rp = q-0.5*z;
  //  return (0.5*(Potential->V(r)+Potential->V(rp)));
  return 0.0;
}



bool PAzeroFitClass::Read (IOSectionClass &in,
				double smallestBeta, int NumBetas)
{
  //  SmallestBeta = smallestBeta;
  // Read Particles;
    assert(in.OpenSection("Fits"));
    assert(in.OpenSection("Particle1"));
    Particle1.Read(in);
    in.CloseSection();
    assert(in.OpenSection("Particle2"));
    Particle2.Read(in);
    in.CloseSection();
  //  lambda = Particle1.lambda + Particle2.lambda;
  //  assert (lambda == 0.0);

  // Read Potential;
  //  assert(in.OpenSection("Potential"));
  //  Potential = ReadPH(in);
  //  in.CloseSection();
  //  in.CloseSection();
  return true;
}



