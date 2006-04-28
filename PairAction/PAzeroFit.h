#ifndef PA_ZERO_FIT_H
#define PA_ZERO_FIT_H
#include "PAFitBase.h"


class PAzeroFitClass : public PairActionFitClass
{
public:

  void ReadParams  (IOSectionClass &inSection) {};
  void WriteBetaIndependentInfo (IOSectionClass &outSection) {};
  /// Returns weighter RMS error
  void Error (Rho &rho, double &Uerror, double &dUerror) {};
  void DoFit (Rho &rho) {};
  void WriteFit(IOSectionClass &outSection) {};
  void ReadInput (IOSectionClass &inSection) {};
  
  bool Read  (IOSectionClass &inSection, double lowestBeta,
	      int NumBetas);
  double U(double q, double z, double s2, int level);
  double dU(double q, double z, double s2, int level);
  void Derivs (double q, double z, double s2, int level,
	       double &d_dq, double &d_dz, double &d_ds);
  double V(double r);
  bool IsLongRange();
  PAzeroFitClass()
  { 

  }
};

#endif
