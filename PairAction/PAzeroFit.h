#ifndef PA_ZERO_FIT_H
#define PA_ZERO_FIT_H
#include "PAFitBase.h"


class PAzeroFitClass : public PairActionFitClass
{
public:
  #ifdef MAKE_FIT
    void ReadParams  (IOSectionClass &inSection) {};
    void WriteBetaIndependentInfo (IOSectionClass &outSection) {};
    /// Returns weighter RMS error
    void Error (Rho &rho, double &Uerror, double &dUerror) {};
    void AddFit (Rho &rho) {};
    void WriteFits(IOSectionClass &outSection) {};
    void ReadInput (IOSectionClass &inSection) {};
  #endif
  bool Read  (IOSectionClass &inSection, double lowestBeta,
	      int NumBetas);
  double U(double q, double z, double s2, int level);
  double dU(double q, double z, double s2, int level);
  bool IsLongRange();
  PAzeroFitClass()
  { 

  }
};

#endif
