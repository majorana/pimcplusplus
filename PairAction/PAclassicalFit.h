#ifndef PA_CLASSICAL_FIT_H
#define PA_CLASSICAL_FIT_H
#include "PAFitBase.h"
#include "../Splines/BicubicSpline.h"

class PAclassicalFitClass : public PairActionFitClass
{
private:
  double Vlong_k (double boxVol, double k, int level);
  double Vlong (double q, int level);
  // Real space cutoff parameter;
  double alpha;
public:
#ifdef MAKE_FIT
  void ReadParams  (IOSectionClass &inSection);
  void WriteBetaIndependentInfo (IOSectionClass &outSection);
  /// Returns weighter RMS error
  void Error (Rho &rho, double &Uerror, double &dUerror);
  void AddFit (Rho &rho);
  void WriteFits(IOSectionClass &outSection);
  void ReadInput (IOSectionClass &inSection);
#endif
  bool Read  (IOSectionClass &inSection, double lowestBeta,
	      int NumBetas);
  double U(double q, double z, double s2, int level);
  double dU(double q, double z, double s2, int level);

  /////////////////////////
  /// Long-ranged stuff ///
  /////////////////////////
  /// The diagonal action only -- used for long-range breakup
  double Udiag(double q, int level);
  /// The q-derivative of the above
  double Udiag_p(double q, int level);
  /// The q-derivative of the above
  double Udiag_pp(double q, int level);
  /// The beta-derivative of the diagonal action
  double dUdiag    (double q, int level);
  /// The q-derivative of the above
  double dUdiag_p  (double q, int level);
  /// The q-derivative of the above
  double dUdiag_pp (double q, int level);
  /// The potential to which this action corresponds.
  double V  (double r);
  /// The q-derivative of the above
  double Vp (double r);
  /// The q-derivative of the above
  double Vpp(double r);

  bool IsLongRange();
  //  void DoBreakup(const dVec &box, const Array<dVec,1> &kVecs);

  PAclassicalFitClass()
  { 

  }
};

#endif
