#ifndef PA_COULOMBBC_FIT_H
#define PA_COULOMBBC_FIT_H
#include "PAFitBase.h"
#include "../Splines/BicubicSpline.h"

class PAcoulombBCFitClass : public PairActionFitClass
{
private:
  bool GridsAreMine;
  double Vlong_k (double boxVol, double k, int level);
  double dVlong_k (double boxVol, double k, int level);
  double Vlong (double q, int level);
  double dVlong (double q, int level);
  // Real space cutoff parameter;
  double alpha;
public:
  Grid *qgrid, *tgrid;
  Array<BicubicSpline,1> Usplines, dUsplines;
#ifdef MAKE_FIT
  void ReadParams  (IOSectionClass &inSection);
  void WriteBetaIndependentInfo (IOSectionClass &outSection);
  /// Returns weighter RMS error
  void Error (Rho &rho, double &Uerror, double &dUerror);
  void AddFit (Rho &rho);
  void WriteFits(IOSectionClass &outSection);
#endif
  void Write (IOSectionClass &outSection);
  bool Read  (IOSectionClass &inSection, double lowestBeta,
	      int NumBetas);
  double U(double q, double z, double s2, int level);
  double dU(double q, double z, double s2, int level);
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

  /////////////////////////
  /// Long-ranged stuff ///
  /////////////////////////
  bool IsLongRange();
  void DoBreakup(const dVec &box, const Array<dVec,1> &kVecs);

  PAcoulombBCFitClass()
  { 
    GridsAreMine = false; 
    NumBetas=0;
  }
  ~PAcoulombBCFitClass()
  { if (GridsAreMine){ delete qgrid; delete tgrid; } }
};

#endif
