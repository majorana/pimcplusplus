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
  void DoFits();
  /** These arrays hold the coefficients of the expansion of the tail
      of the diagonal part of the action/potential in inverse powers
      of q or r.  The first index is for the inverse power, and the
      second gives the level.  That is, outside rcut, 
      \f[ U_\text{diag}(q, \text{level}) \approx Ucoef(0,level)/q +
      Ucoef(1,level)/r^2 + Ucoef(2,leve)/r^3 \f] */
  Array<double,2> Ucoefs, dUcoefs;
  Array<double,1> Vcoefs;
  /// The cutoff radius for the long range breakup.
  double rcut;  
public:
  Grid *qgrid, *tgrid;
  Array<BicubicSpline,1> Usplines, dUsplines;
#ifdef MAKE_FIT
  void ReadParams  (IOSectionClass &inSection);
  void WriteBetaIndependentInfo (IOSectionClass &outSection);
  /// Returns weighter RMS error
  void Error (Rho &rho, double &Uerror, double &dUerror);
  void DoFit (Rho &rho);
  void WriteFit(IOSectionClass &outSection);
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
  /// Set the cutoff radius for the long range breakup
  void Setrc(double rc);
  double Xk_U  (double k, int level);
  double Xk_dU (double k, int level);
  double Xk_V  (double k);
  double Vk    (double k);

  /////////////////////////
  /// Long-ranged stuff ///
  /////////////////////////
  bool IsLongRange();
  ///  void DoBreakup(const Vec3 &box, const Array<Vec3,1> &kVecs);

  PAcoulombBCFitClass()
  { 
    GridsAreMine = false; 
    NumBetas=0;
  }
  ~PAcoulombBCFitClass()
  { if (GridsAreMine){ delete qgrid; delete tgrid; } }
};

#endif
