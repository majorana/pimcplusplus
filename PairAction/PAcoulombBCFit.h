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
  // Product of the two charges.
  double Z1Z2;
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
