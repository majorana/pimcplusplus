#ifndef PA_COULOMBBC_FIT_H
#define PA_COULOMBBC_FIT_H

#include "PAFitBase.h"

class PAcoulombBCFitClass : public PairActionFitClass
{
private:
  bool GridsAreMine;
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
  PAcoulombBCFitClass()
  { 
    GridsAreMine = false; 
    NumBetas=0;
  }
  ~PAcoulombBCFitClass()
  { if (GridsAreMine){ delete qgrid; delete tgrid; } }
};

#endif
