#ifndef PA_SZ_FIT_H
#define PA_SZ_FIT_H

#include "PAFitBase.h"

class PAszFitClass : public PairActionFitClass
{
private:
  bool GridIsMine;
  Array<double, 1> Coefs;
public:
  int Order;
  Array<double,1> smax;  // Array index is level
  Grid *qgrid;
  Array<MultiCubicSpline,1> Ukj;
  Array<MultiCubicSpline,1> dUkj;
#ifdef MAKE_FIT
  void ReadParams  (IOSectionClass &inSection);
  void WriteBetaIndependentInfo (IOSectionClass &outSection);
  void AddFit (Rho &rho);
  void WriteFits(IOSectionClass &outSection);
  void Error(Rho &rho, double &Uerror, double &dUerror);
#endif
  void Write (IOSectionClass &outSection);
  bool Read  (IOSectionClass &inSection, double lowestBeta,
	      int NumBetas);
  double U(double r, double rp, double costheta, int level);
  double dU(double r, double rp, double costheta, int level);
  PAszFitClass()
  { 
    GridIsMine = false; 
    NumBetas=0;
  }
  ~PAszFitClass()
  { if (GridIsMine) delete qgrid;  }
};

#endif
