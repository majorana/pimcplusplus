#ifndef PA_S_FIT_H
#define PA_S_FIT_H
#include "PAFitBase.h"
#include "../Splines/BicubicSpline.h"

class PAsFitClass : public PairActionFitClass
{
private:
  bool GridsAreMine;
  Array<double,1> Coefs;
public:
  Grid *qgrid, *ygrid;
  Array<MultiBicubicSpline,1> Usplines, dUsplines;
  Array<double,1> sMax;
  int Order;
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
  PAsFitClass()
  { 
    GridsAreMine = false; 
    NumBetas=0;
  }
  ~PAsFitClass()
  { if (GridsAreMine){ delete qgrid; delete ygrid; } }
};

#endif
