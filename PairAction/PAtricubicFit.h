#ifndef PA_TRICUBIC_FIT_H
#define PA_TRICUBIC_FIT_H
#include "PAFitBase.h"
#include "../Splines/TricubicSpline.h"
#ifdef MAKE_FIT
#include "../MPI/Communication.h"
#endif

class PAtricubicFitClass : public PairActionFitClass
{
private:
  bool GridsAreMine;
#ifdef MAKE_FIT
  CommunicatorClass Comm;
#endif
public:
  Grid *qgrid, *ygrid, *tgrid;
  Array<TricubicSpline,1> Usplines, dUsplines;
  Array<double,1> sMax, sMaxInv;
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
  PAtricubicFitClass()
  { 
#ifdef MAKE_FIT
    Comm.SetWorld();
#endif
    GridsAreMine = false; 
    NumBetas=0;
  }
  ~PAtricubicFitClass()
  { if (GridsAreMine){ delete qgrid; delete ygrid; } }
};

#endif
