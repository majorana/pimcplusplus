#ifndef PA_TRICUBIC_FIT_H
#define PA_TRICUBIC_FIT_H
#include "PAFitBase.h"
#include "../Splines/MyTricubicSpline.h"
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
  Array<MyTricubicSpline,1> Usplines, dUsplines;
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

  bool IsLongRange();

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
