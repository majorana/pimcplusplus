#ifndef DAVID_LONG_RANGE_CLASS_H
#define DAVID_LONG_RANGE_CLASS_H

#include "ActionBase.h"
#include "../Common/PairAction/PAFit.h"




/// The LongRangeClass is an action class responsible for the long
/// wavelength components of the action that are summed in k-space.
/// This class reads in and uses an optimized breakup from a file that
/// David supplies. 

class DavidLongRangeClass : public ActionBaseClass
{
protected:
  //  Array<PairActionFitClass*,2> &PairMatrix;
  //  Array<PairActionFitClass*,1> &PairArray;
  //  LinearGrid LongGrid;


  /// This calculates the quantity 
  /// \f$ X_k \equiv -\frac{4\pi}{\Omega k} \int_{r_c}^\infty dr \, r \sin(kr) V(r).\f$
  //  double CalcXk (int paIndex, int level, double k, double rc, JobType type);
  // This must be called after all of the OptimizedBreakup_x's

  //  int Level, ki;

  ///The values of the action will be uk*\rho_k*\rho_{-k}
  Array<double,1> uk;
  Array<double,1> duk;
public:
  // void Init(IOSectionClass &in);
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2,  int level);
  DavidLongRangeClass (PathDataClass &pathData);

};

#endif
