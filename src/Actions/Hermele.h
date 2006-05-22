#ifndef HERMELE_CLASS_H
#define HERMELE_CLASS_H

#include "ActionBase.h"
#include <Common/PairAction/PAFit.h>
#include <fftw3.h>

/// The JosephsonClass is an action class which stores the part of
/// the potential pair action that is summed in real space.  If the
/// potential is short range, it contains the whole potential action.
/// This action, in general, contains off diagaonal contributions.
class HermeleClass : public ActionBaseClass
{
protected:
  int TotalTime;
public:
  double alpha;
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  double T;
  double J;
  fftw_complex *inArray, *out;
  fftw_plan p;
  void Phi2Omega();
  double ComputeFourierAction();
  double c;
  double K_s;
  double omega_c2;

  HermeleClass (PathDataClass &pathData);
  ~HermeleClass()
  {
    fftw_destroy_plan(p);
    fftw_free(inArray);
    fftw_free(out);
  }

};

#endif
