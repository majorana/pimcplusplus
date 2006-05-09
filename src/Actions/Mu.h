#ifndef MU_CLASS_H
#define MU_CLASS_H

#include "ActionBase.h"




/// The ShortRangeClass is an action class which stores the part of
/// the potential pair action that is summed in real space.  If the
/// potential is short range, it contains the whole potential action.
/// This action, in general, contains off diagaonal contributions.
class MuClass : public ActionBaseClass
{
protected:
public:
  void Read (IOSectionClass &in);
  bool PadWorm();
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta(int x, int y, int z)
  {
    return 0.0;
  }
  double Mu;
  MuClass (PathDataClass &pathData);
    
};

#endif
