#ifndef JOSEPHSON_CLASS_H
#define JOSEPHSON_CLASS_H

#include "ActionBase.h"
#include <Common/PairAction/PAFit.h>


/// The JosephsonClass is an action class which stores the part of
/// the potential pair action that is summed in real space.  If the
/// potential is short range, it contains the whole potential action.
/// This action, in general, contains off diagaonal contributions.
class JosephsonClass : public ActionBaseClass
{
protected:
  int TotalTime;
public:
  double alpha;
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  JosephsonClass (PathDataClass &pathData);
};

#endif
