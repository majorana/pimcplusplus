#ifndef KINETIC_CLASS_H
#define KINETIC_CLASS_H

#include "ActionBase.h"

class KineticClass : public ActionBaseClass
{
public:
  void Read (IOSectionClass &in);
  double Action (int slice1, int slice2, 
		 const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  KineticClass (PathDataClass &pathData);
};


#endif
