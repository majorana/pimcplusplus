#ifndef STRUCTURE_REJECT_H
#define STRUCTURE_REJECT_H

#include "ActionBase.h"

class StructureRejectClass : public ActionBaseClass
{
  Array<double,1> Sk;
  int Species1;
  int Species2;
  int TotalCounts;
public:
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  StructureRejectClass (PathDataClass &pathData);
};


#endif
