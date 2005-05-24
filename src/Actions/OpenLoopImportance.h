#ifndef OPENLOOP_IMPORTANCE_H
#define OPENLOOP_IMPORTANCE_H

#include "ActionBase.h"

typedef enum {NOIMP,DISTIMP,DISPXIMP,CONSTSHIFT} SampleChoice;

class OpenLoopImportanceClass : public ActionBaseClass
{
public:
  void Read (IOSectionClass &in);
  double Action (int slice1, int slice2, 
		 const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  OpenLoopImportanceClass (PathDataClass &pathData);
  SampleChoice ImpChoice;
  double Shift;
};


#endif
