#ifndef OPENLOOP_IMPORTANCE_H
#define OPENLOOP_IMPORTANCE_H

#include "ActionBase.h"

typedef enum {NOIMP,DISTIMP,DISPXIMP,CONSTSHIFT,TUNEDFUNCTION,RETUNEDFUNCTION, POLYNOMIAL,EXPONENTIAL} SampleChoice;

typedef enum {MIN_IMAGE_DISP,MIN_IMAGE_DIST,DISP,DIST} Xvalue;

class OpenLoopImportanceClass : public ActionBaseClass
{
public:
  Array<double,1> Polynom;
  double a;
  double alpha;
  double s;
  void Read (IOSectionClass &in);
  double SingleAction (int slice1, int slice2, 
		       const Array<int,1> &activeParticles, int level);
  double d_dBeta (int slice1, int slice2, int level);
  OpenLoopImportanceClass (PathDataClass &pathData);
  SampleChoice ImpChoice;
  Xvalue Xis;
  double Shift;
};


#endif
