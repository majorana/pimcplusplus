#ifndef PWSCF_DFT_ACTION_CLASS_H
#define PWSCF_DFT_ACTION_CLASS_H

#include "ActionBase.h"
#include <vector>
#include <fstream>
#include "popen2.h"

class pwscfActionClass: public ActionBaseClass
{
  int qCount;
	public:
  // hack
  ofstream out;

	double SingleAction(int slice1,int slice2,
		      const Array<int,1> &activeParticles,int level);
	double ComputeEnergy(int slice1,int slice2,
					const Array<int,1> &activeParticles,int level);
  double d_dBeta (int slice1, int slice2, int level);
  std::string GetName();
	void SetPtclPos(int slice, const char* filename);
  void Read (IOSectionClass &in);
  
  pwscfActionClass(PathDataClass &pathData);
};

#endif
