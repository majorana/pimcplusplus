#ifndef QBOX_DFT_ACTION_CLASS_H
#define QBOX_DFT_ACTION_CLASS_H

#include "ActionBase.h"
#include <vector>
#include <fstream>

#ifdef USE_QBOX

#include "/home/gergely/build/qbox-1.44.0/src/qbInterface.h"
class QBoxActionClass: public ActionBaseClass
{
	Array<string, 1> ptclSet;
  // interface to qbox
  qbox::qbLib* qb;
	int age;
	double prevEnergy, newEnergy, oldEnergy;
	int steps;
	double tol, w;
	bool feedback;
  string command;
  int toread, towrite;
	//ofstream toqbox;
	//ifstream fromqbox;
	Array<string,1> ptclID;
  
	public:
  // hack
  ofstream out;
  bool isAction;

	double SingleAction(int slice1,int slice2,
		      const Array<int,1> &activeParticles,int level);
	double ComputeEnergy(int slice1,int slice2,
					const Array<int,1> &activeParticles,int level);
  double d_dBeta (int slice1, int slice2, int level);
  std::string GetName();
	double Collect();
	void SetPtclPos(int slice, Array<int,1> activeParticles);
  void Read (IOSectionClass &in);
  void AcceptCopy (int slice1, int slice2);
  void RejectCopy (int slice1, int slice2);
  
  QBoxActionClass(PathDataClass &pathData);
};
#endif

#endif
