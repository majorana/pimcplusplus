#ifndef PHONON_FREE_ENERGY_H
#define PHONON_FREE_ENERGY_H

#include "../Blitz.h"

class PhononFreeEnergy
{
  int Vorder, Torder;
  Array<double,1> FCoefs, UCoefs;
  void EvalBasis (double V, double T,
		  Array<double,1> basis);
  Array<double,1> Btmp;
  const double kB;  // Hartrees per Kelvin
public:
  double F_VT(double V, double T);
  double F_PT(double P, double T);
  double F_PV(double P, double V);

  double P_VT(double V, double T);
  double P_VT_FD(double V, double T);
  double dP_dV (double V, double T);
  double dP_dV_FD (double V, double T);
  double dP_dT (double V, double T);
  double dP_dT_FD (double V, double T);
  
  void FitF (int vorder, int vorder,
	     Array<double,1> &F,
	     Array<double,1> &V,
	     Array<double,1> &T);
  void FitU (int vorder, int vorder,
	     Array<double,1> &U,
	     Array<double,1> &V,
	     Array<double,1> &T);
  PhononFreeEnergy() : kB(3.16681526543384e-06)
  {
  }
};

#endif
