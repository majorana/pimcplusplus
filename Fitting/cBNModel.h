#ifndef CBN_MODEL_H
#define CBN_MODEL_H

#include "VinetEOS2.h"
#include "NonlinearFit.h"
#include <vector>
#include "PhononFreeEnergy.h"
#include "DebyeModel.h"
#include "RamanFreq.h"


class cBNModel
{
public:
  VinetEOSClass Static;
  SplineFreeEnergy Phonon;
  inline double P(double V, double T)
  { return Static.Pressure(V) + Phonon.P(V,T); }
  double FindV(double P, double T);

  void SetStatic (string fname);
  void SetPhonon (string fname);
};

#endif
