#ifndef NODAL_ACTION_CLASS_H
#define NODAL_ACTION_CLASS_H

#include "ActionBase.h"
#include "../Common/Splines/CubicSpline.h"

class PathClass;


class FreeNodalActionClass : public ActionBaseClass
{
private:
  PathClass &Path;

  /// These splines will hold the free-particle action for
  /// periodic boundary conditions.  The array is over time-slice
  /// separation from the reference slice.
  Array<TinyVector<CubicSpline,NDIM>,1> ActionSplines;
  TinyVector<LinearGrid,NDIM> ActionGrids;
  void SetupFreeActions();
  double ActionImageSum (double L, double lambdaTau, double disp);

  Array<double,2> DetMatrix, Cofactors;
  Array<dVec,1> GradVec;
  double Det();
  void GradientDet (int slice, double &det, Array<dVec,1> &gradient);
  void GradientDetFD (int slice, double &det, Array<dVec,1> &gradient);
  double NodalDist (int slice);
  int SpeciesNum;
public:
  double Action (int slice1, int slice2,
		 const Array<int,1> &activeParticles,
		 int level);
  double d_dBeta (int slice1, int slice2, 
		  int level);
  void Read (IOSectionClass &in);
  FreeNodalActionClass (PathDataClass &pathData, int speciesNum);
};


#endif
