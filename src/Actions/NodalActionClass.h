#ifndef NODAL_ACTION_CLASS_H
#define NODAL_ACTION_CLASS_H

#include "ActionBase.h"
#include "../Common/Splines/CubicSpline.h"

class PathClass;

class NodalActionClass : public ActionBaseClass
{
public:
  virtual bool IsPositive (int slice) = 0;
  NodalActionClass (PathDataClass &pathData) :
    ActionBaseClass (pathData)
  {

  }
};

/// FreeNodalActionClass implements the nodal action corresponding to
/// the free fermion density matrix in periodic boundary conditions.
/// Currently, the beta-derivative of the action is computed only
/// approximately, but in a way that becomes accurate as \f$ \tau
/// \rightarrow 0 \f$.
class FreeNodalActionClass : public NodalActionClass
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
  double ActionkSum (double L, double lambdaTau, double disp);

  Array<double,2> DetMatrix, Cofactors;
  Array<dVec,1> GradVec, SavePath;
  double Det(int slice);
  void GradientDet (int slice, double &det, Array<dVec,1> &gradient);
  void GradientDetFD (int slice, double &det, Array<dVec,1> &gradient);
  double NodalDist (int slice);
  double HybridDist(int slice, double lambdaTau);
  /// This returns the upper bound on the distance to a node by
  /// returning the minimum distance to particle coincidence.
  double MaxDist(int slice);
  /// This calculates the distance to the node along the line in the
  /// direction of the gradient by a bisection search
  double LineSearchDist (int slice);
  /// This calculates the nodal distance by an iterated Newton-Raphson
  /// approach
  double NewtonRaphsonDist (int slice);
  int SpeciesNum;
  int NumGradDists, NumLineDists;
public:
  double Action (int slice1, int slice2,
		 const Array<int,1> &activeParticles,
		 int level);
  double d_dBeta (int slice1, int slice2, 
		  int level);
  double SimpleAction (int slice1, int slice2,
		       const Array<int,1> &activeParticles,
		       int level);
  /// Returns true if the nodal restriction is satisfied for my
  /// species at timeslice slice.  If slice is the reference slice,
  /// returns true.
  bool IsPositive (int slice);
  void Read (IOSectionClass &in);
  FreeNodalActionClass (PathDataClass &pathData, int speciesNum);
};


#endif
