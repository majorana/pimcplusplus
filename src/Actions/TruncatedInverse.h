#ifndef TRUNCATED_INVERSE__CLASS_H
#define TRUNCATED_INVERSE__CLASS_H

#include "../MirroredClass.h"
#include "NodalActionClass.h"
#include <Common/Splines/CubicSpline.h>


/// VartionalPIClass implements the nodal action corresponding to
/// the free fermion density matrix in periodic boundary conditions.
/// Currently, the beta-derivative of the action is computed only
/// approximately, but in a way that becomes accurate as \f$ \tau
/// \rightarrow 0 \f$.
class TruncatedInverseClass : public NodalActionClass
{
private:
  PathClass &Path;
  void calc_u();
  Array<double,1> u;
  Array<double,1> newCol;
  /// These splines will hold the free-particle action for
  /// periodic boundary conditions.  The array is over time-slice
  /// separation from the reference slice.
  Array<double,2> DetMatrix;
  Array<double,2> SmallDetMatrix;
  Mirrored1DClass<double> DeterminantList;

public:
  double SingleAction (int slice1, int slice2,
		       const Array<int,1> &activeParticles,
		       int level);
  double d_dBeta (int slice1, int slice2, 
		  int level);
  void BuildDeterminantMatrix();
  void BuildSmallDeterminantMatrix();
  void CheckDeterminantMatrix();
  void Read (IOSectionClass &in);
  bool IsGroundState();
  bool IsPositive(int x);
  NodeType Type();
  void AcceptCopy(int slice1, int slice2);
  void RejectCopy(int slice1, int slice2);
  void WriteInfo (IOSectionClass &out);
  int ChangedColumn;
  dVec olddvec;
  dVec newdvec;
  
  TruncatedInverseClass (PathDataClass &pathData);
};

#endif

