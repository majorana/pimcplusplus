#ifndef ACTION_BASE_H
#define ACTION_BASE_H

#include "../Common.h"
#include <Common/IO/IO.h>

using namespace IO;

class PathDataClass;
class PathClass;

/// The ActionBaseClass is an abstract base class from which all the
/// physical action classes are derived.  It's primary methods are
/// Action and it's beta-derivative.
class ActionBaseClass
{
protected:
  PathDataClass &PathData;
  PathClass &Path;
public:
  /// This function takes a range of time slices from slice1 to
  /// slice2, inclusive, and an array of particles which are changing
  /// positions.  The level is used in bisection moves in which, at
  /// higher levels, we skip 2^level intervening time slices, building
  /// up a new path in a recursive style.
  virtual double SingleAction(int slice1, int slice2,
			      const Array<int,1> &activeParticles,
			      int level) = 0;

  /// This is the main function that gets called.  It simply returns
  /// SingleAction if we are not using correlated sampling.  If we are
  /// using correlated sampling, it returns the average of the actions
  /// for the A and B configurations.  At level 0, we must do
  /// something special, so Action should not be called with
  /// Correlated sampling on and level = 0;  Note that this is
  /// virtual, so that the nodal classes can override.
  virtual double Action(int slice1, int slice2,
			const Array<int,1> &activeParticles,
			int level);

  /// This function returns the \f$\beta\f$-derivative of the above
  /// function.  Since we are interested in total energy, it does not
  /// take a list of particles we are moving. 
  virtual double d_dBeta (int slice1, int slice2,
			  int level) = 0;
  // This returns the sum over all time slices, using MPI
  // to sum over processors if necessary.
  virtual void Read (IOSectionClass &in);
  virtual void ShiftData (int slices2Shift);
  ActionBaseClass(PathDataClass &pathData);				   
};

/// The PotentialBaseClass is an abstract base class for storing
/// potentials that go into the energy.  Currently, the ShortRange and
/// LongRange potentials are derived from it.
class PotentialBaseClass
{
protected:
  PathDataClass &PathData;
  PathClass &Path;
public:
  virtual double V (int slice) = 0;
  PotentialBaseClass (PathDataClass &pathData);
};

#endif
