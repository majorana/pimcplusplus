#ifndef ACTIONS_CLASS_H
#define ACTIONS_CLASS_H
#include "ShortRangeClass.h"
#include "ShortRangeApproximateClass.h"
#include "LongRangeClass.h"
#include "LongRangeRPAClass.h"
#include "ShortRangePotClass.h"
#include "LongRangePotClass.h"
#include "KineticClass.h"
#include "NodalActionClass.h"
#include "DavidLongRangeClass.h"
#include "TIP5PWaterClass.h"
#include "ST2WaterClass.h"
#include "OpenLoopImportance.h"
#include "StructureReject.h"

/// ActionsClass is a shell of a class holding all of the necessary
/// ActionBaseClass derivatives representing the different actions.
/// It includes kinetic, short range, long range, long range RPA
/// version, and nodal actions.  It may later contain importance
/// sampling "actions" as well.
class ActionsClass
{
private:
  Array<PairActionFitClass*,1> PairArray;
  Array<PairActionFitClass*,2> PairMatrix;
  PathDataClass &PathData;
  int MaxLevels; //is this the right place for this?
public:

  // Actions
  OpenLoopImportanceClass OpenLoopImportance;
  StructureRejectClass StructureReject;
  /// The Kinetic action
  KineticClass Kinetic;

  /// The short range part of the pair action.  This is the complete
  /// pair action in the case of short-range potententials.  The
  /// short-range action is summed in real space. 
  ShortRangeClass ShortRange;
  ShortRangeApproximateClass ShortRangeApproximate;

  /// The long range part of the action, which is summed in k-space.  
  LongRangeClass LongRange;

  /// The Random Phase Approximation-corrected form of the above.
  LongRangeRPAClass LongRangeRPA;

  ///David's Long Range Class
  DavidLongRangeClass DavidLongRange;

  /// Action for simulations using the TIP5P water model
  TIP5PWaterClass TIP5PWater;

  /// Action for simulations using the ST2 water model
  ST2WaterClass ST2Water;

  /// This array of actions are used for Restricted PIMC for
  /// fermions.  These effective actions ensure that the paths do not
  /// cross the nodes of some trial density matrix with respective to
  /// the reference slice.
  Array<NodalActionClass *,1> NodalActions;
  //DiagonalClass Diagonal;
  //ImportanceSampleClass ImportanceSample;

  // Potentials
  ShortRangePotClass ShortRangePot;
  LongRangePotClass  LongRangePot;
  
  /// Stores whether we use Random Phase Approximation corrections to
  /// the long range action.
  bool UseRPA;

  /// Stores number of images to sum over for kinetic action and energy.
  int NumImages;


  /// Return the all the energies for this processor's segment of
  /// the path.  Must do global sum to get total energy.
  void Energy (double& kinetic, double &duShort, double &duLong, 
	       double &node, double &vShort, double &vLong);

  /// Read the action parameters from the input file and do the
  /// necessary initialization.  This reads the pair actions, and does
  /// long range breakups and RPA corrections if necessary.
  void Read(IOSectionClass &in);
  ActionsClass(PathDataClass &pathData) : 
    ShortRange(pathData,PairMatrix),
    ShortRangeApproximate(pathData,PairMatrix),
    ShortRangePot(pathData, PairMatrix),
    LongRange(pathData,PairMatrix,PairArray), 
    DavidLongRange(pathData),
    LongRangeRPA(pathData, PairMatrix, PairArray),
    LongRangePot(pathData, PairMatrix),
    OpenLoopImportance(pathData),
    Kinetic(pathData),
    PathData(pathData),
    TIP5PWater(pathData),
    ST2Water(pathData),
    StructureReject(pathData),
    NumImages(1)
  {
    ///Do nothing for now
  }





};

#endif
