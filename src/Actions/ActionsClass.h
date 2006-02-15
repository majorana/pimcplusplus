#ifndef ACTIONS_CLASS_H
#define ACTIONS_CLASS_H
#include "ShortRangeClass.h"
#include "ShortRangeOnClass.h"
#include "ShortRangeApproximateClass.h"
#include "LongRangeClass.h"
#include "LongRangeRPAClass.h"
#include "ShortRangePotClass.h"
#include "LongRangePotClass.h"
#include "KineticClass.h"
#include "NodalActionClass.h"
#include "FreeNodalActionClass.h"
#include "GroundStateNodalActionClass.h"
#include "FixedPhaseActionClass.h"
#include "DavidLongRangeClass.h"
#include "TIP5PWaterClass.h"
#include "ST2WaterClass.h"
#include "OpenLoopImportance.h"
#include "StructureReject.h"
#include "KineticSphereClass.h"

/// ActionsClass is a shell of a class holding all of the necessary
/// ActionBaseClass derivatives representing the different actions.
/// It includes kinetic, short range, long range, long range RPA
/// version, and nodal actions.  It may later contain importance
/// sampling "actions" as well.
class ActionsClass
{
private:
  // This stores pointers to the pair actions.
  Array<PairActionFitClass*,1> PairArray;
  PathDataClass &PathData;
  int MaxLevels; //is this the right place for this?
  void ReadNodalActions (IOSectionClass &in);
public:
  // This stores pointers to pair action fits for each pair of species.
  Array<PairActionFitClass*,2> PairMatrix;

  /// Used to keep track of the total action
  double TotalA, TotalB;

  /// Specifies whether to use long range breakups
  bool UseLongRange;


  // Actions
  OpenLoopImportanceClass OpenLoopImportance;
  StructureRejectClass StructureReject;
  /// The Kinetic action
  KineticClass Kinetic;
  KineticSphereClass KineticSphere;

  /// The short range part of the pair action.  This is the complete
  /// pair action in the case of short-range potententials.  The
  /// short-range action is summed in real space. 
  ShortRangeClass ShortRange;
  ShortRangeOnClass ShortRangeOn;
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

  Potential &GetPotential (int species1, int species2);
  /// Return the all the actions for this processor's segment of
  /// the path.  Must do global sum to get total action.
  void GetActions (double& kinetic, double &duShort, double &duLong, 
		   double &node);

  /// This function adds to F the current forces calculated from the
  /// gradient of the action.  Note that you must do an AllSum over
  /// the clone processors to get the total.
  void GetForces(const Array<int,1> &ptcls, 
		 Array<dVec,1> &Fshort, Array<dVec,1> &Flong);
  /// Finite difference version for testing.
  void GetForcesFD(const Array<int,1> &ptcls, Array<dVec,1> &F);

  /// Return the all the energies for this processor's segment of
  /// the path.  Must do global sum to get total energy.
  void Energy (double& kinetic, double &duShort, double &duLong, 
	       double &node, double &vShort, double &vLong);

  /// Read the action parameters from the input file and do the
  /// necessary initialization.  This reads the pair actions, and does
  /// long range breakups and RPA corrections if necessary.
  void Read(IOSectionClass &in);

  /// This routine does any necessary shifting of stored data for the
  /// paths.  It just calls the ShiftData
  void ShiftData (int slicesToShift);
  void MoveJoin (int oldJoinPos, int newJoinPos);

  void AcceptCopy (int startSlice, int endSlice,
		   const Array<int,1> &activeParticles);
  void RejectCopy (int startSlice, int endSlice,
		   const Array<int,1> &activeParticles);
  /// This should be called after the paths have been constructed to
  /// initialize any cached data;
  void Init();

  void Setk(Vec3 k);

  bool HaveLongRange();

  /// This function writes any pertinent information related to the
  /// actions to the output file.
  void WriteInfo (IOSectionClass &out);

  void UpdateNodalActions();

  inline int GetMaxLevels() { return MaxLevels; }

  ActionsClass(PathDataClass &pathData) : 
    ShortRange(pathData,PairMatrix),
    ShortRangeOn(pathData,PairMatrix),
    ShortRangeApproximate(pathData,PairMatrix),
    ShortRangePot(pathData, PairMatrix),
    LongRange(pathData,PairMatrix,PairArray), 
    DavidLongRange(pathData),
    LongRangeRPA(pathData, PairMatrix, PairArray),
    LongRangePot(pathData, PairMatrix),
    OpenLoopImportance(pathData),
    Kinetic(pathData),
    KineticSphere(pathData),
    PathData(pathData),
    TIP5PWater(pathData),
    ST2Water(pathData),
    StructureReject(pathData),
    NumImages(1),
    UseLongRange(true)
  {
    ///Do nothing for now
  }





};

#endif
