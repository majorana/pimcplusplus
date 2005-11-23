#include "ActionBase.h"
#include "../PathDataClass.h"

void 
ActionBaseClass::Read (IOSectionClass &in)
{
  // Do nothing for now
}

ActionBaseClass::ActionBaseClass(PathDataClass &pathData) : 
  PathData(pathData), Path(pathData.Path)
{
  /* Do nothing */
}

void 
ActionBaseClass::ShiftData (int slicesToShift)
{
  // Do nothing 
}

PotentialBaseClass::PotentialBaseClass(PathDataClass &pathData) : 
  PathData(pathData), Path(pathData.Path)
{
  /* Do nothing */
}

double
ActionBaseClass::Action(int slice1, int slice2,
			const Array<int,1> &activeParticles, int level)
{
//   if (PathData.UseCorrelatedSampling()) {
//     if (level == 0) {
//       cerr << "Action should not be called at level 0 with "
// 	   << "correlated sampling on.\n";  
//       abort();
//     }
//     double actionA, actionB;
//     PathData.SetIonConfig(0);
//     actionA = SingleAction(slice1, slice2, activeParticles, level);
//     PathData.SetIonConfig(1);
//     actionB = SingleAction(slice1, slice2, activeParticles, level);
//     return 0.5*(actionA + actionB);
//   }
//   else 
    return SingleAction(slice1, slice2, activeParticles, level);
}

void
ActionBaseClass::GradAction(int slice1, int slice2,
			    const Array<int,1> &ptcls, int level,
			    Array<dVec,1> &gradVec)
{
  /// The default implementation does nothing for those actions which
  /// we do not wish to contribute to the force calculation.
}
