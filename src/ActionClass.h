#ifndef ACTION_CLASS
#define ACTION_CLASS

#include "Common/Splines/CubicSpline.h"


#include "SpeciesClass.h"
#include "PathClass.h"
//#include "DistanceTablePBCClass.h"
//#include "DistanceTableFreeClass.h"
#include "DistanceTableClass.h"
#include "Common/PairAction/PAFit.h"


/// This is the class that controls all of the actions and is in
/// charge of calculating them. When this is initialized a pointer needs
/// to be sent that has the memoizedData and SpeciesClass 
class ActionClass
{
private:
public:
  DistanceTableClass *DistanceTable;
  /// This holds all of the Pair Action Classes
  Array<PairActionFitClass*,1> PairActionVector;
  /// Holds indices to which PairActionClass in the PairAcctionVector
  /// you use for a given pair of particles indexed by
  /// (species1,species2) 
  Array<int,2> PairMatrix;
  PathClass &Path;
  /// Temperature
  double tau;
  /// The maximum number of levels we can have in a bisection move.
  int MaxLevels;
  void Read(IOSectionClass &IOSection);
  /// Calculates the total action.
  double calcTotalAction(int startSlice, int endSlice, 
			 Array<int,1> changedParticles,
			 int level);

  /// Function to calculate the total action.
  void calcTotalAction();
  ActionClass(PathClass  &p_path) : Path(p_path) 
  {
  }
};



#endif
