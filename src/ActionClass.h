#ifndef ACTION_CLASS
#define ACTION_CLASS

#include "PathClass.h"
#include "Common/PairAction/PAFit.h"
#include "Common/Ewald/Ewald.h"

/// This is the class that controls all of the actions and is in
/// charge of calculating them. When this is initialized a pointer needs
/// to be sent that has the memoizedData and SpeciesClass 
class PathDataClass;

class ActionClass
{
private:
  PathDataClass &PathData;
  PathClass &Path;
  double CalcLRAction(int slice, int level);
public:
  /// This holds all of the Pair Action Classes
  Array<PairActionFitClass*,1> PairActionVector;
  Array<EwaldClass*,1> EwaldVector;
  /// Holds indices to which PairActionClass in the PairAcctionVector
  /// you use for a given pair of particles indexed by
  /// (species1,species2) 
  Array<int,2> PairMatrix;
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
  void PrintDensityMatrix();
  ActionClass(PathDataClass  &pathdata);
};



#endif
