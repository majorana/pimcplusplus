#ifndef ACTION_CLASS
#define ACTION_CLASS

#include "PathClass.h"
#include "Common/PairAction/PAFit.h"
#include "Common/Ewald/Ewald.h"

class PathDataClass;
typedef enum {DO_U, DO_DU, DO_V} TaskType; 

/// This is the class that controls all of the actions and is in
/// charge of calculating them. When this is initialized a pointer needs
/// to be sent that has the memoizedData and SpeciesClass 
class ActionClass
{
private:
  PathDataClass &PathData;
  PathClass &Path;
  LinearGrid LongGrid;
  /// This calculates the quantity 
  /// \f$ X_k \equiv -\frac{4\pi}{\Omega k} \int_{r_c}^\infty dr \, r \sin(kr) V(r).\f$
  double CalcXk (int paIndex, int level, double k, double rc, TaskType type);
  // This must be called after all of the OptimizedBreakup_x's
  Array<double,1> RPAIntegrand(double t, const Array<double,1> &Uvec);
  void SetupRPA();
  void TestRPA();
  bool RPATaskIsU;
  int Level, ki;
  bool UseRPA;
public:
  inline Array<double,1> operator()(double t, Array<double,1> uwvec)
  { return RPAIntegrand(t, uwvec); }
  double LongRange_U(int slice, int level);
  double LongRange_U_RPA(int slice, int level);
  double LongRange_dU(int slice, int level);
  double LongRange_dU_RPA(int slice, int level);
  double LongRange_V(int slice);
  /// This holds all of the Pair Action Classes
  Array<PairActionFitClass*,1> PairActionVector;

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
  double UAction (int startSlice, int endSlice, 
		  const Array<int,1> &changedParticles, int level);
  double KAction (int startSlice, int endSlice, 
		  const Array<int,1> &changedParticles, int level);
  double TotalAction(int startSlice, int endSlice, 
		     const Array<int,1> &changedParticles, int level);

  void Energy (int slice1, int level,
	       double &spring, double &dU);
  double PotentialEnergy (int slice);

  /// This computes the optimized breakups for the pair actions stored
  /// in PairActionVector.  The parameters are the number of knots in
  /// the "spline" representation of the long-range action and the
  /// k-space cutoff.  
  /// Only \f$\mathbf{k}\f$ with \f$|\mathbf{k}| < k_c$\f will be
  /// included in the simulation sum.
  void OptimizedBreakup_U(int numKnots);
  void OptimizedBreakup_dU(int numKnots);
  void OptimizedBreakup_V(int numKnots);
  ActionClass(PathDataClass  &pathdata);
};



#endif
