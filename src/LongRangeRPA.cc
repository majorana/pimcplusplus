#include "Common/Integration/RungeKutta.h"
#include "ActionClass.h"

inline void order (int species1, int species2)
{
  if (species1 > species2) {
    int tmp = species1;
    species1 = species2;
    species2 = tmp;
  }
}


inline int uindex (int species1, int species2, int numSpecies)
{
  order(species1, species2);
  return (species2*(species2+1)/2 + species1);
}

inline int windex (int species1, int species2, int numSpecies)
{
  int offset = numSpecies*(numSpecies+1)/2;
  return (offset+uindex(species1, species2, numSpecies));
}


/// Note:  we must set Level and ki 
/// uwvec is indexed in the following way, 
/// given s1= species1 and s2 = species2, with s1<=s2
/// u_(s1,s2) = uwvec(
Array<double,1> ActionClass::RPAIntegrand(double t, 
					  const Array<double,1> &uwvec)
{
  double levelTau = tau;
  for (int i=0; i<Level; i++)
    levelTau *= 2.0;

    Array<double,1> duwvec(uwvec.size());
  int numSpecies = Path.NumSpecies();
  dVec &k = Path.kVecs(ki);
  double k2 = dot (k, k);
  // First, calculate the \f$\dot{u}\f$'s.
  for (int species1=0; species1<numSpecies; species1++)
    for (int species2=species1; species2<numSpecies; species2++) {
      PairActionFitClass &pa12 = 
	*PairActionVector(PairMatrix(species1,species2));

      double vlongk = pa12.Ulong_k(Level, ki)/levelTau;
      //cerr << "vlong_k = " << vlongk << endl;
      int i12 = uindex(species1,species2,numSpecies);
      duwvec(i12) = 
	-pa12.lambda*k2*uwvec(i12) + vlongk;
      //cerr << "duwvec = " << duwvec(i12) << endl;
      for (int species3=0; species3<numSpecies; species3++) {
	int N3 = Path.Species(species3).NumParticles;
	double lambda3 = Path.Species(species3).lambda;
	double u23 = uwvec(uindex(species2, species3, numSpecies));
	double u13 = uwvec(uindex(species1, species3, numSpecies));
	double w23 = uwvec(windex(species2, species3, numSpecies));
	double w13 = uwvec(windex(species1, species3, numSpecies));
	duwvec(i12) -= 0.5*k2*lambda3*(u23*u13 /*+ w23*u13 + u23*w13*/);
      }
    }
  // Next, calculate the \f$\dot{w}\f$'s.
  for (int species1=0; species1<numSpecies; species1++)
    for (int species2=species1; species2<numSpecies; species2++) {
      PairActionFitClass &pa12 = 
	*PairActionVector(PairMatrix(species1,species2));
      double vshortk = pa12.Vk(sqrt(k2)) - pa12.Ulong_k(Level, ki)/levelTau;
      int i12 = windex(species1,species2,numSpecies);
      duwvec(i12) = 
	-pa12.lambda*k2*uwvec(i12) + vshortk;
      for (int species3=0; species3<numSpecies; species3++) {
	int N3 = Path.Species(species3).NumParticles;
	double lambda3 = Path.Species(species3).lambda;
	double u23 = uwvec(uindex(species2, species3, numSpecies));
	double u13 = uwvec(uindex(species1, species3, numSpecies));
	duwvec(i12) -= 0.5*k2*lambda3*(u23*u13);
      }
    }
  return duwvec;    
}

// All of the OptimizedBreakups must be computed before this is called.
void ActionClass::SetupRPA()
{
  for (int paIndex=0; paIndex<PairActionVector.size(); paIndex++) {
    PairActionFitClass &pa = *PairActionVector(paIndex);
    pa.U_RPA_long_k.resize(MaxLevels,Path.kVecs.size());
    pa.dU_RPA_long_k.resize(MaxLevels,Path.kVecs.size());
  }


  const int numPoints = 2000;
  double levelTau = tau;
  int m = Path.NumSpecies()*(Path.NumSpecies()+1);
  Array<double,2> uwvec(numPoints, m);

  // Setup integrator
  RungeKutta2<ActionClass> integrator(*this);

  // Calculated RPA for U
  for (Level=0; Level<MaxLevels; Level++) {
    LinearGrid tGrid(0.0, levelTau, numPoints);
    for (ki=0; ki<Path.kVecs.size(); ki++) {
      cerr << "ki = " << ki << endl;
      /// Set initial conditions
      for (int i=0; i<m; i++)
	uwvec(0, i) = 0.0;
      integrator.Integrate(tGrid, 0, numPoints-1, uwvec);
      for (int species1=0; species1<Path.NumSpecies(); species1++)
	for (int species2=species1; species2<Path.NumSpecies(); species2++) {
	  PairActionFitClass &pa = 
	    *PairActionVector(PairMatrix(species1, species2));
	  pa.U_RPA_long_k(Level, ki) = 
	    uwvec(numPoints-1,uindex(species1, species2, Path.NumSpecies()));
	  cerr << "Numerical value = " << pa.U_RPA_long_k(Level,ki) << endl;
	  double k2 = dot(Path.kVecs(ki), Path.kVecs(ki));
	  double lambda = pa.lambda;
	  double N = Path.Species(species1).NumParticles;
	  double vklong = pa.Ulong_k(Level, ki)/levelTau;
	  cerr << "uklong =          " << pa.Ulong_k(Level, ki) << endl;
	  double Uanalytic = (-1.0+sqrt(1.0+4.0*N*vklong/(lambda*k2)))/(2.0*N);
	  cerr << "One component ground state analytic = " << Uanalytic 
	       << endl;
	}
    
    levelTau *= 2.0;
    }
  }

//   // Calculated RPA for dU
//   for (Level=0; Level<MaxLevels; Level++) {
//     LinearGrid tGrid(0.0, levelTau, numPoints);
//     for (ki=0; ki<Path.kVecs.size(); ki++) {
//       /// Set initial conditions
//       for (int i=0; i<m; i++)
// 	uwvec(0, i) = 0.0;
//       integrator.Integrate(tGrid, 0, numPoints-1, uwvec);
//       for (int species1=0; species1<Path.NumSpecies(); species1++)
// 	for (int species2=species1; species2<Path.NumSpecies(); species2++) {
// 	  PairActionFitClass &pa = 
// 	    *PairActionVector(PairMatrix(species1, species2));
// 	  pa.U_RPA_long_k(Level, ki) = 
// 	    uwvec(numPoints-1,uindex(species1, species2, Path.NumSpecies()));
// 	}
    
//     levelTau *= 2.0;
//     }
//   }
}

