#include "LongRangePotClass.h"
#include "../PathDataClass.h"

LongRangePotClass::LongRangePotClass 
(PathDataClass &pathData, Array<PairActionFitClass*,2> &pairMatrix) :
  PotentialBaseClass (pathData), PairMatrix(pairMatrix)
{
  // Do nothing 
}

inline double mag2(complex<double> z)
{
  return z.real()*z.real()+z.imag()*z.imag();
}

double LongRangePotClass::V(int slice)
{
  double homo = 0.0;
  double hetero = 0.0;
  double background = 0.0;
  double k0Terms = 0.0;
  if (PathData.Actions.HaveLongRange()) {
    PathClass &Path = PathData.Path;
    // First, do the homologous (same species) terms
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      PairActionFitClass &pa = *PairMatrix(species,species);
      if (pa.IsLongRange()) {
	for (int ki=0; ki<Path.kVecs.size(); ki++) {
	  double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	  homo += 0.5 * 2.0* rhok2 * pa.Vlong_k(ki);
	}
      }
      int N = Path.Species(species).NumParticles;
      // We can't forget the Madelung term.
      homo -= 0.5 * N * pa.Vlong_r0;
      // Or the neutralizing background term
      background -= 0.5*N*N*pa.Vshort_k0;
      k0Terms += 0.5*N*pa.Vlong_k0;
    }
    
    // Now do the heterologous terms
    for (int species1=0; species1<Path.NumSpecies(); species1++)
      for (int species2=species1+1; species2<Path.NumSpecies(); species2++) {
	PairActionFitClass &pa= *PairMatrix(species1, species2);
	if (pa.IsLongRange()) {
	  for (int ki=0; ki<Path.kVecs.size(); ki++) {
	    double rhorho = 
	      Path.Rho_k(slice, species1, ki).real() *
	      Path.Rho_k(slice, species2, ki).real() + 
	      Path.Rho_k(slice, species1, ki).imag() *
	      Path.Rho_k(slice, species2, ki).imag();
	    hetero += 2.0 * rhorho * pa.Vlong_k(ki);
	  }
	  int N1 = Path.Species(species1).NumParticles;
	  int N2 = Path.Species(species2).NumParticles;
	  background -= N1*N2*pa.Vshort_k0;
	  k0Terms += N1*N2*pa.Vlong_k0;
	}
      }
  }
  return (homo+hetero+k0Terms/*+background*/);
}
