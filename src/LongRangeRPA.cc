#include "Common/Integration/RungeKutta.h"
#include "ActionClass.h"

inline void order (int &species1, int &species2)
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

  double boxVol = Path.GetVol();
  Array<double,1> duwvec(uwvec.size());
  int numSpecies = Path.NumSpecies();
  dVec &k = Path.kVecs(ki);
  double k2 = dot (k, k);
  // First, calculate the \f$\dot{u}\f$'s.
  for (int species1=0; species1<numSpecies; species1++)
    for (int species2=species1; species2<numSpecies; species2++) {
      PairActionFitClass &pa12 = 
	*PairActionVector(PairMatrix(species1,species2));

      double vlongk;
      if (RPATaskIsU)
	vlongk = pa12.Ulong_k(Level, ki)/levelTau;
      else
	vlongk = pa12.dUlong_k(Level,ki);
      int i12 = uindex(species1,species2,numSpecies);
      double lambda12 = 
	0.5* (Path.Species(species1).lambda + Path.Species(species2).lambda);
      duwvec(i12) = -lambda12*k2*uwvec(i12) + vlongk;
      for (int species3=0; species3<numSpecies; species3++) {
	int N3 = Path.Species(species3).NumParticles;
	double lambda3 = Path.Species(species3).lambda;
	double u23 = uwvec(uindex(species2, species3, numSpecies));
	double u13 = uwvec(uindex(species1, species3, numSpecies));
	double w23 = uwvec(windex(species2, species3, numSpecies));
	double w13 = uwvec(windex(species1, species3, numSpecies));
	duwvec(i12) -= 0.5*k2*N3*lambda3*(u23*u13 + w23*u13 + u23*w13);
      }
    }
  // Next, calculate the \f$\dot{w}\f$'s.
  for (int species1=0; species1<numSpecies; species1++)
    for (int species2=species1; species2<numSpecies; species2++) {
      PairActionFitClass &pa12 = 
	*PairActionVector(PairMatrix(species1,species2));
      double k = sqrt(k2);
      double vlongk;
      if (RPATaskIsU)
	vlongk = pa12.Ulong_k(Level, ki)/levelTau;
      else
	vlongk = pa12.dUlong_k(Level,ki);
      double vshortk = pa12.Vk(k)/boxVol - vlongk;
      int i12 = windex(species1,species2,numSpecies);
      double lambda12 = 
	0.5* (Path.Species(species1).lambda + Path.Species(species2).lambda);
      duwvec(i12) = -lambda12*k2*uwvec(i12) + vshortk;
      for (int species3=0; species3<numSpecies; species3++) {
	int N3 = Path.Species(species3).NumParticles;
	double lambda3 = Path.Species(species3).lambda;
	double w23 = uwvec(windex(species2, species3, numSpecies));
	double w13 = uwvec(windex(species1, species3, numSpecies));
	duwvec(i12) -= 0.5*k2*N3*lambda3*(w23*w13);
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

  double boxVol = Path.GetVol();
  // Calculated RPA for U
  RPATaskIsU = true;
  for (Level=0; Level<MaxLevels; Level++) {
    LinearGrid tGrid(0.0, levelTau, numPoints);
    for (ki=0; ki<Path.kVecs.size(); ki++) {
      //cerr << "ki = " << ki << endl;
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
	  double k2 = dot(Path.kVecs(ki), Path.kVecs(ki));
	  double k = sqrt(k2);
	  //cerr << "k = " << k << endl;
	  //cerr << "Numerical value = " << pa.U_RPA_long_k(Level,ki) << endl;
	  double lambda = pa.lambda;
	  double N = Path.Species(species1).NumParticles;
	  double vklong = pa.Ulong_k(Level, ki)/levelTau;
	  //cerr << "uklong  =         " << pa.Ulong_k(Level, ki) << endl;
	  //cerr << "vshort  =         " << 
	  //  pa.Vk(k)/boxVol - pa.Ulong_k(Level, ki)/levelTau << endl; 
	  //cerr << "ukshort =         " << 
	  //  uwvec(numPoints-1,windex(species1,species2,Path.NumSpecies()))
	  //     << endl;
	  double Uanalytic = (-1.0+sqrt(1.0+4.0*N*vklong/(lambda*k2)))/(2.0*N);
	  //cerr << "One component ground state analytic = " << Uanalytic 
	  //     << endl;
	}
    
    levelTau *= 2.0;
    }
  }

  // Calculated RPA for dU
  RPATaskIsU = false;
  for (Level=0; Level<MaxLevels; Level++) {
    LinearGrid tGrid(0.0, levelTau, numPoints);
    for (ki=0; ki<Path.kVecs.size(); ki++) {
      /// Set initial conditions
      for (int i=0; i<m; i++)
	uwvec(0, i) = 0.0;
      integrator.Integrate(tGrid, 0, numPoints-1, uwvec);
      Array<double,1> duvec(m);
      duvec = RPAIntegrand (levelTau, uwvec(numPoints-1,Range::all()));
      for (int species1=0; species1<Path.NumSpecies(); species1++)
	for (int species2=species1; species2<Path.NumSpecies(); species2++) {
	  PairActionFitClass &pa = 
	    *PairActionVector(PairMatrix(species1, species2));
	  pa.dU_RPA_long_k(Level, ki) = 
	    duvec (uindex(species1, species2, Path.NumSpecies()));
	}
    } 
  }
  TestRPA();
}


void ActionClass::TestRPA()
{
  FILE *fout = fopen ("RPAtest.dat", "w");
  double minL = Path.GetBox()[0];
  for (int i=1; i<NDIM; i++)
    minL = min(minL, Path.GetBox()[i]);
  int numPoints = 1000;
  LinearGrid rgrid(0.0, 0.5*minL, numPoints);
  int m = (Path.NumSpecies()*(Path.NumSpecies()+1))/2;
  Array<double,2> u(m,numPoints), du(m,numPoints), 
    uRPA(m,numPoints), duRPA(m, numPoints);
  u = 0.0; du = 0.0; uRPA = 0.0; duRPA = 0.0;
  for (int ri=0; ri<rgrid.NumPoints; ri++) {
    double r = rgrid(ri);
    dVec vr;
    vr[0] = r; vr[1] = 0.0; vr[2] = 0.0;
    for (int ki=0; ki<Path.kVecs.size(); ki++) {
      dVec &vk = Path.kVecs(ki);
      double coskr = cos(dot(vk, vr));
      
      for (int species1=0; species1<Path.NumSpecies(); species1++) 
	for (int species2 = species1; species2<Path.NumSpecies(); species2++) {
	  PairActionFitClass &pa = 
	    *PairActionVector(PairMatrix(species1, species2));
	  int si = uindex(species1, species2, Path.NumSpecies());
	  u(si, ri)     += pa.Ulong_k(0, ki) * 2.0*(1.0+coskr);
	  du(si, ri)    += pa.dUlong_k(0, ki) * 2.0*(1.0+coskr);
	  uRPA(si, ri)  += pa.U_RPA_long_k(0, ki) * 2.0*(1.0+coskr);
	  duRPA(si, ri) += pa.dU_RPA_long_k(0, ki) * 2.0*(1.0+coskr);
	}
    }
    fprintf (fout, "%1.12e", r);
    for (int si=0; si<m; si++) 
      fprintf (fout, " %1.12e %1.12e %1.12e %1.12e", 
	       u(si, ri), uRPA(si, ri), du(si, ri), duRPA(si,ri));
    fprintf (fout, "\n");
    
  }
  fclose (fout);
}
