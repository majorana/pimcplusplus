#include "LongRangeClass.h"
#include "../PathDataClass.h"
#include "../Common/Ewald/OptimizedBreakup.h"


#include "../Common/Ewald/OptimizedBreakup.h"
#include "../Common/Integration/GKIntegration.h"

class CoulombXkIntegrand
{
private:
  PairActionFitClass &PA;
  int Level;
  double k;
  double beta;
  JobType Task;

  inline double Uintegrand(double r)
  {
    double U = PA.Udiag(r, Level);
    U -= beta * PA.Z1Z2/r;
    return r * sin(k*r)*U;
  }
  inline double dUintegrand(double r)
  {
    double dU = PA.dUdiag(r, Level);
    dU -= PA.Z1Z2/r;
    return r * sin(k*r)*dU;
  }
  inline double Vintegrand(double r)
  {
    double V = PA.V(r);
    V -= PA.Z1Z2/r;
    return r * sin(k*r)*V;
  }
  
public:

  inline double operator()(double r)
  {
    if (Task == JOB_U)
      return Uintegrand(r);
    else if (Task == JOB_DU)
      return dUintegrand(r);
    else
      return Vintegrand(r);
  }

  CoulombXkIntegrand (PairActionFitClass &pa, int level, double k_,
		      JobType task) :
    PA(pa), Level(level), k(k_), Task(task)
  { 
    beta = pa.SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
  }
};

class XkIntegrand
{
private:
  PairActionFitClass &PA;
  int Level;
  double k;
  JobType Task;
  inline double Uintegrand(double r)
  {
    double U = PA.Udiag(r, Level);
    return r * sin(k*r)*U;
  }
  inline double dUintegrand(double r)
  {
    double dU = PA.dUdiag(r, Level);
    return r * sin(k*r)*dU;
  }
  inline double Vintegrand(double r)
  {
    double V = PA.V(r);
    return r * sin(k*r)*V;
  } 

public:
  inline double operator()(double r) 
  {
    if (Task == JOB_U)
      return Uintegrand(r);
    else if (Task == JOB_DU)
      return dUintegrand(r);
    else
      return Vintegrand(r);
 
  }
  XkIntegrand (PairActionFitClass &pa, int level, double k_,
	       JobType task) :
    PA(pa), Level(level), k(k_), Task(task)
  { /* do nothing else*/  }
};



/// This calculates the quantity 
/// \f$ X_k \equiv -\frac{4 \pi}{\Omega k} \int_{r_c}^\infty dr \, r \sin(kr) V(r).\f$
double LongRangeClass::CalcXk (int paIndex, int level, double k, double rc,
			    JobType task)
{
  PathClass &Path=PathData.Path;
  double absTol = 1.0e-7;
  double relTol = 1.0e-5;

  PairActionFitClass &pa = *PairArray(paIndex);
  if (task == JOB_U)
    absTol *= pa.SmallestBeta*pow(2.0, (double)level);

  double Xk;
  if (pa.Z1Z2 == 0.0) {
    XkIntegrand integrand(pa, level, k, task);
    GKIntegration<XkIntegrand, GK31> integrator(integrand);
    Xk = 4.0*M_PI/(Path.GetVol()*k) * 
      integrator.Integrate(rc, 20.0*rc, absTol, relTol, false);
    return Xk;
  }
  else {
    CoulombXkIntegrand integrand(pa, level, k, task);
    GKIntegration<CoulombXkIntegrand, GK31> integrator(integrand);
    // integrator.SetRelativeErrorMode();
    
    
    if (false/*task != JOB_V*/) {
      integrator.SetRelativeErrorMode();
      Xk = -4.0*M_PI/(Path.GetVol()*k) * 
	integrator.Integrate (rc, 20.0*rc, relTol);    
    }
    else
      Xk = -4.0*M_PI/(Path.GetVol()*k) * 
	integrator.Integrate(rc, 20.0*rc, absTol, relTol, false);
    
    /// Add in the analytic part that I ignored
    /// Multiply analytic term by tau only for U -- do not multiply
    /// for dU or V.
    
    double coef;
    if (task == JOB_U) {
      coef = pa.SmallestBeta;
      for (int i=0; i<level; i++)
	coef *= 2.0;
    }
    else
      coef = 1.0;

    Xk -= coef*4.0*M_PI*pa.Z1Z2/(Path.GetVol()*k*k)*cos(k*rc);
    return (Xk);
  }
}


class UshortIntegrand
{
private:
  PairActionFitClass &PA;
  int Level;
  JobType Task;
  inline double Uintegrand(double r)
  {
    double Ushort = PA.Udiag(r, Level) - PA.Ulong(Level)(r);
    return r*r*Ushort;
  }
  inline double dUintegrand(double r)
  {
    double dUshort = PA.dUdiag(r, Level) - PA.dUlong(Level)(r);
    return r*r*dUshort;
  }
  inline double Vintegrand(double r)
  {
    double Vshort = PA.V(r) - PA.Vlong(r);
    return r*r*Vshort;
  } 

public:
  inline double operator()(double r) 
  {
    if (Task == JOB_U)
      return Uintegrand(r);
    else if (Task == JOB_DU)
      return dUintegrand(r);
    else
      return Vintegrand(r);
 
  }
  UshortIntegrand (PairActionFitClass &pa, int level,
	       JobType task) :
    PA(pa), Level(level), Task(task)
  { /* do nothing else*/  }
};



void LongRangeClass::Read(IOSectionClass& in)
{
  //do nothing for now
}


LongRangeClass::LongRangeClass(PathDataClass &pathData,
			       Array<PairActionFitClass* ,2> &pairMatrix,
			       Array<PairActionFitClass*, 1> &pairArray) : 
  ActionBaseClass (pathData),
  PairMatrix(pairMatrix),
  PairArray(pairArray)
{
}

inline double mag2 (const complex<double> &z)
{
  return (z.real()*z.real() + z.imag()*z.imag());
}

double LongRangeClass::Action (int slice1, int slice2,
			       const Array<int,1> &changedParticles,
			       int level)
{
  double homo = 0.0;
  double hetero = 0.0;
  int skip = (1<<level);
  double levelTau = Path.tau * (double)skip;
  for (int slice=slice1; slice<=slice2; slice+=skip) {
    double factor;
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    // First, do the homologous (same species) terms
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      PairActionFitClass &pa = *PairMatrix(species,species);
      if (pa.IsLongRange()) {
	for (int ki=0; ki<Path.kVecs.size(); ki++) {
	  double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	  homo += factor * 0.5 * 2.0* rhok2 * pa.Ulong_k(level,ki);
	}
      }
      int N = Path.Species(species).NumParticles;
      // We can't forget the Madelung term.
      homo -= factor * 0.5 * N * pa.Ulong_r0(level);
      // Or the neutralizing background term
      homo -= factor * 0.5*N*N*pa.Ushort_k0(level);
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
	    hetero += factor * 2.0 * rhorho * pa.Ulong_k(level,ki);
	  }
	  int N1 = Path.Species(species1).NumParticles;
	  int N2 = Path.Species(species2).NumParticles;
	  hetero -= factor * N1*N2*pa.Ushort_k0(level);
	}
      }
  }
  return (homo+hetero);
}


double LongRangeClass::d_dBeta (int slice1, int slice2,  int level)
{
  double homo = 0.0;
  double hetero = 0.0;
  int skip = (1<<level);
  double levelTau = Path.tau * (double)skip;
  for (int slice=slice1; slice<=slice2; slice+=skip) {
    double factor;
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    // First, do the homologous (same species) terms
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      PairActionFitClass &PA = *PairMatrix(species,species);
      if (PA.IsLongRange()) {
	for (int ki=0; ki<Path.kVecs.size(); ki++) {
	  double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	  homo += factor* 0.5 * 2.0* rhok2 * PA.dUlong_k(level,ki);
	}
      }
      int N = Path.Species(species).NumParticles;
      // We can't forget the Madelung term.
      homo -= 0.5 * N * PA.dUlong_r0(level);
      // Or the neutralizing background term
      homo -= 0.5*N*N*PA.dUshort_k0(level);
    }
    
    // Now do the heterologous terms
    for (int species1=0; species1<Path.NumSpecies(); species1++)
      for (int species2=species1+1; species2<Path.NumSpecies(); species2++) {
	PairActionFitClass &PA = *PairMatrix(species1,species2);
	if (PA.IsLongRange()) {
	  for (int ki=0; ki<Path.kVecs.size(); ki++) {
	    double rhorho = 
	      Path.Rho_k(slice, species1, ki).real() *
	      Path.Rho_k(slice, species2, ki).real() + 
	      Path.Rho_k(slice, species1, ki).imag() *
	      Path.Rho_k(slice, species2, ki).imag();
	    hetero += factor* 2.0 * rhorho * PA.dUlong_k(level,ki);
	  }
	  int N1 = Path.Species(species1).NumParticles;
	  int N2 = Path.Species(species2).NumParticles;
	  hetero -= N1*N2*PA.dUshort_k0(level);
	}
      }
  }
  return (homo+hetero);
    
}

void LongRangeClass::Init(IOSectionClass &in)
{
  if (PathData.Path.Getkc() == 0.0) {
    cerr << "Missing kCutoff in System section.  Aborting.\n";
    abort();
  }
  int numKnots;
  assert(in.ReadVar ("NumBreakupKnots", numKnots));
  cerr << "Doing optimized long range breakups...\n";
  OptimizedBreakup_U(numKnots);
  OptimizedBreakup_dU(numKnots);
  OptimizedBreakup_V(numKnots);
}




/// This computes the optimized breakups for the pair actions stored
/// in PairArray.  The parameters are the number of knots in
/// the "spline" representation of the long-range action and the
/// k-space cutoff.  
/// Only \f$\mathbf{k}\f$ with \f$|\mathbf{k}| < k_c$\f will be
/// included in the simulation sum.
void LongRangeClass::OptimizedBreakup_U(int numKnots)
{
  /// BUG: Long Range Optimized Breakups only work for NDIM=3
#if NDIM==3
  PathClass &Path=PathData.Path;
  const double tolerance = 1.0e-7;
  double kCut = Path.Getkc();
  dVec box = Path.GetBox();
  double boxVol = box[0]*box[1]*box[2];
  double rc = 0.5*box[0];
  for (int i=1; i<NDIM; i++)
    rc = min (rc, 0.5*box[i]);
  double kvol = Path.GetkBox()[0];
  for (int i=1; i<NDIM; i++)
    kvol *= Path.GetkBox()[i];
  double kavg = pow(kvol,1.0/3.0);

  

  LPQHI_BasisClass basis;
  basis.Set_rc(rc);
  basis.SetBox(box);
  basis.SetNumKnots (numKnots);

  // We try to pick kcont to keep reasonable number of k-vectors
  double kCont = 50.0 * kavg;
  double delta = basis.GetDelta();
  double kMax = 20.0*M_PI/delta;
  cerr << "kCont = " << kCont 
       << " kMax = " << kMax << endl;

  OptimizedBreakupClass breakup(basis);
  breakup.SetkVecs (kCut, kCont, kMax);
  int numk = breakup.kpoints.size();
  int N = basis.NumElements();
  Array<double,1> t(N);
  Array<bool,1>   adjust (N);
  Array<double,1> Xk(numk);

  // Would be 0.5, but with two timeslice distdisp, it could be a
  // little longer
  double rmax = 0.75 * sqrt (dot(box,box));
  const int numPoints = 1000;
  LongGrid.Init (0.0, rmax, numPoints);
  Array<double,1> Ulong_r(numPoints);

  for (int paIndex=0; paIndex<PairArray.size(); paIndex++) {
    PairActionFitClass &pa = *PairArray(paIndex);
    pa.Setrc (rc);
    pa.Ulong.resize(pa.NumBetas);
    pa.Ulong_k.resize(pa.NumBetas,Path.kVecs.size());
    pa.Ulong_k = 0.0;
    pa.Ulong_r0.resize(pa.NumBetas);
    pa.Ushort_k0.resize(pa.NumBetas);
    for (int level=0; level<pa.NumBetas; level++) {
      Ulong_r = 0.0;
      
      // Calculate Xk's
      cerr << "Calculating Xk's for U...\n";
      for (int ki=0; ki<numk; ki++) {
	//Xk(ki) = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc, JOB_U);
	//double oldXk = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc, JOB_U);
	double k = breakup.kpoints(ki)[0];
	Xk(ki) = pa.Xk_U (k, level) / boxVol;
      }
      cerr << "Done.\n";
      

      // Set boundary conditions at rc:  Force value and first and
      // second derivatives of long-range potential to match the full
      // potential at rc.
      adjust = true;
      /// Warning:  the following constraints may cause instabilities!!!!
      //  t(N-3) = pa.Udiag(rc, level);                 adjust(N-3) = false;
      //  t(N-2) = pa.Udiag_p(rc, level)*delta;         adjust(N-2) = false;
      //  t(N-1) = pa.Udiag_pp(rc, level)*delta*delta;  adjust(N-1) = false;
      //  t(1) = 0.0;                                   adjust(1)   = false;
      
      // Now, do the optimal breakup:  this gives me the coefficents
      // of the basis functions, h_n in the array t.
      cerr << "Doing U breakup...\n";
      breakup.DoBreakup (Xk, t, adjust);
      cerr << "Done.\n";
      
      // Now, we must put this information into the pair action
      // object.  First do real space part
      pa.Ulong_r0(level)=0.0;
      for (int n=0; n<N; n++)
	pa.Ulong_r0(level) += t(n)*basis.h(n,0.0);
      for (int i=0; i<LongGrid.NumPoints; i++) {
	double r = LongGrid(i);
	if (r <= rc) {
	  // Sum over basis functions
	  for (int n=0; n<N; n++) 
	    Ulong_r(i) += t(n) * basis.h(n, r);
	}
	else
	  Ulong_r(i) = pa.Udiag (r, level);
      }
      pa.Ulong(level).Init(&LongGrid, Ulong_r);

      // Calculate FT of Ushort at k=0
      UshortIntegrand integrand(pa, level, JOB_U);
      GKIntegration<UshortIntegrand, GK31> integrator(integrand);
      integrator.SetRelativeErrorMode();
      pa.Ushort_k0(level) = 4.0*M_PI/boxVol * 
	integrator.Integrate(0.0, rc, tolerance);
      cerr << "Ushort_k0(" << level << ") = " << pa.Ushort_k0(level) << endl;

      // Now do k-space part
      for (int ki=0; ki < Path.kVecs.size(); ki++) {
	const dVec &kv = Path.kVecs(ki);
	double k = sqrt (dot(kv,kv));
	// Sum over basis functions
	for (int n=0; n<N; n++)
	  pa.Ulong_k(level,ki) += t(n) * basis.c(n,k);
	// Now add on part from rc to infinity
	// pa.Ulong_k(level,ki) -= CalcXk(paIndex, level, k, rc, JOB_U);
	pa.Ulong_k(level,ki) -= pa.Xk_U(k, level)/boxVol;
      }
//       // HACK HACK HACK HACK
//       FILE *fout = fopen ("Vlongk.dat", "w");
//       for (double k=0; k<50.0; k+=0.01) {
// 	double U = 0.0;
// 	// Sum over basis functions
// 	for (int n=0; n<N; n++)
// 	  U += t(n) * basis.c(n,k);
// 	fprintf (fout, "%1.16e %1.16e ", k, U);
// 	// Now add on part from rc to infinity
// 	U -= CalcXk(paIndex, level, k, rc);
// 	fprintf (fout, "%1.16e \n", U);
//       }
//       fclose (fout);
    }
  }
#endif
}

void LongRangeClass::OptimizedBreakup_dU(int numKnots)
{
  ///BUG: Optimized Breakup only works when NDIM==3
#if NDIM==3
  PathClass &Path=PathData.Path;
  const double tolerance = 1.0e-7;
  double kCut = Path.Getkc();
  dVec box = Path.GetBox();
  double boxVol = box[0]*box[1]*box[2];
  double rc = 0.5*box[0];
  for (int i=1; i<NDIM; i++)
    rc = min (rc, 0.5*box[i]);
  double kvol = Path.GetkBox()[0];
  for (int i=1; i<NDIM; i++)
    kvol *= Path.GetkBox()[i];
  double kavg = pow(kvol,1.0/3.0);

  LPQHI_BasisClass basis;
  basis.Set_rc(rc);
  basis.SetBox(box);
  basis.SetNumKnots (numKnots);

  // We try to pick kcont to keep reasonable number of k-vectors
  double kCont = 50.0 * kavg;
  double delta = basis.GetDelta();
  double kMax = 20.0*M_PI/delta;
  cerr << "kCont = " << kCont 
       << " kMax = " << kMax << endl;

  OptimizedBreakupClass breakup(basis);
  breakup.SetkVecs (kCut, kCont, kMax);
  int numk = breakup.kpoints.size();
  int N = basis.NumElements();
  Array<double,1> t(N);
  Array<bool,1>   adjust (N);
  Array<double,1> Xk(numk);

  // Would be 0.5, but with two timeslice distdisp, it could be a
  // little longer
  double rmax = 0.75 * sqrt (dot(box,box));
  const int numPoints = 1000;
  LongGrid.Init (0.0, rmax, numPoints);
  Array<double,1> dUlong_r(numPoints);

  for (int paIndex=0; paIndex<PairArray.size(); paIndex++) {
    PairActionFitClass &pa = *PairArray(paIndex);
    pa.Setrc (rc);
    pa.dUlong.resize(pa.NumBetas);
    pa.dUlong_k.resize(pa.NumBetas,Path.kVecs.size());
    pa.dUlong_k = 0.0;
    pa.dUlong_r0.resize(pa.NumBetas);
    pa.dUshort_k0.resize(pa.NumBetas);
    for (int level=0; level<pa.NumBetas; level++) {
      dUlong_r = 0.0;

      // Calculate Xk's
      cerr << "Calculating Xk's for dU...\n";
      for (int ki=0; ki<numk; ki++) {	
	// Xk(ki) = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc, JOB_DU);
	double k = breakup.kpoints(ki)[0];
	Xk(ki) = pa.Xk_dU (k, level) / boxVol;
      }
      cerr << "Done.\n";
      // Set boundary conditions at rc:  Force value and first and
      // second derivatives of long-range potential to match the full
      // potential at rc.
      adjust = true;
      double delta = basis.GetDelta();
      /// Warning:  the following constraints may cause instabilities!!!!
//       t(N-3) = pa.dUdiag(rc, level);                 adjust(N-3) = false;
//       t(N-2) = pa.dUdiag_p(rc, level)*delta;         adjust(N-2) = false;
//       t(N-1) = pa.dUdiag_pp(rc, level)*delta*delta;  adjust(N-1) = false;
//       t(1) = 0.0;                                    adjust(1)   = false;

      // Now, do the optimal breakup:  this gives me the coefficents
      // of the basis functions, h_n in the array t.
      cerr << "Doing dU breakup...\n";
      breakup.DoBreakup (Xk, t, adjust);
      cerr << "Done.\n";
      
      // Now, we must put this information into the pair action
      // object.  First do real space part
      pa.dUlong_r0(level)=0.0;
      for (int n=0; n<N; n++)
	pa.dUlong_r0(level) += t(n)*basis.h(n,0.0);
      for (int i=0; i<LongGrid.NumPoints; i++) {
	double r = LongGrid(i);
	if (r <= rc) {
	  // Sum over basis functions
	  for (int n=0; n<N; n++) 
	    dUlong_r(i) += t(n) * basis.h(n, r);
	}
	else
	  dUlong_r(i) = pa.dUdiag (r, level);
      }
      pa.dUlong(level).Init(&LongGrid, dUlong_r);

      // Calculate FT of Ushort at k=0
      UshortIntegrand integrand(pa, level, JOB_DU);
      GKIntegration<UshortIntegrand, GK31> integrator(integrand);
      integrator.SetRelativeErrorMode();
      pa.dUshort_k0(level) = 4.0*M_PI/boxVol * 
	integrator.Integrate(0.0, rc, tolerance);
      cerr << "dUshort_k0(" << level << ") = " << pa.dUshort_k0(level) << endl;

      // Now do k-space part
      for (int ki=0; ki < Path.kVecs.size(); ki++) {
	const dVec &kv = Path.kVecs(ki);
	double k = sqrt (dot(kv,kv));
	// Sum over basis functions
	for (int n=0; n<N; n++)
	  pa.dUlong_k(level,ki) += t(n) * basis.c(n,k);
	// Now add on part from rc to infinity
	//pa.dUlong_k(level,ki) -= CalcXk(paIndex, level, k, rc, JOB_DU);
	pa.dUlong_k(level,ki) -= pa.Xk_dU(k, level) / boxVol;
      }
    }
  }
#endif
}



// void LongRangeClass::OptimizedBreakup_dU(int numKnots)
// {
//   double kCut = Path.Getkc();
//   dVec box = Path.GetBox();
//   double rc = 0.5*box[0];
//   for (int i=1; i<NDIM; i++)
//     rc = min (rc, 0.5*box[i]);
//   double kvol = Path.GetkBox()[0];
//   for (int i=1; i<NDIM; i++)
//     kvol *= Path.GetkBox()[i];
//   double kavg = pow(kvol,1.0/3.0);
//   // We try to pick kcont to keep reasonable number of k-vectors
//   double kCont = 50.0 * kavg;
//   double kMax = 100 * kavg;
//   cerr << "kCont = " << kCont 
//        << " kMax = " << kMax << endl;

//   LPQHI_BasisClass basis;
//   basis.Set_rc(rc);
//   basis.SetBox(box);
//   basis.SetNumKnots (numKnots);

//   OptimizedBreakupClass breakup(basis);
//   breakup.SetkVecs (kCut, kCont, kMax);
//   int numk = breakup.kpoints.size();
//   int N = basis.NumElements();
//   Array<double,1> t(N);
//   Array<bool,1>   adjust (N);
//   Array<double,1> Xk(numk);

//   // Would be 0.5, but with two timeslice distdisp, it could be a
//   // little longer
//   double rmax = 0.75 * sqrt (dot(box,box));
//   const int numPoints = 1000;
//   LongGrid.Init (0.0, rmax, numPoints);
//   Array<double,1> dUlong_r(numPoints);

//   for (int paIndex=0; paIndex<PairArray.size(); paIndex++) {
//     PairActionFitClass &pa = *PairArray(paIndex);
//     pa.dUlong.resize(pa.NumBetas);
//     pa.dUlong_k.resize(pa.NumBetas,Path.kVecs.size());
//     pa.dUlong_k = 0.0;
//     pa.dUlong_r0.resize(pa.NumBetas);
//     for (int level=0; level<pa.NumBetas; level++) {
//       dUlong_r = 0.0;
//       // Calculate Xk's
//       for (int ki=0; ki<numk; ki++)
// 	Xk(ki) = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc,
// 			JOB_DU);

//       // Set boundary conditions at rc:  Force value and first and
//       // second derivatives of long-range potential to match the full
//       // potential at rc.
//       adjust = true;
//       double delta = basis.GetDelta();
//       cerr << "dUdiag(rc) = " << pa.dUdiag(rc, level) << endl;
//       cerr << "V(rc) = " << pa.V(rc) << endl;
//       cerr << "dUdiag_p(rc) = " << pa.dUdiag_p(rc, level) << endl;
//       cerr << "Vp(rc) = " << pa.Vp(rc) << endl;
//       cerr << "dUdiag_pp(rc) = " << pa.dUdiag_pp(rc, level) << endl;
//       cerr << "Vpp(rc) = " << pa.Vpp(rc) << endl;

      /// Warning:  the following constraints may cause instabilities!!!!
//       // t(N-3) = pa.dUdiag(rc, level);                 adjust(N-3) = false;
//       //t(N-2) = pa.dUdiag_p(rc, level)*delta;         adjust(N-2) = false;
//       //t(N-1) = pa.Vpp(rc)*delta*delta;  adjust(N-1) = false;
//       //t(1) = 0.0;                                    adjust(1)   = false;

// //       t(N-3) = pa.dUdiag(rc, level);     adjust(N-3) = false;
// //       t(N-2) = pa.dUdiag_p(rc, level);   adjust(N-2) = false;
// //       t(N-1) = pa.dUdiag_pp(rc, level);  adjust(N-1) = false;
// //       t(1) = 0.0;                        adjust(1)   = false;

//       // Now, do the optimal breakup:  this gives me the coefficents
//       // of the basis functions, h_n in the array t.
//       breakup.DoBreakup (Xk, t, adjust);
      
//       // Now, we must put this information into the pair action
//       // object.  First do real space part
//       pa.dUlong_r0(level)=0.0;
//       for (int n=0; n<N; n++)
// 	pa.dUlong_r0(level) += t(n)*basis.h(n,0.0);
//       for (int i=0; i<LongGrid.NumPoints; i++) {
// 	double r = LongGrid(i);
// 	if (r <= rc) {
// 	  // Sum over basis functions
// 	  for (int n=0; n<N; n++) 
// 	    dUlong_r(i) += t(n) * basis.h(n, r);
// 	}
// 	else
// 	  dUlong_r(i) = pa.dUdiag (r, level);
//       }
//       pa.dUlong(level).Init(&LongGrid, dUlong_r);

//       // Now do k-space part
//       for (int ki=0; ki < Path.kVecs.size(); ki++) {
// 	const dVec &kv = Path.kVecs(ki);
// 	double k = sqrt (dot(kv,kv));
// 	// Sum over basis functions
// 	for (int n=0; n<N; n++)
// 	  pa.dUlong_k(level,ki) += t(n) * basis.c(n,k);
// 	// Now add on part from rc to infinity
// 	pa.dUlong_k(level,ki) -= CalcXk(paIndex, level, k, rc, JOB_DU);
//       }
//     }
//   }
// }


void LongRangeClass::OptimizedBreakup_V(int numKnots)
{
  ///BUG: Optimized breakup only works when NDIM==3
#if NDIM==3
  const double tolerance = 1.0e-7;
  double kCut = Path.Getkc();
  dVec box = Path.GetBox();
  double boxVol = box[0]*box[1]*box[2];
  double rc = 0.5*box[0];
  for (int i=1; i<NDIM; i++)
    rc = min (rc, 0.5*box[i]);
  double kvol = Path.GetkBox()[0];
  for (int i=1; i<NDIM; i++)
    kvol *= Path.GetkBox()[i];
  double kavg = pow(kvol,1.0/3.0);
//   // We try to pick kcont to keep reasonable number of k-vectors
//   double kCont = 50.0 * kavg;
//   double kMax = 100 * kavg;
//   cerr << "kCont = " << kCont 
//        << " kMax = " << kMax << endl;

  LPQHI_BasisClass basis;
  basis.Set_rc(rc);
  basis.SetBox(box);
  basis.SetNumKnots (numKnots);

  // We try to pick kcont to keep reasonable number of k-vectors
  double kCont = 50.0 * kavg;
  double delta = basis.GetDelta();
  double kMax = 20.0*M_PI/delta;
  cerr << "kCont = " << kCont 
       << " kMax = " << kMax << endl;


  OptimizedBreakupClass breakup(basis);
  breakup.SetkVecs (kCut, kCont, kMax);
  int numk = breakup.kpoints.size();
  int N = basis.NumElements();
  Array<double,1> t(N);
  Array<bool,1>   adjust (N);
  Array<double,1> Xk(numk);

  // Would be 0.5, but with two timeslice distdisp, it could be a
  // little longer
  double rmax = 0.75 * sqrt (dot(box,box));
  const int numPoints = 1000;
  LongGrid.Init (0.0, rmax, numPoints);
  Array<double,1> Vlong_r(numPoints);

  for (int paIndex=0; paIndex<PairArray.size(); paIndex++) {
    PairActionFitClass &pa = *PairArray(paIndex);
    pa.Setrc (rc);
    pa.Vlong_k.resize(Path.kVecs.size());
    pa.Vlong_k = 0.0;
    Vlong_r = 0.0;

    // Calculate Xk's
    for (int ki=0; ki<numk; ki++) {
      Xk(ki) = CalcXk(paIndex, 0, breakup.kpoints(ki)[0], rc, JOB_V);
      double k = breakup.kpoints(ki)[0];
      Xk(ki) = pa.Xk_V (k) / boxVol;
    }
    
    // Set boundary conditions at rc:  Force value and first and
    // second derivatives of long-range potential to match the full
    // potential at rc.
    adjust = true;
    double delta = basis.GetDelta();
    t(N-3) = pa.V(rc);                 adjust(N-3) = false;
    t(N-2) = pa.Vp(rc)*delta;          adjust(N-2) = false;
    t(N-1) = pa.Vpp(rc)*delta*delta;   adjust(N-1) = false;
    t(1) = 0.0;                        adjust(1)   = false;
    
    // Now, do the optimal breakup:  this gives me the coefficents
    // of the basis functions, h_n in the array t.
    cerr << "Doing V breakup...\n";
    breakup.DoBreakup (Xk, t, adjust);
    cerr << "Done.\n";
    
    // Now, we must put this information into the pair action
    // object.  First do real space part
    pa.Vlong_r0=0.0;
    for (int n=0; n<N; n++)
      pa.Vlong_r0 += t(n)*basis.h(n,0.0);
    for (int i=0; i<LongGrid.NumPoints; i++) {
      double r = LongGrid(i);
      if (r <= rc) {
	// Sum over basis functions
	for (int n=0; n<N; n++) 
	  Vlong_r(i) += t(n) * basis.h(n, r);
      }
      else
	Vlong_r(i) = pa.V (r);
    }
    pa.Vlong.Init(&LongGrid, Vlong_r);
    // Calculate FT of Ushort at k=0
    UshortIntegrand integrand(pa, 0, JOB_V);
    GKIntegration<UshortIntegrand, GK31> integrator(integrand);
    integrator.SetRelativeErrorMode();
    pa.Vshort_k0 = 4.0*M_PI/boxVol * 
      integrator.Integrate(0.0, rc, tolerance);
    cerr << "Vshort_k0 = " << pa.Vshort_k0 << endl;


    // Now do k-space part
    for (int ki=0; ki < Path.kVecs.size(); ki++) {
      const dVec &kv = Path.kVecs(ki);
      double k = sqrt (dot(kv,kv));
      // Sum over basis functions
      for (int n=0; n<N; n++)
	pa.Vlong_k(ki) += t(n) * basis.c(n,k);
      // Now add on part from rc to infinity
      //pa.Vlong_k(ki) -= CalcXk(paIndex, 0, k, rc, JOB_V);
      pa.Vlong_k(ki) -= pa.Xk_V(k) / boxVol;
    }
  }
#endif
}
