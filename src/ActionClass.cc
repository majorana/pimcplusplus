#include "PathDataClass.h"

ActionClass::ActionClass(PathDataClass  &pathdata) : 
  PathData(pathdata), Path(pathdata.Path)
{
}


void ActionClass::Read(IOSectionClass& inSection)
{ 
  assert(inSection.ReadVar ("tau", tau));
  assert(inSection.ReadVar ("MaxLevels", MaxLevels));
  cerr << "MaxLevels = " << MaxLevels << endl;

  Array<string,1> PAFiles;
  assert (inSection.ReadVar ("PairActionFiles", PAFiles));
  int numPairActions = PAFiles.size();
  PairActionVector.resize(numPairActions);
  PairMatrix.resize(Path.NumSpecies(),Path.NumSpecies());
  // Initialize to a nonsense value so we can later check in the table
  // element was filled in.

  PairMatrix = -1;
  // Read pair actions files
  IOSectionClass PAIO;
  for (int i=0; i<numPairActions; i++) {
    // cerr << "i = " << i << endl;
    assert(PAIO.OpenFile (PAFiles(i)));
    PairActionVector(i) = ReadPAFit (PAIO, tau, MaxLevels);
    int type1 = Path.SpeciesNum(PairActionVector(i)->Particle1.Name);
    int type2 = Path.SpeciesNum(PairActionVector(i)->Particle2.Name);
    if (type1==-1) {
      cerr << "Unrecognized type \""
	   << PairActionVector(i)->Particle1.Name << "\".\n";
      abort();
    }
    if (type2==-1) {
      cerr << "Unrecognized type \""
	   << PairActionVector(i)->Particle2.Name << "\".\n";
      abort();
    }
    PairMatrix(type1,type2) = i;
    PairMatrix(type2,type1) = i;
    PAIO.CloseFile();
  }



  // Now check to make sure all PairActions that we need are defined.
  for (int species1=0; species1<Path.NumSpecies(); species1++)
    for (int species2=0; species2<Path.NumSpecies(); species2++)
      if (PairMatrix(species1,species2) == -1) {
	if ((species1 != species2) || 
	    (Path.Species(species1).NumParticles > 1)) {
	  cerr << "We're missing a PairAction for species1 = "
	       << Path.Species(species1).Name << " and species2 = "
	       << Path.Species(species2).Name << endl;
	  exit(1);
	}
      }
  cerr << "Finished reading the action.\n"; 
}


double ActionClass::UAction (int startSlice, int endSlice, 
			     const Array<int,1> &changedParticles, int level)
{
  // First, sum the pair actions
  for (int counter=0;counter<Path.DoPtcl.size();counter++){
    Path.DoPtcl(counter)=true;
  }
  double TotalU = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = tau* (1<<level);
  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = changedParticles(ptcl1Index);
    Path.DoPtcl(ptcl1) = false;
    int species1=Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++){
      if (Path.DoPtcl(ptcl2)){
	int PairIndex = PairMatrix(species1,
				   Path.ParticleSpeciesNum(ptcl2));

	for (int slice=startSlice;slice<endSlice;slice+=skip){
	  dVec r, rp;
	  double rmag, rpmag;

	  PathData.Path.DistDisp(slice, slice+skip, ptcl1, ptcl2,
				 rmag, rpmag, r, rp);

	  double s2 = dot (r-rp, r-rp);
	  double q = 0.5 * (rmag + rpmag);
	  double z = (rmag - rpmag);

	  double U;
	  U = PairActionVector(PairIndex)->U(q,z,s2, level);
	  U -= PairActionVector(PairIndex)->Ulong(level)(q);
	  TotalU += U;
	}
      }
    }
  }
  return (TotalU);
}

double ActionClass::KAction (int startSlice, int endSlice, 
			     const Array<int,1> &changedParticles, int level)
{
  double TotalK = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = tau* (1<<level);

  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl = changedParticles(ptclIndex);
    int species=Path.ParticleSpeciesNum(ptcl);
    double FourLambdaTauInv=1.0/(4.0*Path.Species(species).lambda*levelTau);
    for (int slice=startSlice; slice < endSlice;slice+=skip) {
      dVec vel;
      vel = PathData.Path.Velocity(slice, slice+skip, ptcl);
      double GaussProd = 1.0;
      for (int dim=0; dim<NDIM; dim++) {
	int NumImage=1;
	double GaussSum=0.0;
	for (int image=-NumImage; image<=NumImage; image++) {
	  double dist = vel[dim]+(double)image*Path.GetBox()[dim];
	  GaussSum += exp(-dist*dist*FourLambdaTauInv);
	}
	GaussProd *= GaussSum;
      }
      TotalK -= log(GaussProd);    
      //TotalK += dot(vel,vel)*FourLambdaTauInv; 
    }
  }
  //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
  return (TotalK);
}


double ActionClass::TotalAction(int startSlice, int endSlice, 
				const Array<int,1> &changedParticles,
				int level)
{
  return UAction(startSlice,endSlice,changedParticles,level)+
    KAction(startSlice,endSlice,changedParticles,level);
}


  

inline double mag2 (const complex<double> &z)
{
  return (z.real()*z.real() + z.imag()*z.imag());
}


/// Calculates the long-range part of the action at a given timeslice  
double ActionClass::LongRangeAction(int slice, int level)
{
  double homo = 0.0;
  double hetero = 0.0;

  // First, do the homologous (same species) terms
  for (int species=0; species<Path.NumSpecies(); species++) {
    Path.CalcRho_ks_Fast(slice,species);
    int paIndex = PairMatrix(species,species);
    PairActionFitClass &PA = *PairActionVector(paIndex);
    if (PA.IsLongRange()) {
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
	double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	homo += 0.5 * 2.0* rhok2 * PA.Ulong_k(level,ki);
      }
    }
    // We can't forget the Madelung term.
    homo -= 0.5 * Path.Species(species).NumParticles * PA.Ulong_0(level);
  }

  // Now do the heterologous terms
  for (int species1=0; species1<Path.NumSpecies(); species1++)
    for (int species2=species1+1; species2<Path.NumSpecies(); species2++) {
      int paIndex = PairMatrix(species1, species2);
      PairActionFitClass &PA = *PairActionVector(paIndex);
      if (PA.IsLongRange()) {
	for (int ki=0; ki<Path.kVecs.size(); ki++) {
	  double rhorho = 
	    Path.Rho_k(slice, species1, ki).real() *
	    Path.Rho_k(slice, species2, ki).real() + 
	    Path.Rho_k(slice, species1, ki).imag() *
	    Path.Rho_k(slice, species2, ki).imag();
	  hetero += 2.0 * rhorho * PA.Ulong_k(level,ki);
	}
      }
    }
  return (homo+hetero);
}


#include "Common/Ewald/OptimizedBreakup.h"
#include "Common/Integration/GKIntegration.h"

class CoulombXkIntegrand
{
private:
  PairActionFitClass &PA;
  int Level;
  double k;
  double beta;
public:
  inline double operator()(double r)
  {
    double U = PA.Udiag(r, Level);
    U -= beta * PA.Z1Z2/r;
    return r * sin(k*r)*U;
  }
  CoulombXkIntegrand (PairActionFitClass &pa, int level, double k_) :
    PA(pa), Level(level), k(k_)
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
public:
  inline double operator()(double r) 
  {
    return r * sin(k*r) * PA.Udiag(r, Level);
  }
  XkIntegrand (PairActionFitClass &pa, int level, double k_) :
    PA(pa), Level(level), k(k_)
  { /* do nothing else*/  }
};



/// This calculates the quantity 
/// \f$ X_k \equiv -\frac{4 \pi}{\Omega k} \int_{r_c}^\infty dr \, r \sin(kr) V(r).\f$
double ActionClass::CalcXk (int paIndex, int level, double k, double rc)
{
  const double tolerance = 1.0e-10;
  PairActionFitClass &pa = *PairActionVector(paIndex);
  if (pa.Z1Z2 == 0.0) {
    XkIntegrand integrand(pa, level, k);
    GKIntegration<XkIntegrand, GK31> integrator(integrand);
    integrator.SetRelativeErrorMode();
    double Xk = -4.0*M_PI/(Path.GetVol()*k) * 
      integrator.Integrate(rc, 30.0*rc, tolerance,tolerance,false);
    return Xk;
  }
  else {
    CoulombXkIntegrand integrand(pa, level, k);
    GKIntegration<CoulombXkIntegrand, GK31> integrator(integrand);
    integrator.SetRelativeErrorMode();
    double Xk = -4.0*M_PI/(Path.GetVol()*k) * 
      integrator.Integrate(rc, 30.0*rc, tolerance,tolerance,false);
    /// Add in the analytic part that I ignored
    double beta = pa.SmallestBeta;
    for (int i=0; i<level; i++)
      beta *= 2.0;
    //Xk -= beta*4.0*M_PI*pa.Z1Z2/(Path.GetVol()*k*k)*cos(k*rc);
    Xk = -beta*4.0*M_PI*pa.Z1Z2/(Path.GetVol()*k*k)*cos(k*rc);
    return Xk;
  }
  
}


void ActionClass::OptimizedBreakup(int numKnots, double kCut)
{
  dVec box = Path.GetBox();
  double rc = 0.5*box[0];
  for (int i=1; i<NDIM; i++)
    rc = min (rc, 0.5*box[i]);
  double kvol = Path.GetkBox()[0];
  for (int i=1; i<NDIM; i++)
    kvol *= Path.GetkBox()[i];
  double kavg = pow(kvol,1.0/3.0);
  /// We try to pick kcont to keep reasonable number of k-vectors
  double kCont = 40.0 * kavg;
  double kMax = 800 * kavg;

  LPQHI_BasisClass basis;
  basis.Set_rc(rc);
  basis.SetBox(box);
  basis.SetNumKnots (numKnots);

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
  UlongGrid.Init (0.0, rmax, 1000);
  Array<double,1> Ulong_r(1000);

  for (int paIndex=0; paIndex<PairActionVector.size(); paIndex++) {
    PairActionFitClass &pa = *PairActionVector(paIndex);
    pa.Ulong.resize(MaxLevels);
    pa.Ulong_k.resize(MaxLevels,Path.kVecs.size());
    pa.Ulong_k = 0.0;
    pa.Ulong_0.resize(MaxLevels);
    for (int level=0; level<MaxLevels; level++) {
      Ulong_r = 0.0;
      /// Calculate Xk's
      for (int ki=0; ki<numk; ki++)
	Xk(ki) = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc);

      /// Set boundary conditions at rc:  Force value and first and
      /// second derivatives of long-range potential to match the full
      /// potential at rc.
      adjust = true;
      t(N-3) = pa.Udiag(rc, level);     adjust(N-3) = false;
      t(N-2) = pa.Udiag_p(rc, level);   adjust(N-2) = false;
      t(N-1) = pa.Udiag_pp(rc, level);  adjust(N-1) = false;
      t(1) = 0.0;                       adjust(1)   = false;

      /// Now, do the optimal breakup:  this gives me the coefficents
      /// of the basis functions, h_n in the array t.
      breakup.DoBreakup (Xk, t, adjust);
      
      cerr << "t = " << t << endl;
      /// Now, we must put this information into the pair action
      /// object.  First do real space part
      pa.Ulong_0(level)=0.0;
      for (int n=0; n<N; n++)
	pa.Ulong_0(level) += t(n)*basis.h(n,0.0);
      for (int i=0; i<UlongGrid.NumPoints; i++) {
	double r = UlongGrid(i);
	if (r <= rc) {
	  /// Sum over basis functions
	  for (int n=0; n<N; n++) 
	    Ulong_r(i) += t(n) * basis.h(n, r);
	}
	else
	  Ulong_r(i) = pa.Udiag (r, level);
      }
      pa.Ulong(level).Init(&UlongGrid, Ulong_r);

      // Now do k-space part
      for (int ki=0; ki < Path.kVecs.size(); ki++) {
	const dVec &kv = Path.kVecs(ki);
	double k = sqrt (dot(kv,kv));
	/// Sum over basis functions
	for (int n=0; n<N; n++)
	  pa.Ulong_k(level,ki) += t(n) * basis.c(n,k);
	/// Now add on part from rc to infinity
	pa.Ulong_k(level,ki) -= CalcXk(paIndex, level, k, rc);
      }
    }
  }
}
