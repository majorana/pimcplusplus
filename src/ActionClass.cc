#include "PathDataClass.h"

ActionClass::ActionClass(PathDataClass  &pathdata) : 
  PathData(pathdata), Path(pathdata.Path), UseRPA(false)
{
}


void ActionClass::Read(IOSectionClass& inSection)
{ 
  assert(inSection.ReadVar ("tau", tau));
  assert(inSection.ReadVar ("MaxLevels", MaxLevels));
  cerr << "MaxLevels = " << MaxLevels << endl;

  inSection.ReadVar ("UseRPA", UseRPA);
  if (UseRPA) 
    cerr << "Using RPA for long range action.\n";
  else      
    cerr << "Not using RPA for long range action.\n";


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
  bool longRange = false;
  for (int i=0; i<numPairActions; i++) {
    // cerr << "i = " << i << endl;
    assert(PAIO.OpenFile (PAFiles(i)));
    PairActionVector(i) = ReadPAFit (PAIO, tau, MaxLevels);
    longRange |= PairActionVector(i)->IsLongRange();
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

  if (longRange) {
    if (Path.Getkc() == 0.0) {
      cerr << "Missing kCutoff in System section.  Aborting.\n";
      abort();
    }
    int numKnots;
    assert(inSection.ReadVar ("NumBreakupKnots", numKnots));
    cerr << "Doing optimized long range breakups...\n";
    OptimizedBreakup_U(numKnots);
    OptimizedBreakup_dU(numKnots);
    OptimizedBreakup_V(numKnots);
    if (UseRPA) {
      cerr << "Doing RPA correction...\n";
      SetupRPA();
      cerr << "done.\n";
    }
    // Print out some debug info
    for (int i=0; i<numPairActions; i++) {
      string fname = PairActionVector(i)->Particle1.Name + "-" +
	PairActionVector(i)->Particle2.Name + ".dat";
      FILE *fout = fopen (fname.c_str(), "w");
      for (double q=0.0; q<20000.0; q+=1.0) {
	double Udiag = PairActionVector(i)->Udiag(q, 0);
	double Ulong = PairActionVector(i)->Ulong(0)(q);
	double dUdiag = PairActionVector(i)->dUdiag(q, 0);
	double dUlong = PairActionVector(i)->dUlong(0)(q);
	double V = PairActionVector(i)->V(q);
	double Vlong = PairActionVector(i)->Vlong(q);
	fprintf (fout, "%1.16e %1.16e %1.16e %1.16e %1.16e %1.16e %1.16e\n", 
		 q, Udiag, Ulong, dUdiag, dUlong, V, Vlong);
      }
      fclose (fout);
    }
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

double ActionClass::OtherAction(int startSlice, int endSlice, 
		   const Array<int,1> &changedParticles, int level)
{

  int openLink=PathData.Path.OpenLink;
  int openPtcl=PathData.Path.OpenPtcl;
  dVec disp;
  double dist;
  PathData.Path.DistDisp(openLink,openPtcl,PathData.Path.NumParticles(),
			 dist,disp); //This is distance between head and tail!
  return log(0.1+(dist*dist)*(0.94*exp(-dist*dist)+0.06));


}



double ActionClass::UApproximateAction (int startSlice, int endSlice, 
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
    for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++) {
      if (Path.DoPtcl(ptcl2)){
	int PairIndex = PairMatrix(species1,
				   Path.ParticleSpeciesNum(ptcl2));

	for (int slice=startSlice;slice<endSlice;slice+=skip){
	  dVec r, rp;
	  double rmag, rpmag;
	  double U;
	  PathData.Path.DistDisp(slice, slice+skip, ptcl1, ptcl2,
				 rmag, rpmag, r, rp);
	  
	  if (level>0){
	    double V1=PairActionVector(PairIndex)->V(rmag);
	    double V2=PairActionVector(PairIndex)->V(rpmag);
	    U=(PairActionVector(PairIndex)->U(rmag,0,0,level)+
	       PairActionVector(PairIndex)->U(rpmag,0,0,level))/2.0;
	    

	  }
	  else{
	    double s2 = dot (r-rp, r-rp);
	    double q = 0.5 * (rmag + rpmag);
	    double z = (rmag - rpmag);


	    U = PairActionVector(PairIndex)->U(q,z,s2, level);
	    // Subtract off long-range part from short-range action
	    if (PairActionVector(PairIndex)->IsLongRange())
	      U -= 0.5* (PairActionVector(PairIndex)->Ulong(level)(rmag) +
			 PairActionVector(PairIndex)->Ulong(level)(rpmag));

	  }
	  TotalU += U;
	}
      }
    }
  }
  if (PathData.Path.LongRange){
    // Now add in the long-range part of the action
    // In primitive form, end slices get weighted by 1/2.  Others by 1.
    for (int slice=startSlice;slice<=endSlice;slice+=skip) 
      if ((slice==startSlice) || (slice == endSlice))
	TotalU += 0.5*LongRange_U (slice, level);
      else
	TotalU += LongRange_U (slice, level);
  }
  return (TotalU);
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
    for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++) {
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
	  // Subtract off long-range part from short-range action
	  if (PairActionVector(PairIndex)->IsLongRange())
	    U -= 0.5* (PairActionVector(PairIndex)->Ulong(level)(rmag) +
		       PairActionVector(PairIndex)->Ulong(level)(rpmag));
	  TotalU += U;
	}
      }
    }
  }

  if (PathData.Path.LongRange){
    // Now add in the long-range part of the action
    // In primitive form, end slices get weighted by 1/2.  Others by 1.
    for (int slice=startSlice;slice<=endSlice;slice+=skip) 
      if ((slice==startSlice) || (slice == endSlice))
	TotalU += 0.5*LongRange_U (slice, level);
      else
	TotalU += LongRange_U (slice, level);
  }
  return (TotalU);

  if (PathData.Path.LongRange){
    // Now add in the long-range part of the action
    // In primitive form, end slices get weighted by 1/2.  Others by 1.
    if (UseRPA) {
      for (int slice=startSlice;slice<=endSlice;slice+=skip) 
	if ((slice==startSlice) || (slice == endSlice))
	  TotalU += 0.5*LongRange_U_RPA (slice, level);
	else
	  TotalU += LongRange_U_RPA (slice, level);
    }
    else {
      for (int slice=startSlice;slice<=endSlice;slice+=skip) 
	if ((slice==startSlice) || (slice == endSlice))
	  TotalU += 0.5*LongRange_U (slice, level);
	else
	  TotalU += LongRange_U (slice, level);
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

  if (PathData.Path.OpenPaths){
    return UApproximateAction(startSlice,endSlice,changedParticles,level)+
      KAction(startSlice,endSlice,changedParticles,level)+
      OtherAction(startSlice,endSlice,changedParticles,level);
  }
  else {
    return UAction(startSlice,endSlice,changedParticles,level)+
      KAction(startSlice,endSlice,changedParticles,level);
  }

}


  

inline double mag2 (const complex<double> &z)
{
  return (z.real()*z.real() + z.imag()*z.imag());
}


/// Calculates the long-range part of the action at a given timeslice  
double ActionClass::LongRange_U(int slice, int level)
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
    int N = Path.Species(species).NumParticles;
    // We can't forget the Madelung term.
    homo -= 0.5 * N * PA.Ulong_r0(level);
    // Or the neutralizing background term
    homo -= 0.5*N*N*PA.Ushort_k0(level);
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
	int N1 = Path.Species(species1).NumParticles;
	int N2 = Path.Species(species2).NumParticles;
	hetero -= N1*N2*PA.Ushort_k0(level);
      }
    }
//   cerr << "homo = " << homo << endl;
//   cerr << "hetero = " << hetero << endl;
  return (homo+hetero);
}

/// Calculates the long-range part of the action at a given timeslice  
double ActionClass::LongRange_dU(int slice, int level)
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
	homo += 0.5 * 2.0* rhok2 * PA.dUlong_k(level,ki);
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
      int paIndex = PairMatrix(species1, species2);
      PairActionFitClass &PA = *PairActionVector(paIndex);
      if (PA.IsLongRange()) {
	for (int ki=0; ki<Path.kVecs.size(); ki++) {
	  double rhorho = 
	    Path.Rho_k(slice, species1, ki).real() *
	    Path.Rho_k(slice, species2, ki).real() + 
	    Path.Rho_k(slice, species1, ki).imag() *
	    Path.Rho_k(slice, species2, ki).imag();
	  hetero += 2.0 * rhorho * PA.dUlong_k(level,ki);
	}
	int N1 = Path.Species(species1).NumParticles;
	int N2 = Path.Species(species2).NumParticles;
	hetero -= N1*N2*PA.dUshort_k0(level);
      }
    }
  return (homo+hetero);
}

////////////////////
/// RPA Versions ///
////////////////////
/// Calculates the long-range part of the action at a given timeslice  
double ActionClass::LongRange_U_RPA(int slice, int level)
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
	homo += 0.5 * 2.0* rhok2 * PA.U_RPA_long_k(level,ki);
      }
    }
    int N = Path.Species(species).NumParticles;
    // We can't forget the Madelung term.
    homo -= 0.5 * N * PA.Ulong_r0(level);
    // Or the neutralizing background term
    homo -= 0.5*N*N*PA.Ushort_k0(level);
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
	  hetero += 2.0 * rhorho * PA.U_RPA_long_k(level,ki);
	}
	int N1 = Path.Species(species1).NumParticles;
	int N2 = Path.Species(species2).NumParticles;
	hetero -= N1*N2*PA.Ushort_k0(level);
      }
    }
//   cerr << "homo = " << homo << endl;
//   cerr << "hetero = " << hetero << endl;
  return (homo+hetero);
}

/// Calculates the long-range part of the action at a given timeslice  
double ActionClass::LongRange_dU_RPA(int slice, int level)
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
	homo += 0.5 * 2.0* rhok2 * PA.dU_RPA_long_k(level,ki);
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
      int paIndex = PairMatrix(species1, species2);
      PairActionFitClass &PA = *PairActionVector(paIndex);
      if (PA.IsLongRange()) {
	for (int ki=0; ki<Path.kVecs.size(); ki++) {
	  double rhorho = 
	    Path.Rho_k(slice, species1, ki).real() *
	    Path.Rho_k(slice, species2, ki).real() + 
	    Path.Rho_k(slice, species1, ki).imag() *
	    Path.Rho_k(slice, species2, ki).imag();
	  hetero += 2.0 * rhorho * PA.dU_RPA_long_k(level,ki);
	}
	int N1 = Path.Species(species1).NumParticles;
	int N2 = Path.Species(species2).NumParticles;
	hetero -= N1*N2*PA.dUshort_k0(level);
      }
    }
  return (homo+hetero);
}





/// Calculates the long-range part of the action at a given timeslice  
double ActionClass::LongRange_V(int slice)
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
	homo += 0.5 * 2.0* rhok2 * PA.Vlong_k(ki);
      }
    }
    int N = Path.Species(species).NumParticles;
    // We can't forget the Madelung term.
    homo -= 0.5 * N * PA.Vlong_r0;
    // Or the neutralizing background term
    homo -= 0.5*N*N*PA.Vshort_k0;
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
	  hetero += 2.0 * rhorho * PA.Vlong_k(ki);
	}
	int N1 = Path.Species(species1).NumParticles;
	int N2 = Path.Species(species2).NumParticles;
	hetero -= N1*N2*PA.Vshort_k0;
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
  TaskType Task;

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
    if (Task == DO_U)
      return Uintegrand(r);
    else if (Task == DO_DU)
      return dUintegrand(r);
    else
      return Vintegrand(r);
  }

  CoulombXkIntegrand (PairActionFitClass &pa, int level, double k_,
		      TaskType task) :
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
  TaskType Task;
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
    if (Task == DO_U)
      return Uintegrand(r);
    else if (Task == DO_DU)
      return dUintegrand(r);
    else
      return Vintegrand(r);
 
  }
  XkIntegrand (PairActionFitClass &pa, int level, double k_,
	       TaskType task) :
    PA(pa), Level(level), k(k_), Task(task)
  { /* do nothing else*/  }
};



/// This calculates the quantity 
/// \f$ X_k \equiv -\frac{4 \pi}{\Omega k} \int_{r_c}^\infty dr \, r \sin(kr) V(r).\f$
double ActionClass::CalcXk (int paIndex, int level, double k, double rc,
			    TaskType task)
{
  double absTol = 1.0e-7;
  double relTol = 1.0e-5;

  PairActionFitClass &pa = *PairActionVector(paIndex);
  if (task == DO_U)
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
    
    
    if (false/*task != DO_V*/) {
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
    if (task == DO_U) {
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
  TaskType Task;
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
    if (Task == DO_U)
      return Uintegrand(r);
    else if (Task == DO_DU)
      return dUintegrand(r);
    else
      return Vintegrand(r);
 
  }
  UshortIntegrand (PairActionFitClass &pa, int level,
	       TaskType task) :
    PA(pa), Level(level), Task(task)
  { /* do nothing else*/  }
};


/// This computes the optimized breakups for the pair actions stored
/// in PairActionVector.  The parameters are the number of knots in
/// the "spline" representation of the long-range action and the
/// k-space cutoff.  
/// Only \f$\mathbf{k}\f$ with \f$|\mathbf{k}| < k_c$\f will be
/// included in the simulation sum.
void ActionClass::OptimizedBreakup_U(int numKnots)
{
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

  for (int paIndex=0; paIndex<PairActionVector.size(); paIndex++) {
    PairActionFitClass &pa = *PairActionVector(paIndex);
    pa.Setrc (rc);
    pa.Ulong.resize(MaxLevels);
    pa.Ulong_k.resize(MaxLevels,Path.kVecs.size());
    pa.Ulong_k = 0.0;
    pa.Ulong_r0.resize(MaxLevels);
    pa.Ushort_k0.resize(MaxLevels);
    for (int level=0; level<MaxLevels; level++) {
      Ulong_r = 0.0;
      
      // Calculate Xk's
      cerr << "Calculating Xk's for U...\n";
      for (int ki=0; ki<numk; ki++) {
	//Xk(ki) = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc, DO_U);
	//double oldXk = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc, DO_U);
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
      UshortIntegrand integrand(pa, level, DO_U);
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
	// pa.Ulong_k(level,ki) -= CalcXk(paIndex, level, k, rc, DO_U);
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
}

void ActionClass::OptimizedBreakup_dU(int numKnots)
{
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

  for (int paIndex=0; paIndex<PairActionVector.size(); paIndex++) {
    PairActionFitClass &pa = *PairActionVector(paIndex);
    pa.Setrc (rc);
    pa.dUlong.resize(MaxLevels);
    pa.dUlong_k.resize(MaxLevels,Path.kVecs.size());
    pa.dUlong_k = 0.0;
    pa.dUlong_r0.resize(MaxLevels);
    pa.dUshort_k0.resize(MaxLevels);
    for (int level=0; level<MaxLevels; level++) {
      dUlong_r = 0.0;

      // Calculate Xk's
      cerr << "Calculating Xk's for dU...\n";
      for (int ki=0; ki<numk; ki++) {	
	// Xk(ki) = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc, DO_DU);
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
      UshortIntegrand integrand(pa, level, DO_DU);
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
	//pa.dUlong_k(level,ki) -= CalcXk(paIndex, level, k, rc, DO_DU);
	pa.dUlong_k(level,ki) -= pa.Xk_dU(k, level) / boxVol;
      }
    }
  }
}



// void ActionClass::OptimizedBreakup_dU(int numKnots)
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

//   for (int paIndex=0; paIndex<PairActionVector.size(); paIndex++) {
//     PairActionFitClass &pa = *PairActionVector(paIndex);
//     pa.dUlong.resize(MaxLevels);
//     pa.dUlong_k.resize(MaxLevels,Path.kVecs.size());
//     pa.dUlong_k = 0.0;
//     pa.dUlong_r0.resize(MaxLevels);
//     for (int level=0; level<MaxLevels; level++) {
//       dUlong_r = 0.0;
//       // Calculate Xk's
//       for (int ki=0; ki<numk; ki++)
// 	Xk(ki) = CalcXk(paIndex, level, breakup.kpoints(ki)[0], rc,
// 			DO_DU);

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
// 	pa.dUlong_k(level,ki) -= CalcXk(paIndex, level, k, rc, DO_DU);
//       }
//     }
//   }
// }


void ActionClass::OptimizedBreakup_V(int numKnots)
{
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

  for (int paIndex=0; paIndex<PairActionVector.size(); paIndex++) {
    PairActionFitClass &pa = *PairActionVector(paIndex);
    pa.Setrc (rc);
    pa.Vlong_k.resize(Path.kVecs.size());
    pa.Vlong_k = 0.0;
    Vlong_r = 0.0;

    // Calculate Xk's
    for (int ki=0; ki<numk; ki++) {
      Xk(ki) = CalcXk(paIndex, 0, breakup.kpoints(ki)[0], rc, DO_V);
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
    UshortIntegrand integrand(pa, 0, DO_V);
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
      //pa.Vlong_k(ki) -= CalcXk(paIndex, 0, k, rc, DO_V);
      pa.Vlong_k(ki) -= pa.Xk_V(k) / boxVol;
    }
  }
}



/// Make sure you move the join to the end so we don't have to worry
/// about permutations before calling this.
void ActionClass::Energy(int slice1, int level,
			 double &spring, double &dU)
{
  double levelTau=tau;
  int numPtcls = PathData.NumParticles();
  int slice2 = slice1 + (1<<level);
  for (int i=0; i<level; i++) 
    levelTau *= 2.0;
  // Add constant part.  Note: we should really check the number of
  // dimensions. 
  spring = dU = 0.0;
  const int NumImage=1;
  for (int ptcl=0; ptcl<numPtcls; ptcl++)
    if (PathData.Path.ParticleSpecies(ptcl).lambda != 0.0)
      spring += 1.5/levelTau;

  for (int ptcl1=0; ptcl1<numPtcls; ptcl1++) {
    // Do free-particle part
    int species1 = PathData.Path.ParticleSpeciesNum(ptcl1);
    double lambda = PathData.Path.ParticleSpecies(ptcl1).lambda;
    if (lambda != 0.0) {
      double FourLambdaTauInv = 
	1.0/(4.0*PathData.Path.Species(species1).lambda*levelTau);
      dVec vel;
      vel = PathData.Path.Velocity(slice1, slice2, ptcl1);
      double Z = 1.0;
      dVec GaussSum=0.0;
      for (int dim=0; dim<NDIM; dim++) {
	for (int image=-NumImage; image<=NumImage; image++) {
	  double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	    GaussSum[dim] += exp(-dist*dist*FourLambdaTauInv);
	}
	Z *= GaussSum[dim];
      }
      dVec numSum=0.0;
      for (int dim=0;dim<NDIM;dim++){
	for (int image=-NumImage;image<=NumImage;image++){
	  double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	  numSum[dim] += 
	    (-dist*dist*FourLambdaTauInv/levelTau)*exp(-dist*dist*FourLambdaTauInv);
	}
      }
      double scalarnumSum=0.0;
      for (int dim=0;dim<NDIM;dim++){
	dVec numProd=1.0;
	for (int dim2=0;dim2<NDIM;dim2++){
	  if (dim2!=dim){
	    numProd[dim] *= GaussSum[dim2];
	  }
	  else {
	    numProd[dim] *=  numSum[dim2];
	  }
	  
	}
	scalarnumSum += numProd[dim];
      }
      spring += scalarnumSum/Z; 
    }
    
    
    for (int ptcl2=0; ptcl2<ptcl1; ptcl2++) {
      dVec r, rp;
      double rmag, rpmag;
      PathData.Path.DistDisp(slice1, slice2, ptcl1, ptcl2, rmag,rpmag,r,rp); 
      
      double s2 = dot(r-rp, r-rp);
      double q = 0.5*(rmag+rpmag);
      double z = (rmag-rpmag);

      int PairIndex = 
	PairMatrix(species1, 
				   PathData.Path.ParticleSpeciesNum(ptcl2));
      dU += PairActionVector(PairIndex)->dU(q, z, s2, level);
      // Subtract off long-range part from short-range action
      if (PairActionVector(PairIndex)->IsLongRange())
	dU -= 0.5*(PairActionVector(PairIndex)->dUlong(level)(rmag)+
		   PairActionVector(PairIndex)->dUlong(level)(rpmag));
		   
    }
  }
  if (PathData.Path.LongRange){
    // Add long range part of energy
    dU += 0.5*(LongRange_dU (slice1, level)+LongRange_dU(slice2,level));
  }
}




double ActionClass::PotentialEnergy (int slice)
{
  double vSum=0.0;
  for (int ptcl1=0; ptcl1<Path.NumParticles(); ptcl1++){
    int species1=PathData.Path.ParticleSpeciesNum(ptcl1);
    for (int ptcl2=0;ptcl2<ptcl1;ptcl2++){
      int species2=PathData.Path.ParticleSpeciesNum(ptcl2);
      int PairIndex =PairMatrix(species1,species2);
      dVec r;
      double rmag;
      PathData.Path.DistDisp(slice, ptcl1, ptcl2,rmag, r); 
      vSum +=(PairActionVector(PairIndex))->V(rmag);
      if (PairActionVector(PairIndex)->IsLongRange())
	vSum -= PairActionVector(PairIndex)->Vlong(rmag);
    }
  }
  if (PathData.Path.LongRange){
    vSum += LongRange_V(slice);
  }
  return vSum;
} 
 
