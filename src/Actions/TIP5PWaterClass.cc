#include "TIP5PWaterClass.h"
#include "../PathDataClass.h"

TIP5PWaterClass::TIP5PWaterClass (PathDataClass &pathData) :
  ActionBaseClass (pathData)
{
  //Do  nothing for now
}

double TIP5PWaterClass::Action (int startSlice, int endSlice, 
	       const Array<int,1> &activeParticles, int level)
{
//  cerr << "I'm calculating the Action" << endl;
  for (int counter=0;counter<Path.DoPtcl.size();counter++){
    Path.DoPtcl(counter)=true;
  }

// obnoxious physical constants
  double elementary_charge = 1.602*pow(10.0,-19);
  double N_Avogadro = 6.022*pow(10.0,23.0);
  double kcal_to_joule = 4184;
  double epsilon_not = 8.85*pow(10.0,-12);
  double angstrom_to_m = pow(10.0,10);
  double SI = 1/(4*M_PI*epsilon_not);
  double k_B = 1.3807*pow(10.0,-23);
  double erg_to_eV = 1.6*pow(10.0,12);
  double joule_to_eV = pow(10.0,19)/1.602;

  double TotalU = 0.0;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;
  int speciesO=Path.SpeciesNum("O");
  int speciesp=Path.SpeciesNum("p");
  int speciese=Path.SpeciesNum("e");
  double CUTOFF = 7.0; // setcutoff: spherical cutoff in angstroms

  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = activeParticles(ptcl1Index);
/*   need to exclude intramolecular interactions - compare with line commented out below
      for(int i = 0;i<=4;++i){
        Path.DoPtcl(ptcl1+4*i) = false;
      }
*/
    Path.DoPtcl(ptcl1) = false;
    int species1=Path.ParticleSpeciesNum(ptcl1);
//    cerr << "Choosing particle " << ptcl1 << " which is of species " << species1 << endl;

    if (species1==speciesO){
      int counter = 0;
//    Calculate o-o interaction
      for (int ptcl2=Path.Species(speciesO).FirstPtcl;ptcl2<=Path.Species(speciesO).LastPtcl;ptcl2++) {///loop over oxygen
	if (Path.DoPtcl(ptcl2)){
//          int PairIndex = PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2));
	  for (int slice=startSlice;slice<=endSlice;slice+=skip){
	    dVec r, rp;
	    double rmag, rpmag;
	    PathData.Path.DistDisp(slice, ptcl1, ptcl2, rmag, r);
//  implement spherical cutoff  
            if (rmag <= CUTOFF){
              counter ++;
              double sigma_over_r = PathData.Species(species1).Sigma/rmag;
              double sigma_over_cutoff = PathData.Species(species1).Sigma/CUTOFF;
              double offset = pow(sigma_over_cutoff,12) - pow(sigma_over_cutoff,6);
	      double lj = 4*PathData.Species(species1).Epsilon*(pow(sigma_over_r,12)-pow(sigma_over_r,6) - offset); // this is in kcal/mol 
	      TotalU += lj;
//	      cerr << "lj  " << lj << " at distance " << rmag << endl;
            }
            else{
//              cerr << rmag << " is outside cutoff" << endl;
            }
	  }
	}
      }
      ///end calculating o-o interactions
//      cerr << "I calculated Lennard-Jones interactions with: " << counter << " particles" << endl;
    }
    else{
      /// calculating coulomb interactions
      for (int ptcl2=Path.Species(speciesp).FirstPtcl;ptcl2<=Path.Species(speciesp).LastPtcl;ptcl2++) {
        ///loop over protons
//  don't compute intramolecular interactions
	if (Path.DoPtcl(ptcl2)&&Path.MolRef(ptcl1)!=Path.MolRef(ptcl2)){
	  for (int slice=startSlice;slice<=endSlice;slice+=skip){
	    double rmag;
	    dVec r;
	    PathData.Path.DistDisp(slice,ptcl1,ptcl2,rmag,r);
// implement spherical cutoff
            double Ormag = OOSeparation(slice,ptcl1,ptcl2);
            if (Ormag <= CUTOFF){
	      double coulomb_const = SI*angstrom_to_m*elementary_charge*elementary_charge*PathData.Species(species1).Charge*PathData.Species(speciesp).Charge*N_Avogadro/kcal_to_joule;
              double coulomb = coulomb_const*(1.0/rmag - 1.0/CUTOFF);
	      TotalU += coulomb;
//	      cerr << "protons " << coulomb << " at distance " << rmag  << " with offset " << coulomb_const/CUTOFF << endl;
            }
          }
	}
      }
      for (int ptcl2=Path.Species(speciese).FirstPtcl;ptcl2<=Path.Species(speciese).LastPtcl;ptcl2++) {
        ///loop over electrons
//  don't compute intramolecular interactions
	if (Path.DoPtcl(ptcl2)&&Path.MolRef(ptcl1)!=Path.MolRef(ptcl2)){
	  for (int slice=startSlice;slice<=endSlice;slice+=skip){
	    double rmag;
	    dVec r;
	    PathData.Path.DistDisp(slice,ptcl1,ptcl2,rmag,r);
// implement spherical cutoff
            double Ormag = OOSeparation(slice,ptcl1,ptcl2);
            if (Ormag <= CUTOFF){
	      double coulomb_const = SI*angstrom_to_m*elementary_charge*elementary_charge*PathData.Species(species1).Charge*PathData.Species(speciese).Charge*N_Avogadro/kcal_to_joule;
              double coulomb = coulomb_const*(1.0/rmag - 1.0/CUTOFF);
	      TotalU += coulomb;
//	      cerr << "electrons " << coulomb << " at distance " << rmag << " with offset " << coulomb_const/CUTOFF << endl;
            }
	  }
	}
      }
    
	
    }
   /// end calculating coulomb interactions 
  
  
  }
  double TotalU_times_tau = TotalU*PathData.Path.tau;
//  cerr << TotalU << " and times tau " << TotalU_times_tau << " at temp " << 1.0/PathData.Path.tau << endl;
//  cerr << "I'm returning TIP5P action " << TotalU_times_tau << endl;
  return (TotalU_times_tau);
}

double TIP5PWaterClass::d_dBeta (int startSlice, int endSlice,  int level)
{
  double thermal = 6/PathData.Path.tau;
  Array<int,1> activeParticles(PathData.Path.NumParticles());
  for (int i=0;i<PathData.Path.NumParticles();i++){
    activeParticles(i)=i;
  }
 
  for (int counter=0;counter<Path.DoPtcl.size();counter++){
    Path.DoPtcl(counter)=true;
  }

// obnoxious physical constants
  double elementary_charge = 1.602*pow(10.0,-19);
  double N_Avogadro = 6.022*pow(10.0,23.0);
  double kcal_to_joule = 4184;
  double epsilon_not = 8.85*pow(10.0,-12);
  double angstrom_to_m = pow(10.0,10);
  double SI = 1/(4*M_PI*epsilon_not);
  double k_B = 1.3807*pow(10.0,-23);
  double erg_to_eV = 1.6*pow(10.0,12);
  double joule_to_eV = pow(10.0,19)/1.602;

  double TotalU = 0.0;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;
  int speciesO=Path.SpeciesNum("O");
  int speciesp=Path.SpeciesNum("p");
  int speciese=Path.SpeciesNum("e");
  double CUTOFF = 7.0; // setcutoff: spherical cutoff in angstroms

  for (int ptcl1Index=0; ptcl1Index<numChangedPtcls; ptcl1Index++){
    int ptcl1 = activeParticles(ptcl1Index);
/*   need to exclude intramolecular interactions - compare with line commented out below
      for(int i = 0;i<=4;++i){
        Path.DoPtcl(ptcl1+4*i) = false;
      }
*/
    Path.DoPtcl(ptcl1) = false;
    int species1=Path.ParticleSpeciesNum(ptcl1);
//    cerr << "Choosing particle " << ptcl1 << " which is of species " << species1 << endl;

    if (species1==speciesO){
      int counter = 0;
//    Calculate o-o interaction
      for (int ptcl2=Path.Species(speciesO).FirstPtcl;ptcl2<=Path.Species(speciesO).LastPtcl;ptcl2++) {///loop over oxygen
	if (Path.DoPtcl(ptcl2)){
//          int PairIndex = PairMatrix(species1, Path.ParticleSpeciesNum(ptcl2));
	  
	  for (int slice=startSlice;slice<=endSlice;slice+=skip){
	    dVec r, rp;
	    double rmag, rpmag;
	    PathData.Path.DistDisp(slice, ptcl1, ptcl2,
				   rmag, r);
//  implement spherical cutoff  
            if (rmag <= CUTOFF){
              counter ++;
              double sigma_over_r = PathData.Species(species1).Sigma/rmag;
	      double lj = 4*PathData.Species(species1).Epsilon*(pow(sigma_over_r,12)-pow(sigma_over_r,6)); // this is in kcal/mol 
	      TotalU += lj;
//	      cerr << "lj  " << lj << " at distance " << rmag << endl;
            }
	  }
	}
      }
      ///end calculating o-o interactions
//      cerr << "I calculated Lennard-Jones interactions with: " << counter << " particles" << endl;
    }
    else{
      /// calculating coulomb interactions
      for (int ptcl2=Path.Species(speciesp).FirstPtcl;ptcl2<=Path.Species(speciesp).LastPtcl;ptcl2++) {///loop over protons
//  don't compute intramolecular interactions
	if (Path.DoPtcl(ptcl2)&&Path.MolRef(ptcl1)!=Path.MolRef(ptcl2)){
	  for (int slice=startSlice;slice<=endSlice;slice+=skip){
	    double rmag;
	    dVec r;
	    PathData.Path.DistDisp(slice,ptcl1,ptcl2,rmag,r);
// implement spherical cutoff
            double Ormag = OOSeparation(slice,ptcl1,ptcl2);
            if (Ormag <= CUTOFF){
	      double coulomb_const = SI*angstrom_to_m*elementary_charge*elementary_charge*PathData.Species(species1).Charge*PathData.Species(speciesp).Charge*N_Avogadro/kcal_to_joule;
              double coulomb = coulomb_const/rmag;
	      TotalU += coulomb;
//	      cerr << "protons " << coulomb << " at distance " << rmag  << " with offset " << coulomb_const/CUTOFF << endl;
            }
          }
	}
      }
      for (int ptcl2=Path.Species(speciese).FirstPtcl;ptcl2<=Path.Species(speciese).LastPtcl;ptcl2++) {///loop over electrons
//  don't compute intramolecular interactions
	if (Path.DoPtcl(ptcl2)&&Path.MolRef(ptcl1)!=Path.MolRef(ptcl2)){
	  for (int slice=startSlice;slice<=endSlice;slice+=skip){
	    double rmag;
	    dVec r;
	    PathData.Path.DistDisp(slice,ptcl1,ptcl2,rmag,r);
// implement spherical cutoff
            double Ormag = OOSeparation(slice,ptcl1,ptcl2);
            if (Ormag <= CUTOFF){
	      double coulomb_const = SI*angstrom_to_m*elementary_charge*elementary_charge*PathData.Species(species1).Charge*PathData.Species(speciese).Charge*N_Avogadro/kcal_to_joule;
              double coulomb = coulomb_const/rmag;
	      TotalU += coulomb;
//	      cerr << "electrons " << coulomb << " at distance " << rmag << " with offset " << coulomb_const/CUTOFF << endl;
            }
	  }
	}
      }
    
	
    }
   /// end calculating coulomb interactions 
  
  
  }
//  double TotalU_times_beta = TotalU*kcal_to_joule/(N_Avogadro*k_B)*PathData.Path.tau;
//  cerr << TotalU << " and times beta " << TotalU_times_beta << " at temp " << 1.0/PathData.Path.tau << endl;
//  cerr << "Energy function is returning " << TotalU << endl;

  double energy_per_molecule = TotalU/PathData.Path.numMol;
//  return energy_per_molecule; // + thermal;
  return TotalU;
}

double TIP5PWaterClass::OOSeparation (int slice,int ptcl1,int ptcl2)
{
  int speciesO=Path.SpeciesNum("O");
  int Optcl1, Optcl2;
  dVec Or;
  double Ormag;

  Optcl1 = Path.Species(speciesO).FirstPtcl + Path.MolRef(ptcl1);
  Optcl2 = Path.Species(speciesO).FirstPtcl + Path.MolRef(ptcl2);
  PathData.Path.DistDisp(slice, Optcl1, Optcl2, Ormag, Or);
//  cerr << "We have particles " << ptcl1 << " and " << ptcl2 << " belonging to molecules " << Path.MolRef(ptcl1) << " and " << Path.MolRef(ptcl2) << " and calculate the distance between O " << Optcl1 << " and " << Optcl2 << " from " << Path.MolRef(Optcl1) << " and " << Path.MolRef(Optcl2) << endl;
  return Ormag;
}

void TIP5PWaterClass::Read (IOSectionClass &in)
{
  //do nothing for now
}

double TIP5PWaterClass::RotationalKinetic(int startSlice, int endSlice, const Array<int,1> &activeParticles,int level)
{
  double RotK = 0.0;
  int numChangedPtcls = activeParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl = activeParticles(ptclIndex);
    int species=Path.ParticleSpeciesNum(ptcl);
    double FourLambdaTauInv=1.0/(4.0*lambda_p*levelTau);
    for (int slice=startSlice; slice < endSlice;slice+=skip) {
      // Load the coordinates we'll be manipulating; these are coords for adjacent time slices of a proton.
      dVec coord1 = PathData.Path(slice,ptcl);
      dVec coord2 = PathData.Path(slice+skip,ptcl);
//cerr << "loaded " << coord1 << " and " << coord2 << " for slice " << slice << " and ptcl " << ptcl << endl;
      // Identify the index for the corresponding oxygen, i.e. COM
      int Optcl = FindCOM(ptcl);
      dVec Ocoord1 = PathData.Path(slice,Optcl);
      dVec Ocoord2 = PathData.Path(slice+skip,Optcl);
//cerr << "identified COM " << Optcl << " with coords " << Ocoord1 << " and " << Ocoord2 << endl;
      // We want to measure rotations relative to the COM ONLY, so we need to redefine coords WRT their respective COMs so that we can calculate polar and azimuthal angles WRT the COM.
      // Sync COMs for coord1 and coord2
//      dVec COMdisp = Displacement(slice,slice+skip,Optcl,Optcl);
//      coord2 -= COMdisp;
      // Redefine coordinates WRT COM
      coord1 -= Ocoord1;
      coord2 -= Ocoord2;
//cerr << "now WRT COM coords are " << coord1 << " and " << coord2 << endl;
      // Get the corresponding polar and azimuthal angles
      double theta1,theta2,phi1,phi2;
      GetAngles(coord1,theta1,phi1);
      GetAngles(coord2,theta2,phi2); 
//cerr << "corresponding angles are " << theta1 << ", " << phi1 << ", " << theta2 << ", " << phi2 << endl;
      // Calculate the rotational kinetic energy
      double vel_squared = CalcEnergy(theta1,theta2-theta1,phi2-phi1);
//cerr << "from which I calculate vel_squared " << vel_squared << endl;

      double GaussProd = 1.0;
//    for (int dim=0; dim<NDIM; dim++) {
//  	int NumImage=1;
      double GaussSum=0.0;
//	for (int image=-NumImage; image<=NumImage; image++) {
//	  double dist = vel[dim]+(double)image*Path.GetBox()[dim];
      GaussSum += exp(-vel_squared*FourLambdaTauInv);
//      }

      GaussProd *= GaussSum;
//      }
      RotK -= log(GaussProd);    
    //RotK += dot(vel,vel)*FourLambdaTauInv; 
    }
  }
  //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
//  //cerr << "I'm returning kinetic action " << RotK << endl;
  return (RotK);
}

double TIP5PWaterClass::RotationalEnergy(int startSlice, int endSlice, int level)
{
  double spring=0.0;
  // ldexp(double x, int n) = x*2^n
  double levelTau=ldexp(Path.tau, level);
//   for (int i=0; i<level; i++) 
//     levelTau *= 2.0;
  spring  = 0.0;  
  int skip = 1<<level;
  const int NumImage=1;  
  double Z = 0.0;
  double FourLambdaTauInv = 1.0/(4.0*lambda_p*levelTau);
  for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++) {
    // Do free-particle part
    int speciesNum  = Path.ParticleSpeciesNum(ptcl);
    if (speciesNum == PathData.Path.SpeciesNum("p")){
      SpeciesClass &species = Path.Species(speciesNum);
      for (int slice=startSlice; slice<endSlice; slice+=skip) {
	spring += (0.5*NDIM)/levelTau;

      // Load the coordinates we'll be manipulating; these are coords for adjacent time slices of a proton.
      dVec coord1 = PathData.Path(slice,ptcl);
      dVec coord2 = PathData.Path(slice+skip,ptcl);
//cerr << "loaded " << coord1 << " and " << coord2 << " for slice " << slice << " and ptcl " << ptcl << endl;
      // Identify the index for the corresponding oxygen, i.e. COM
      int Optcl = FindCOM(ptcl);
      dVec Ocoord1 = PathData.Path(slice,Optcl);
      dVec Ocoord2 = PathData.Path(slice+skip,Optcl);
//cerr << "identified COM " << Optcl << " with coords " << Ocoord1 << " and " << Ocoord2 << endl;
      // We want to measure rotations relative to the COM ONLY, so we need to redefine coords WRT their respective COMs so that we can calculate polar and azimuthal angles WRT the COM.
      // Sync COMs for coord1 and coord2
//      dVec COMdisp = Displacement(slice,slice+skip,Optcl,Optcl);
//      coord2 -= COMdisp;
      // Redefine coordinates WRT COM
      coord1 -= Ocoord1;
      coord2 -= Ocoord2;
//cerr << "now WRT COM coords are " << coord1 << " and " << coord2 << endl;
      // Get the corresponding polar and azimuthal angles
      double theta1,theta2,phi1,phi2;
      GetAngles(coord1,theta1,phi1);
      GetAngles(coord2,theta2,phi2); 
//cerr << "corresponding angles are " << theta1 << ", " << phi1 << ", " << theta2 << ", " << phi2 << endl;
      // Calculate the rotational kinetic energy
      double vel_squared = CalcEnergy(theta1,theta2-theta1,phi2-phi1);
//cerr << "from which I calculate vel_squared " << vel_squared << endl;

        double GaussSum;
        double numSum;
//      dVec GaussSum=0.0;
//      dVec numSum=0.0;
//      for (int dim=0; dim<NDIM; dim++) {
//        for (int image=-NumImage; image<=NumImage; image++) {
//	    double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	double d2overFLT = vel_squared*FourLambdaTauInv;
	double expPart = exp(-d2overFLT);
        GaussSum = expPart;
        numSum = -d2overFLT/levelTau* expPart; 
//	GaussSum[dim] += expPart;
//	numSum[dim] += -d2overFLT/levelTau* expPart;
//	}
//	Z *= GaussSum[dim];
        Z += GaussSum;
//      }
/*      double scalarnumSum=0.0;
        for (int dim=0;dim<NDIM;dim++) {
          dVec numProd=1.0;
	  for (int dim2=0;dim2<NDIM;dim2++) {
	    if (dim2!=dim)
	      numProd[dim] *= GaussSum[dim2];
	    else 
	      numProd[dim] *=  numSum[dim2];
	  }
	  scalarnumSum += numProd[dim];
        }
*/
	//cerr << "Z = " << Z << " scalarnumSum = " << scalarnumSum << endl;
        spring += numSum; 
      }
    }
  }
  spring = spring/Z;
//  cerr << "spring = " << spring << endl;
  return spring;
}

void TIP5PWaterClass::GetAngles(dVec disp, double &theta, double &phi)
{
  double x = disp(0);
  double y = disp(1);
  double z = disp(2);
  double R = O_H_moment_arm;
  theta = acos(z/R);
  double SineTheta = sin(theta);  
  phi = acos(x/(R*SineTheta));
  double phicheck = asin(y/(R*SineTheta));
  if(y<0){
    phi = 2*M_PI - phi;
  }
}

double TIP5PWaterClass::CalcEnergy(double reftheta,double dtheta, double dphi)
{
  double R = O_H_moment_arm;
  double SineTheta = sin(reftheta);
  double omega_squared = dtheta*dtheta + dphi*dphi*SineTheta*SineTheta;
  double vel_squared = omega_squared*R*R;
  return vel_squared;
}

double TIP5PWaterClass::SineOfPolar(dVec coords)
{
  double x = coords(0);
  double z = coords(2);
  double R = O_H_moment_arm;
  double theta = acos(z/R);  
cerr << "I find theta = " << theta << " for coords " << x << ", " << coords(1) << ", " << z << ";";
  double SineTheta = sin(theta);
cerr << "returning " << SineTheta << endl;
  return SineTheta;
}

dVec TIP5PWaterClass::COMVelocity (int slice1,int slice2,int ptcl)
{
  int speciesO=PathData.Path.SpeciesNum("O");
  int Optcl;
  dVec Ovel;

  Optcl = Path.Species(speciesO).FirstPtcl + Path.MolRef(ptcl);
  Ovel = PathData.Path.Velocity(slice1, slice2, Optcl);
//cerr << "I'm correcting velocity for ptcl " << ptcl << " of species " << Path.ParticleSpeciesNum(ptcl) << ".  Found COM oxygen at " << Optcl;
//cerr << "Returning COM velocity " << Ovel << endl;
  return Ovel;
}

dVec TIP5PWaterClass::COMCoords (int slice, int ptcl)
{
  int speciesO=PathData.Path.SpeciesNum("O");
  int Optcl;
  dVec relative_coords;

  Optcl = Path.Species(speciesO).FirstPtcl + Path.MolRef(ptcl);
  relative_coords = PathData.Path(slice,ptcl) - PathData.Path(slice,Optcl);
//cerr << "I'm correcting velocity for ptcl " << ptcl << " of species " << Path.ParticleSpeciesNum(ptcl) << ".  Found COM oxygen at " << Optcl;
//cerr << "Returning COM velocity " << Ovel << endl;
  return relative_coords;
}

dVec TIP5PWaterClass::Displacement(int slice1, int slice2, int ptcl1, int ptcl2)
{
  dVec disp;
  disp = PathData.Path(slice1,ptcl1) - PathData.Path(slice2,ptcl2);
  return disp;
}

int TIP5PWaterClass::FindCOM(int ptcl)
{
  int speciesO=PathData.Path.SpeciesNum("O");
  int Optcl;
  Optcl = Path.Species(speciesO).FirstPtcl + Path.MolRef(ptcl);
  return Optcl;
}

double TIP5PWaterClass::ProtonKineticAction (int slice1, int slice2, const Array<int,1> &changedParticles, int level)
{
  double TotalK = 0.0;
  int numChangedPtcls = changedParticles.size();
  int skip = 1<<level;
  double levelTau = Path.tau* (1<<level);
  for (int ptclIndex=0; ptclIndex<numChangedPtcls; ptclIndex++){
    int ptcl = changedParticles(ptclIndex);
    int species=Path.ParticleSpeciesNum(ptcl);
    double lambda = lambda_p;
    if (lambda != 0){
      double FourLambdaTauInv=1.0/(4.0*lambda*levelTau);
      for (int slice=slice1; slice < slice2;slice+=skip) {
      // This is the same as Kinetic.Action except that coordinates are WRT their COMs.
        dVec coord1 = PathData.Path(slice,ptcl);
        dVec coord2 = PathData.Path(slice+skip,ptcl);
        int Optcl = FindCOM(ptcl);
        dVec Ocoord1 = PathData.Path(slice,Optcl);
        dVec Ocoord2 = PathData.Path(slice+skip,Optcl);
        coord1 -= Ocoord1;
        coord2 -= Ocoord2;
        dVec vel;
        vel = coord1 - coord2;
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
      }
    }
  }
  return (TotalK);
}

double TIP5PWaterClass::ProtonKineticEnergy (int slice1, int slice2, int level)
{
  double spring=0.0;
  double levelTau=ldexp(Path.tau, level);
  spring  = 0.0;  
  int skip = 1<<level;
  const int NumImage=1;  
  for (int ptcl=0; ptcl<Path.NumParticles(); ptcl++) {
    int speciesNum  = Path.ParticleSpeciesNum(ptcl);
    SpeciesClass &species = Path.Species(speciesNum);
    double lambda = lambda_p;
    if (lambda != 0.0) {
      double FourLambdaTauInv = 1.0/(4.0*lambda*levelTau);
      for (int slice=slice1; slice<slice2; slice+=skip) {
	spring += (0.5*NDIM)/levelTau;
      // This is the same as Kinetic.Action except that coordinates are WRT their COMs.
        dVec coord1 = PathData.Path(slice,ptcl);
        dVec coord2 = PathData.Path(slice+skip,ptcl);
        int Optcl = FindCOM(ptcl);
        dVec Ocoord1 = PathData.Path(slice,Optcl);
        dVec Ocoord2 = PathData.Path(slice+skip,Optcl);
        coord1 -= Ocoord1;
        coord2 -= Ocoord2;
        dVec vel;
        vel = coord1 - coord2;
	double Z = 1.0;
	dVec GaussSum=0.0;
	dVec numSum=0.0;
	for (int dim=0; dim<NDIM; dim++) {
	  for (int image=-NumImage; image<=NumImage; image++) {
	    double dist = vel[dim]+(double)image*PathData.Path.GetBox()[dim];
	    double d2overFLT = dist*dist*FourLambdaTauInv;
	    double expPart = exp(-d2overFLT);
	    GaussSum[dim] += expPart;
	    numSum[dim] += -d2overFLT/levelTau* expPart;
	  }
	  Z *= GaussSum[dim];
	}
	double scalarnumSum=0.0;
	for (int dim=0;dim<NDIM;dim++) {
	  dVec numProd=1.0;
	  for (int dim2=0;dim2<NDIM;dim2++) {
	    if (dim2!=dim)
	      numProd[dim] *= GaussSum[dim2];
	    else 
	      numProd[dim] *=  numSum[dim2];
	  }
	  scalarnumSum += numProd[dim];
	} //cerr << "Z = " << Z << " scalarnumSum = " << scalarnumSum << endl;
	spring += scalarnumSum/Z; 
      }
    }
  }
  //  cerr << "spring = " << spring << endl;
  return spring;
}
