#include "ActionClass.h"
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>
#include "Common/PairAction/DavidPAClass.h"

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
  string fileName;
  //  if(inSection.ReadVar("Daviddmfile",fileName)){
  //    DavidPAClass myPAClass;
  //    myPAClass.Print();
  //    myPAClass.ReadDavidSquarerFile(fileName.c_str());
    //    assert(inSection.ReadVar("type1",type1));
    //    assert(inSection.ReadVar("type2",type2));
  //  }

  
  



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

double ActionClass::calcTotalAction(int startSlice, int endSlice, 
				    Array<int,1> changedParticles,
				    int level)//, double &PE, double &KE)
{
  double PE, KE;
  //  cerr<<"My changed particle is "<<changedParticles(0)<<endl;
  // First, sum the pair actions
  for (int counter=0;counter<Path.DoPtcl.size();counter++){
    Path.DoPtcl(counter)=true;
  }
  PE = 0.0;
  KE = 0.0;
  double TotalU = 0.0;
  double TotalK = 0.0;
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
	//cerr<<"PairIndex: "<<PairIndex<<endl;
	for (int slice=startSlice;slice<endSlice;slice+=skip){
	  dVec r1=Path(slice,ptcl1);
	  dVec r2=Path(slice,ptcl2);
	  dVec rp1=Path(slice+skip,ptcl1);
	  dVec rp2=Path(slice+skip,ptcl2);

	  dVec r, rp;
	  double rmag, rpmag;

	  DistanceTable->DistDisp(slice, slice+skip, ptcl1, ptcl2,
				  rmag, rpmag, r, rp);
	  //	  //	  r=r2-r1;
	  //	  //	  rp=rp2-rp1;
	  //	  //	  rmag=sqrt(dot(r,r));
	  //	  //	  rpmag=sqrt(dot(rp,rp));
	  //	  cerr<<"rmag "<<rmag<<endl;
	  //	  cerr<<"rpmag "<<rpmag<<endl;
 
	  double s2 = dot (r-rp, r-rp);
	  double q = 0.5 * (rmag + rpmag);
	  double z = (rmag - rpmag);
	  double U, dU, V;
	  U = PairActionVector(PairIndex)->U(q,z,s2, level);//, U, dU, V);
	  if (U!=0){
	    //	    cerr<<"q, z, s2, U: "<<q<<" "<<" "<<z<<" "<<s2<<" "<<U<<endl;
	    //	    cerr<<r<<rp<<endl;
	  }
	  //	  if (((ptcl1==1) && (ptcl2==2)) || ((ptcl1==0) && (ptcl2==3)))
	  TotalU += U;
	  PE += V;
	  KE -= dU;

	}
      }
      
    }
    double FourLambdaTauInv=1.0/(4.0*Path.Species(species1).lambda*levelTau);
    for (int slice=startSlice; slice < endSlice;slice+=skip) {
      dVec vel;
      vel = DistanceTable->Velocity(slice, slice+skip, ptcl1);
      double GaussProd = 1.0;
      for (int dim=0; dim<NDIM; dim++) {
	int NumImage=1;
	double GaussSum=0.0;
	for (int image=-NumImage; image<=NumImage; image++) {
	  double dist = vel[dim]+(double)image*Path.Box[dim];
	  GaussSum += exp(-dist*dist*FourLambdaTauInv);
	}
	GaussProd *= GaussSum;
      }
      TotalK -= log(GaussProd);
      //vel = Path(slice+skip,ptcl1) - Path(slice,ptcl1);
      //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
      //TotalK += dot(vel,vel)*FourLambdaTauInv; 
    }
  }
  KE += TotalK / levelTau;
  //  static int count=0;
  //  count++;
  //  if (count % 5000==0){
  //    cerr<<"TotalK:  "<<TotalK<<endl;
  //    cerr<<"TotalU:  "<<TotalU<<endl;
  //  }
  return (TotalK + TotalU);
  
  
}


  
