#include "DavidLongRangeClass.h"
#include "../PathDataClass.h"


void DavidLongRangeClass::Read(IOSectionClass &in)
{
  double myNum;
  cerr<<"The current size of the thing is "<<Path.kVecs.size()<<endl;
  uk.resize(Path.kVecs.size());
  duk.resize(Path.kVecs.size());
  for (int counter=0;counter<duk.size();counter++){
    duk(counter)=0.0;
  }

  ifstream infile;
  ///BUG: Currently hardcoded for actual file
  infile.open("bc70r10.lr");
  cerr<<"Beginning now"<<endl;
  for (int lvl=0;lvl<2;lvl++)
    for (int isEnergy=0;isEnergy<3;isEnergy++)
      ///BUG: 
      ///Currently hard coded for 20. Ugly 
      for (int kVec=0;kVec<20;kVec++){
	infile>>myNum;
	if (lvl==1 && isEnergy==1){
	  uk(kVec)=myNum;
	  cout<<myNum<<endl;
	}
	if (lvl==1 && isEnergy!=1){
	  duk(kVec)+=myNum;
	}

      }
  cerr<<"Done"<<endl;
  infile.close();
  
}

inline double mag2 (const complex<double> &z)
{
  return (z.real()*z.real() + z.imag()*z.imag());
}


///Calculates the long range part of the action using David's breakup.
///The short range part must be supplied as a dm file without the long
///range part in it.  It ignores active particles.
double DavidLongRangeClass::Action (int slice1, int slice2, 
	       const Array<int,1> &activeParticles, int level)
{
  //  cerr<<"My level is "<<level<<endl;
  double total=0;
  double factor;
  for (int slice=slice1;slice<=slice2;slice++){
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;

    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      //      PairActionFitClass &pa = *PairMatrix(species,species);
      //      if (pa.IsLongRange()) {
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
	double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	//	cerr<<"My ki is "<<ki<<endl;
	//	cerr<<"The spot I'm acessing is "<<Path.MagKint(ki)<<endl;
	//	cerr<<"The value of this spot is "<<uk(Path.MagKint(ki))<<endl;
	total +=  factor*rhok2 * uk(Path.MagKint(ki));
	
      }
    }
  }
  //  cerr<<"I am being called"<<endl;
  //  cerr<<"My total is "<<total;
  return total;

}

  ///Not really d_dbeta but total energy
double DavidLongRangeClass::d_dBeta (int slice1, int slice2,  int level)
{


  double total=0;
  double factor;
  for (int slice=slice1;slice<=slice2;slice++){
    if ((slice == slice1) || (slice==slice2))
      factor = 0.5;
    else
      factor = 1.0;
    for (int species=0; species<Path.NumSpecies(); species++) {
      Path.CalcRho_ks_Fast(slice,species);
      //      PairActionFitClass &pa = *PairMatrix(species,species);
      //      if (pa.IsLongRange()) {
      for (int ki=0; ki<Path.kVecs.size(); ki++) {
	double rhok2 = mag2(Path.Rho_k(slice,species,ki));
	//	cerr<<"My ki is "<<ki<<endl;
	//	cerr<<"The spot I'm acessing is "<<Path.MagKint(ki)<<endl;
	//	cerr<<"The value of this spot is "<<uk(Path.MagKint(ki))<<endl;
	total +=  factor*rhok2 * duk(Path.MagKint(ki));
	
      }
    }
  }
  //  cerr<<"I am being called"<<endl;
  return total;

}

DavidLongRangeClass::DavidLongRangeClass (PathDataClass &pathData) :
  ActionBaseClass (pathData)
{
  //Do  nothing for now
}
