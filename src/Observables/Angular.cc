#include "Angular.h"

///////////////////////////////////////////////////////
//Angular Correlation Integral[f(t)*f(0),{t,0,beta}]///
//where f(t) = Sum[dot[R(t),V(t)] over ptcl & slices]//
///////////////////////////////////////////////////////




void AngularClass::WriteBlock()
{
  //  Array<double,1> CorSum(Correlation.size());
  //  Path.Communicator.Sum(Correlation, CorSum);

  Correlation=Correlation/TotalCounts;
  CorVar.Write(Correlation);
  CorVar.Flush();
  Correlation=0;
  TotalCounts=0;

}


void AngularClass::Read(IOSectionClass &in)
{  

  ObservableClass::Read(in);
  string speciesName;
  Species=-1;
  assert(in.ReadVar("Species1",speciesName));
  for (int spec=0;spec<PathData.NumSpecies();spec++){ //???what is Species ?
    if (PathData.Species(spec).Name==speciesName){
      Species=spec;
    }
  }
  if (PathData.Path.Communicator.MyProc()==0){
    IOSection.WriteVar("Type","CorrelationFunction");
    IOSection.WriteVar("Cumulative", false);
  }
  Correlation.resize(3,PathData.Path.NumTimeSlices()+1);
  Correlation=0;
  TotalCounts=0;
  /// Now write the one-time output variables
//   if (PathData.Path.Communicator.MyProc()==0)
//     WriteInfo();

}

void AngularClass::Accumulate()
{
  TotalCounts++;
  PathClass &Path= PathData.Path;
  SpeciesClass &species=PathData.Path.Species(Species);
  dVec r;
  dVec angMom1,angMom2;
  dVec v;
  int slicept;
  for (int l=0;l<=2;l++)
  for (int t=0;t<=PathData.NumTimeSlices();t++){
    for (int ptcl=species.FirstPtcl;ptcl<=species.LastPtcl;ptcl++){
      for (int ptcl2=species.FirstPtcl;ptcl2<=species.LastPtcl;ptcl2++)
	for (int slice=0;slice<=PathData.NumTimeSlices()-1;slice++){
	  slicept=(slice+t)%PathData.NumTimeSlices();	
	  dVec r1=Path(slice,ptcl);
	  dVec r2=Path(slicept,ptcl2);
	  double cosGamma=dot(r1,r2)/(dot(r1,r1)*dot(r2,r2));
	  ///Ignoring a 4\Pi normalization here
	  Correlation(l,t)+=(2*l+1)*Legendre(l,cosGamma);	  
      }
    }
  }    
}


