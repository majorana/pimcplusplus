#include "StructureFactor.h"

void StructureFactorClass::Read(IOSectionClass& in)
{
  
  ObservableClass::Read(in);
  string species1Name;
  string species2Name;
  Species1=-1;
  Species2=-1;
  if (in.ReadVar("Species", species1Name)) {
    species2Name=species1Name;
  }
  else {
    assert(in.ReadVar("Species1",species1Name));
    assert(in.ReadVar("Species2",species2Name));
  }
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
  for (int spec=0;spec<PathData.NumSpecies();spec++){
    if (PathData.Species(spec).Name==species1Name){
      Species1=spec;
    }
    if (PathData.Species(spec).Name==species2Name){
      Species2=spec;
    }
  }
  assert(Species1!=-1);
  assert(Species2!=-1);

  TotalCounts=0;
  Sk.resize(PathData.Path.kVecs.size());
  Sk=0;
}



void StructureFactorClass::WriteInfo()
{
  PathClass &Path = PathData.Path;
  Array<dVec,1> &kVecs = PathData.Path.kVecs;
  ObservableClass::WriteInfo();
  Array<double,2> kVecArray(kVecs.size(),NDIM);
  Array<double,1> kMagArray(kVecs.size());

  for (int ki=0; ki < kVecs.size(); ki++) {
    dVec &k = kVecs(ki);
    kMagArray(ki) = sqrt(dot(k,k));
    for (int j=0; j<NDIM; j++)
      kVecArray(ki,j) = k[j];
  }

  IOSection.WriteVar("kVecs", kVecArray);
  
  /// Output data for plotting in analysis code
  IOSection.WriteVar("x", kMagArray);
  IOSection.WriteVar("xlabel", "|k|");
  IOSection.WriteVar("ylabel", "S(k)");
  IOSection.WriteVar("Species1", PathData.Species(Species1).Name);
  IOSection.WriteVar("Species2", PathData.Species(Species2).Name);
  IOSection.WriteVar("Type","CorrelationFunction");
}


void StructureFactorClass::WriteBlock()
{
  Array<dVec,1> &kVecs = PathData.Path.kVecs;
  Array<double,1> SkSum(kVecs.size());
  double norm=0.0;
  int num1 = PathData.Path.Species(Species1).NumParticles;
  int num2 = PathData.Path.Species(Species1).NumParticles;
  norm = TotalCounts * sqrt((double)num1*num2);
  
  PathData.Communicator.Sum(Sk, SkSum);
  if (PathData.Communicator.MyProc()==0) {
    if (FirstTime) {
      FirstTime=false;
      WriteInfo();
      Array<double,2> SofkArray(1,kVecs.size());
      for (int ki=0; ki<kVecs.size(); ki++)
	SofkArray(0,ki) = (double) SkSum(ki) / norm;
      IOSection.WriteVar("y",SofkArray);
      IOVar = IOSection.GetVarPtr("y");
    }
    else {
      Array<double,1> SofkArray(SkSum.size());
      for (int ki=0; ki<kVecs.size(); ki++)
	SofkArray(ki) = (double) SkSum(ki) / norm;
      IOVar->Append(SofkArray);
    }
  }
}




void StructureFactorClass::Accumulate()
{
  Array<dVec,1> &kVecs = PathData.Path.kVecs;

  SpeciesClass &species1=PathData.Path.Species(Species1);
  SpeciesClass &species2=PathData.Path.Species(Species2);

  TimesCalled++;

  if (TimesCalled % DumpFreq==0){
    WriteBlock();
  }
  if ((TimesCalled % Freq)!=0){
    return;
  }

  if (!PathData.Path.LongRange) {
    for (int slice=0; slice < PathData.NumTimeSlices()-1; slice++)
      PathData.Path.CalcRho_ks_Fast(slice, Species1);
    if (Species2 != Species1)
      for (int slice=0; slice < PathData.NumTimeSlices()-1; slice++)
	PathData.Path.CalcRho_ks_Fast(slice, Species2);
  }

  for (int slice=0;slice<PathData.NumTimeSlices()-1;slice++) {
    TotalCounts++;
    for (int ki=0; ki<kVecs.size(); ki++) {
      double a = PathData.Path.Rho_k(slice, Species1, ki).real();
      double b = PathData.Path.Rho_k(slice, Species1, ki).imag();
      double c = PathData.Path.Rho_k(slice, Species2, ki).real();
      double d = PathData.Path.Rho_k(slice, Species2, ki).imag();
      // \f$ Sk(ki) :=  Sk(ki) + \Re(rho^1_k * rho^2_{-k}) \f
      Sk(ki) += a*c + b*d;	
    }
  }
}


void StructureFactorClass::Initialize()
{
  TotalCounts = 0;
  TimesCalled=0;

}


