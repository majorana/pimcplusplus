#include "Common/Blitz.h"
#include "PathDataClass.h"
#ifndef PERMUTE_TABLE_CLASS_H
#define PERMUTE_TABLE_CLASS_H

// class Hclass
// {
//  public:
//   int j;
//   double Htilde_ij;
// };

/// Note:  all indices in this file are with respect to the present
/// species only.  Thus, when accessing the actual path, we must add
/// the PathData.Species(SpeciesNum).FirstPtcl to our indices.

class CycleClass
{
 public:
  int Length;
  TinyVector<int,4> CycleRep;
  double P, C;
  ///Takes the number of particles you have and returns the
  ///representation such that p(j) is the particle that the j'th
  ///particle permutes onto. 
  //  void CanonicalPermRep(Array<int,1> myArray);
  void Apply (PathClass &path, int firstPtcl, int timeSlice);

};




class PermuteTableClass
{
 private:
  int TableSize;
  int SpeciesNum;
  int Slice1, Slice2;
  inline void AddEntry(const CycleClass &cycle);

  PathDataClass &PathData;
  void ConstructHTable();

 public:
  double Norm, NormInv;
  inline int FindEntry(double xi);  


  int NumEntries;
  CycleClass CurrentCycle;
  Array<double,2> HTable;
  Array<CycleClass,1> CycleTable;

  //  void PermuteHTable();

  /// The enhancement factor for each cycle length. That is we
  /// multiply the probability of a 3-cycle by Gamma[2]  
  TinyVector<double,4> Gamma; 
  
  // Smallest value of exp(-s^2/(4*lambda*beta) which we will include
  // in Htable.
  double epsilon;
  void Read(IOSectionClass &inSection);
  void ConstructCycleTable(int speciesNum,int slice1,int slice2);
  void CanonicalPermRep(Array<int,1> P);
  double AttemptPermutation();
  double CalcReverseProb(const PermuteTableClass &forwardTable);
  Array<int,1> CurrentParticles();
  PermuteTableClass(PathDataClass &myPathData) : PathData(myPathData)
  {
    NumEntries=0;
    TableSize = 1000;
    CycleTable.resize(TableSize);
  }

};

inline void PermuteTableClass::AddEntry(const CycleClass &cycle)
{
  if (NumEntries >= (TableSize-1)) {
    TableSize *= 2;
    CycleTable.resizeAndPreserve(TableSize);
  }

  CycleTable(NumEntries) = cycle;
  NumEntries++;
}

//Pass a random number between 0 and 1 to this function and it
//returns the index of the permutation with the appropriate
//probability. 
inline int PermuteTableClass::FindEntry(double xi) 
{
  // Do a binary search
  xi *= Norm; 
  int hi = NumEntries-1;
  int lo = 0;
  int attempt = (hi+lo)>>1;
  while (attempt != lo) {
    attempt = (hi+lo)>>1;
    if (CycleTable(attempt).C > xi)
      hi = attempt;
    else
      lo = attempt;
  }
  return (hi);
}


#endif
