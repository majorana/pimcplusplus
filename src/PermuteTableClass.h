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
  int Ncycles;
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
  double Norm, NormInv;
  inline void AddEntry(const CycleClass &cycle);
  inline int FindEntry(double xi);  
  PathDataClass &PathData;
  void ConstructHTable();


  /// Stores which Htable and PermTable is for the forward move
  /// 0 or 1.  Flips when a permutation is accepted.
 public:
  int NumEntries;

  Array<double,2> HTable;
  Array<CycleClass,1> PermTable;

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
  double CalcReverseProb(const CycleClass &myPerm,
			 const PermuteTableClass &forwardTable);
  PermuteTableClass(PathDataClass &myPathData) : PathData(myPathData)
  {
    NumEntries=0;
    TableSize = 1000;
    PermTable.resize(TableSize);
  }

};

inline void PermuteTableClass::AddEntry(const CycleClass &cycle)
{
  if (NumEntries >= (TableSize-1)) {
    TableSize *= 2;
    PermTable.resizeAndPreserve(TableSize);
  }

  PermTable(NumEntries) = cycle;
  NumEntries++;
}

//Pass a random number between 0 and 1 to this function and it
//returns the index of the permutation with the appropriate
//probability. 
inline int PermuteTableClass::FindEntry(double xi) 
{
  // Do a binary search
  xi *= Norm; 
  int hi = NumEntries;
  int lo = 0;
  int attempt;
  while (hi != lo) {
    attempt = (hi+lo)>>1;
    if (PermTable(attempt).C > xi)
      hi = attempt;
    else
      lo = attempt;
  }
  return (lo);
}


#endif
