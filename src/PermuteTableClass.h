#ifndef PERMUTE_TABLE_CLASS_H
#define PERMUTE_TABLE_CLASS_H

// class Hclass
// {
//  public:
//   int j;
//   double Htilde_ij;
// };

class PermClass
{
 public:
  int Ncycles;
  TinyVector<int,4> Perm;
  double P, C;
}


class PermuteTableClass
{
 private:
  int TableSize;
  int NumEntries;
  int 
  inline void AddEntry(const PermClass &perm);
 public:
  int Species;
  PathDataClass &PathData;
  int Slice1, NumLevels;
  /// The enhancement factor for each cycle length. That is we
  /// multiply the probability of a 3-cycle by Gamma[2]  
  TinyVector<double,4> gamma; 

  // Smallest value of exp(-s^2/(4*lambda*beta) which we will include
  // in Htable.
  double epsilon;
  Array<double,2> HTable;
  Array<PermClass,1> PermTable;

  
  void ConstructHTable(int slice);
  void CalcPermProbs();
  PermuteTableClass()
  {
    TableSize = 1000;
    PermTable.resize(TableSize);
  }
  inline int FindEntry(double xi);
};

inline void PermuteTableClass::AddEntry(PermClass &perm)
{
  if (NumEntries >= (TableSize-1)) {
    TableSize *= 2;
    PermTable.resizeAndPreserve(TableSize);
  }

  PermTable(NumEntries) = perm;
  NumEntries++;
}

inline int PermTableClass::FindEntry(double xi)
{
  // Do a binary search
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
