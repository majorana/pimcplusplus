#include "HBond.h"
#include <iostream>
#include <fstream>
#include <string>

Array<int,2> Table;
Array<int,2> BondCount;
Array<int,1> Histogram;
Array<int,1> LifetimeHist;

  //////////////////////////////////////
 ///  Hydrogen Bond Analysis Class  ///
//////////////////////////////////////

bool HbondClass::IsHBond(int slice, int OA, int OB){
  bool bond = false;
  int numMol = PathData.Path.numMol;
  int PA1 = OA + 3*numMol; 
  int PA2 = OA + 4*numMol; 
  int PB1 = OB + 3*numMol; 
  int PB2 = OB + 4*numMol; 
  if (CheckPair(slice,OA,OB,PB1))
    bond = true;
  else if (CheckPair(slice,OA,OB,PB2))
    bond = true;
  else if (CheckPair(slice,OB,OA,PA1))
    bond = true;
  else if (CheckPair(slice,OB,OA,PA2))
    bond = true;
  return bond;
}

bool HbondClass::CheckPair(int slice, int obond, int ohome, int p){
// Here we use HBond criteria r_OO < 3.5 and HOH bond angle > 145 deg, after Artacho (2004)
  double OOlimit = 3.5;
  double HOHangle = 2*M_PI*145/360;
  bool bond = false;
  double OHmag,OHbondmag,OOmag;
  dVec OH,OHbond,OO;
  PathData.Path.DistDisp(slice,p,obond,OHbondmag,OHbond);
  PathData.Path.DistDisp(slice,p,ohome,OHmag,OH);
  PathData.Path.DistDisp(slice,obond,ohome,OOmag,OO);
  double theta = PathData.Actions.TIP5PWater.GetAngle(OHbond,OH);
  if ((OOmag < OOlimit) && (theta > HOHangle))
    bond = true;
  return bond; 
}

void HbondClass::Read(IOSectionClass& in)
{
cerr << "HBond: Reading..." << endl;  
  ObservableClass::Read(in);
  assert(in.ReadVar("freq",Freq));
  assert(in.ReadVar("dumpFreq",DumpFreq));
cerr << "Freq is " << Freq << " and DumpFreq is " << DumpFreq << endl;
  Table.resize(PathData.Path.numMol,PathData.Path.numMol);
  Table=0;
cerr << "resized Table: " << Table.size() << endl;
  BondCount.resize(PathData.Path.numMol,PathData.Path.numMol);
  BondCount=0;
cerr << "resized BondCount: " << BondCount.size() << endl;
  LifetimeHist.resize(DumpFreq/Freq);
  LifetimeHist = 0;
cerr << "resized LifetimeHist: " << LifetimeHist.size() << endl;
  in.CloseSection();
cerr << "leaving Read" << endl;
}

void HbondClass::WriteBlock()
{
  // IO Stuff; Open files
  std::fstream bondout,lifetimeout;
  bondout.open("./HBondCount.dat", std::ios::out);
  lifetimeout.open("./HBondLifetime.dat", std::ios::out);

  int countNorm = DumpFreq/Freq;
  int lifetimeNorm = 1;

  // Write files
  for (int i = 0; i < Histogram.size(); i++){
    bondout << i << " " << Histogram(i) << endl;
  }
  for (int j = 0; j < LifetimeHist.size(); j++){
    lifetimeout << j << " " << LifetimeHist(j) << endl;
  }

  bondout.close();
  lifetimeout.close();
}

void HbondClass::Accumulate()
{
  // Table accumulates the lifetime of an HBond and is cleared out with period DumpFreq.  BondCount counts the presence or absence of a bond and is cleared with period Freq.

  // Measure and tabulate hbonds at the curent time step

    // Just do this for slice 0 -- classical 
    int slice = 0;
    /// loop over molecules 
    for (int mol1=0;mol1<PathData.Path.numMol;mol1++){
      for (int mol2=mol1+1;mol2<PathData.Path.numMol;mol2++){
        if(IsHBond(slice,mol1,mol2)){
          Table(mol1,mol2)++;
          BondCount(mol1,mol2)++;
        } 
      }
    }

    // Tabulate the number of hbonds and accumulate in Histogram
    int count = 0;
    for (int mol = 0; mol < BondCount.extent(0); mol++){
      // Sum entries over rows
      for (int m = 0; m < BondCount.extent(0); m++)
        count += BondCount(m,mol);
      // Sum entries over columns
      for (int n = 0; n < BondCount.extent(1); n++)
        count += BondCount(mol,n);
      Histogram(count)++;
      count = 0;
    }
    BondCount = 0;
  }  

  // Tabulate hbond lifetime histogram and write data to file
  if ((TimesCalled % DumpFreq) == 0){
    cerr << TimesCalled << ": " << TimesCalled << "; Writing HBond Data" << endl;
    // Tabulate hbond lifetimes and accumulate in LifetimeHist
    int count = 0;
    for (int mol = 0; mol < Table.extent(0); mol++){
      // Sum entries over rows
      for (int m = 0; m < Table.extent(0); m++)
        count += Table(m,mol);
      // Sum entries over columns
      for (int n = 0; n < Table.extent(1); n++)
        count += Table(mol,n);
      LifetimeHist(count)++;
      count = 0;
    }

    // Write files
    WriteBlock();
    Histogram = 0;
    LifetimeHist = 0;

}

void HbondClass::Initialize()
{
cerr << "HBond: Initialize" << endl;
  TimesCalled=0;
  Histogram.resize(10);
  Histogram = 0;
cerr << "Leaving Initialize" << endl;
}
