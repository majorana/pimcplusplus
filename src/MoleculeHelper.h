#ifndef MOLECULE_HELPER_H
#define MOLECULE_HELPER_H

#include "Common.h"
#include <Common/IO/IO.h>

using namespace IO;

class MoleculeManagerClass
{
  int numMolTypes;
  int totalNumMol;
  Array<string,1> names;
  Array<int,1> num;
  Array<Array<int,1>, 1> Members;
  Array<Array<int,1>, 1> ListByType;
  Array<int,1> MolRef;
  Array<string,1> MolLabel;

  public:

  void Read(IOSectionClass& in);
  void Init();
  Array<int,1>& MembersOf(int mol);
  int SizeOf(int mol);
  int NumMol(int type);
  int NumMol(string typeLabel);
  int NumMol();
  string NameOf(int mol);
  Array<int,1>& MolOfType(int type);
  Array<int,1>& MolOfType(string typeLabel);
  int Index(string label);
  inline const int operator() (int ptcl);

  MoleculeManagerClass();
};

inline const int MoleculeManagerClass::operator()(int ptcl)
{
  return MolRef(ptcl);
}

#endif
