#ifndef INPUT_OUTPUT_HDF5_H
#include "InputOutput.h"
#include <iostream>
#include <stack>
#include "hdf5.h"


class VarHDF5Class : public VarClass
{
public:
  bool ReadInto (double &val) { }
  bool ReadInto (int &val) { }
  bool ReadInto (string &val) { }
};


class SectionNumberPairClass
{
public:
  string Name;
  int CurrentNum;
  SectionNumberPairClass()
  { CurrentNum = 0;}
};

class SectionNumberClass
{
private:
  list<SectionNumberPairClass> SecNumList;
public:
  int operator()(string name)
  {
    int num = 0;
    list<SectionNumberPairClass>::iterator iter;
    iter = SecNumList.begin();
    while ((iter != SecNumList.end()) && (iter->Name != name))
      iter++;
    if (iter != SecNumList.end())  // We already have a section of 
      {                            // this name
	iter->CurrentNum++;
	num = iter->CurrentNum;
      }
    else
      {
	SectionNumberPairClass newSec;
	newSec.Name = name;
	newSec.CurrentNum = 0;
	SecNumList.push_back(newSec);
      }
    return (num);
  }
};


class SectionClass
{
public:
  hid_t GroupID;
  SectionNumberClass SectionNumber;
};


class OutputSectionHDF5Class : public OutputSectionClass
{
private:
  hid_t FileID;
  string FileName;
  bool IsOpen;
  stack<SectionClass> SectionStack;
public:
  bool OpenFile (string fileName);
  void OpenSection (string name);
  void CloseSection ();
  void WriteVar (string name, double T);
  void CloseFile();
  OutputSectionHDF5Class() : IsOpen(false)
  { }
};


#endif
