#ifndef INPUT_OUTPUT_HDF5_H
#include "InputOutput.h"
#include <iostream>
#include <stack>
#include "hdf5.h"


/************************************************************
 *                    Input Functions                      *
 ************************************************************/

/// This class holds a reference to a dataset within an
/// HDF5 file that may be loaded at a later time.  The ReadInto
/// functions actually read the data from the HDF5 file into
/// the variable passed by reference.
class VarHDF5Class : public VarClass
{
public:
  hid_t DataSetID;
  H5T_class_t TypeClass;
  int Ndims;
  Array<hsize_t,1> Dimensions;


  bool ReadInto (double &val) { }
  bool ReadInto (int &val) { }
  bool ReadInto (string &val) { }
};


class InputSectionHDF5Class : public InputSectionClass
{
private:
  bool IsOpen, IsRoot;
  hid_t GroupID;

  void ReadGroup (hid_t GroupID, string name, InputSectionClass *parent);
public:
  void PrintTree(InputSectionClass *sec)
  { cerr << "Tree";}
  void GroupIterator (string member_name);
  bool OpenFile (string fileName, InputSectionClass *parent);
  void CloseFile ();
  InputSectionHDF5Class() : IsOpen(false), IsRoot(false)
  { }
};





/************************************************************
 *                    Output Functions                      *
 ************************************************************/

/// In HDF format, we cannot have the two groups with the same name in
/// the same level of hierarchy.  Therefore we implement a kludge in
/// which we append ".#" to each section name, where # is a number
/// which starts from 0 and continues upward as we add more sections
/// with the same name.  This class stores a name and the current
/// number we're at.
class SectionNumberPairClass
{
public:
  string Name;
  int CurrentNum;
  SectionNumberPairClass()
  { CurrentNum = 0;}
};


/// This class holds a list of the above objects.  We need one of
/// these objects for each section.  The () operator takes the name of
/// the section as a string and returns an integer corresponding to
/// the present number we should be using.
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



/// This class quite simply holds the HDF5 group ID and the
/// SectionNumber object for a section.
class HDF5SectionClass
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
  stack<HDF5SectionClass> SectionStack;
public:
  bool OpenFile (string fileName);
  void OpenSection (string name);
  void CloseSection ();
  void WriteVar (string name, double T);
  void WriteVar(string name, Array<double,1> &v);
  void WriteVar(string name, Array<double,2> &v);
  void WriteVar(string name, Array<double,3> &v);
  //void WriteVar (string name, Array<double,2> &m);
  void CloseFile();
  OutputSectionHDF5Class() : IsOpen(false)
  { }
};


#endif
