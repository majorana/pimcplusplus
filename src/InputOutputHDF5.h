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
  Array<hsize_t,1> Dimensions;

  bool ReadInto (double &val);
  bool ReadInto (Array<double,1> &v);
  bool ReadInto (Array<double,2> &v);
  bool ReadInto (Array<double,3> &v);
  bool ReadInto (int &val);
  bool ReadInto (Array<int,1> &v);
  bool ReadInto (Array<int,2> &v);
  bool ReadInto (Array<int,3> &v);
  bool ReadInto (string &val) {return true; }

  ~VarHDF5Class()
  {
    H5Dclose(DataSetID);
  }
};


/// This class stores a section of an HDF5 file.  The boolean value,
/// IsRoot, store whether this particular section is a the root node
/// of an HDF5 file.
class InputSectionHDF5Class : public InputSectionClass
{
private:
  bool IsOpen, IsRoot;
  /// ReadGroup reads a HDF5 group, given by name, from the file.
  /// It reads in all variables and groups within the file, calling
  /// itself recursively for groups within itself.
  void ReadGroup (hid_t parentGroupID, string name, InputSectionClass *parent);
  /// StripName strips the trailing ".#" from a string.  These were
  /// added by the HDF5 writer in order to have multiples sections
  /// with the same name.
  string StripName (string name);
  void PrintTree(int numIndent );
public:
  /// This is the HDF5 handle for the group.
  hid_t GroupID;
  /// This prints the variables and sections below me, mostly for
  /// debugging purposes.

  void PrintTree();
  void GroupIterator (string member_name);
  bool OpenFile (string fileName, string mySectionName,
		 InputSectionClass *parent);
  void Close();
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
