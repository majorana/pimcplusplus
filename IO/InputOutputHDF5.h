#ifndef INPUT_OUTPUT_HDF5_H
#define INPUT_OUTPUT_HDF5_H
#include "InputOutputBase.h"
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
  hid_t DataSetID, DataSpaceID;
  Array<hsize_t,1> Dimensions;

  bool ReadInto (double &val);
  bool ReadInto (Array<double,1> &v);
  bool ReadInto (Array<double,2> &v);
  bool ReadInto (Array<double,3> &v);
  bool ReadInto (int &val);
  bool ReadInto (Array<int,1> &v);
  bool ReadInto (Array<int,2> &v);
  bool ReadInto (Array<int,3> &v);
  bool ReadInto (string &val);
  bool ReadInto (Array<string,1> &val);
  bool ReadInto (Array<string,2> &val);
  bool ReadInto (Array<string,3> &val);
  /*  bool ReadInto (bool &val);
  bool ReadInto (Array<bool,1> &val);
  bool ReadInto (Array<bool,2> &val);
  bool ReadInto (Array<bool,3> &val); */
  bool Append (double val);
  bool Append (Array<double,1> &val);
  bool Append (Array<double,2> &val);
  bool Append (int val);
  bool Append (Array<int,1> &val);
  bool Append (Array<int,2> &val);
  bool Append (string val);
  bool Append (Array<string,1> &strs);
  bool Append (Array<string,2> &strs);
//   bool Append (int val);
//   bool Append (Array<int,1> val);
//   bool Append (Array<int,2 val);
//   bool Append (string val);
//   bool Append (Array<string,1> val);
//   bool Append (Array<string,2 val);

  ~VarHDF5Class()
  {
    H5Dclose(DataSetID);
    H5Sclose(DataSpaceID);
  }
};


/// This class stores a section of an HDF5 file.  The boolean value,
/// IsRoot, store whether this particular section is a the root node
/// of an HDF5 file.
class IOTreeHDF5Class : public IOTreeClass
{
private:
  bool IsOpen;
  /// ReadGroup reads a HDF5 group, given by name, from the file.
  /// It reads in all variables and groups within the file, calling
  /// itself recursively for groups within itself.
  void ReadGroup (hid_t parentGroupID, string name, IOTreeClass *parent);
  /// StripName strips the trailing ".#" from a string.  These were
  /// added by the HDF5 writer in order to have multiples sections
  /// with the same name.

  void PrintTree(int numIndent );
  void StripName (string str,string &newString,
		  int &myInt);
public:
  /// This is the HDF5 handle for the group.
  hid_t GroupID;
  IOTreeClass* NewSection(string name);
  /// This prints the variables and sections below me, mostly for
  /// debugging purposes.
  
  void PrintTree();
  void GroupIterator (string member_name);
  bool OpenFile (string fileName, string mySectionName,
		 IOTreeClass *parent);
  bool NewFile(string fileName,string myName,IOTreeClass* parent);
  void IncludeSection (IOTreeClass *);
  void CloseFile();
  void FlushFile();
  IOTreeHDF5Class() : IOTreeClass()
  {
    IsOpen=false;
    CurrSecNum=0;
  }


  void WriteVar (string name, double val);
  void WriteVar (string name, Array<double,1> &v);
  void WriteVar (string name, Array<double,2> &v);
  void WriteVar (string name, Array<double,3> &v);
  void WriteVar (string name, int val);
  void WriteVar (string name, Array<int,1> &v);
  void WriteVar (string name, Array<int,2> &v);
  void WriteVar (string name, Array<int,3> &v);
  void WriteVar (string name, string str);
  void WriteVar (string name, Array<string,1> &v);
  void WriteVar (string name, Array<string,2> &v);
  void WriteVar (string name, Array<string,3> &v);
};





/************************************************************
 *                    Output Functions                      *
 ************************************************************/





#endif
