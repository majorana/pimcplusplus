#ifndef INPUT_OUTPUT_HDF5_H
#include "InputOutput.h"
#include <iostream>
#include <stack>
#include "hdf5.h"


#ifndef H5_NO_NAMESPACE
using namespace H5;
#endif


class VarHDF5Class : public VarClass
{
public:
  bool ReadInto (double &val) { }
  bool ReadInto (int &val) { }
  bool ReadInto (string &val) { }
};


class OutputSectionHDF5Class : public OutputSectionClass
{
private:
  hid_t FileID;
  string FileName;
  bool IsOpen;
  stack<hid_t> GroupStack;
public:
  bool OpenFile (string fileName);
  void OpenSection (string name);
  void CloseSection ();
  void WriteVar (string name, double &T);
  void CloseFile();
  OutputSectionClass() : IsOpen(false)
  { }
}


#endif
