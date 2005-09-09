#ifndef INPUT_OUTPUT_HDF5_H
#define INPUT_OUTPUT_HDF5_H
#include "IOBase.h"
#include <iostream>
#include <stack>
#include <hdf5.h>


namespace IO {

  /// This class stores a section of an HDF5 file.  The boolean value,
  /// IsRoot, store whether this particular section is a the root node
  /// of an HDF5 file.
  class IOTreeHDF5Class : public IOTreeClass
  {
  private:
    bool IsOpen;
    hid_t BoolType;
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
  
    IOFileType GetType();
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
  };

}


#endif
