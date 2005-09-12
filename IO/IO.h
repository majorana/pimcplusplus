#ifndef IO_H
#define IO_H

#include "IOVar.h"
#include "IOHDF5.h"
#include "IOASCII.h"

namespace IO {

  template<typename T> bool 
  IOTreeClass::WriteVar (string name, T val)
  {
    if (GetFileType() == HDF5_TYPE) {
      IOTreeHDF5Class *h5Tree = dynamic_cast<IOTreeHDF5Class*>(this); 
      if (h5Tree == NULL) {
	cerr << "Error in dynamic_cast in WriteVar.\n";
	abort();
      }
      hid_t groupID = h5Tree->GroupID;
      hid_t boolType = h5Tree->GetBoolType();
      IOVarBase *newVar = NewIOVarHDF5(groupID, name, val, boolType);
      VarList.push_back(newVar);
    }
    else if (GetFileType() == ASCII_TYPE) {
      
    }
    else {
      cerr << "Unknown file type in WriteVar.\n";
      abort();
    }
  }

  template<typename T, int RANK> bool 
  IOTreeClass::WriteVar (string name, const Array<T,RANK> &val)
  {
    if (GetFileType() == HDF5_TYPE) {
      IOTreeHDF5Class *h5Tree = dynamic_cast<IOTreeHDF5Class*>(this);
      if (h5Tree == NULL) {
	cerr << "Error in dynamic_cast in WriteVar.\n";
	abort();
      }
      hid_t groupID = h5Tree->GroupID;
      hid_t boolType = h5Tree->GetBoolType();
      IOVarBase *newVar = NewIOVarHDF5(groupID, name, val, boolType);
      VarList.push_back(newVar);
    }
    else if (GetFileType() == ASCII_TYPE) {
      VarList.push_back(new IOVarASCII<T,RANK>(name, val));
    }
    else {
      cerr << "Unknown file type in WriteVar.\n";
      abort();
    }
  }


  /// In the file name format name.extn, returns the extension.
  /// Actually returns everything after the trailing.
  inline string Extension (string fileName);


  /// This function takes a filename, determines it extension, creates a
  /// new IOTreeASCIIClass or IOTreeHDF5Class based on the
  /// extension, and calls OpenFile on the new object.
  /// Extensions:  
  /// .h5:            HDF5
  /// .xml:           XML
  /// .anything_else  ASCII
  IOTreeClass *ReadTree (string fileName, string myName, IOTreeClass *parent);

  IOTreeClass *NewTree (string fileName, string myName, IOTreeClass *parent);



  ///  Wrapper class for IOTreeClass that gives a nearly identical
  ///  interface as the OutputSectionClass.
  class IOSectionClass
  {
  private:
    IOTreeClass *CurrentSection;
  public:

    /// Opens the file reference by fileName and reads the contents into
    /// the tree in CurrentSection.  Creates a new object based on the
    /// extnesion of the filename.  For ".h5", it creates an
    /// IOTreeHDF5Class.  For ".xml" it creaes an IOTreeXMLClass.
    /// After creating the object, it calls the objects virtual OpenFile
    /// function, reading the contents of the file into the tree.
    bool OpenFile (string fileName);
    string GetName(){ return CurrentSection->Name;}
    inline string GetFileName();
    string GetVarName(int num){ return GetVarPtr(num)->GetName();}
    /// Creates a file at the top level, choosing the appropriate type
    /// based on the file extension.
    bool NewFile (string fileName);

    /// Calls CurrentSections close file and then deletes the
    /// CurrentSection.  
    void CloseFile ();

    /// Flush all buffers to disk for safety
    void FlushFile();

    /// Opens the num'th section with the given name.  The default
    /// value for num is 0.
    bool OpenSection (string name, int num=0);

    /// Opens the num'th section below CurrentSection.
    bool OpenSection (int num);

    /// This mounts a file in the current tree under CurrentSection at
    /// the end of CurrentsSection's SectionList.  It does not change
    /// what CurrentSection points to, ie. it does not descend to the
    /// newly-opened section.
    bool IncludeSection (string name, string fileName);

    /// Creates a new section of the same type as currentSection under
    /// currentSection.  Pushes the new section to the end of the
    /// section list.
    inline void NewSection (string name)
    {  CurrentSection = CurrentSection->NewSection(name); }

    /// This function creates a new file of the appropriate type as
    /// determined by the extension of fileName and mounts it at the end
    /// of the list under CurrentSection.  Returns false if the file
    /// couldn't be created.
    bool NewSection (string name, string fileName);

    /// Closes the current section.  That is, CurrentSection becomes
    /// CurrentSection's parent.
    void CloseSection ();

    /// Template function which reads a variable in the present section
    /// into the passed-by-reference T variable.
    template<class T>
    bool ReadVar(string name, T &var)
    {  return (CurrentSection->ReadVar(name, var)); }

    template<class T>
    bool ReadVar(string name, T &var, T Default)
    { 
      bool success = ReadVar(name, var);
      if (!success)
	var = Default;
      return (success);
    }

    /// Writes a variable under the current section.
    template<typename T> bool
    WriteVar (string name, T val)
    { return CurrentSection->WriteVar(name, val); }

    template<typename T, int RANK> bool
    WriteVar (string name, Array<T,RANK>& val)
    { return CurrentSection->WriteVar(name, val); }
  
    template<class T>
    bool AppendVar(string name, T val)
    { return CurrentSection->AppendVar(name, val); }
  
    inline IOVarBase *GetVarPtr(string name)
    {    return (CurrentSection->GetVarPtr(name)); }

    inline IOVarBase *GetVarPtr(int num)
    { return (CurrentSection->GetVarPtr(num)); }



    /// Returns the number of subsections within the present section
    /// which have the name name.  If called without a name, it returns
    /// the total number of sections.
    inline int CountSections(string name="")
    { return (CurrentSection->CountSections(name)); }
    inline int CountVars()
    {return (CurrentSection->CountVars());}
    /// Calls CurrentSections virtual PrintTree() function.  This is for
    /// debugging purposes.  It spits out a hierarchy of the sections
    /// and variable names.
    void PrintTree()
    { CurrentSection->PrintTree(); }

    IOSectionClass() 
    {
      CurrentSection=NULL;
    }
	
  };


} // Ends namespace IO

#endif
