#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include "InputOutputBase.h"
#include "InputOutputHDF5.h"
// #include "InputOutputASCII.h"

#include <stack>


/// In the file name format name.extn, returns the extension.
/// Actually returns everything after the trailing.
inline string Extension (string fileName)
{
  string extn;
  stack<char> bwExtn;
  int pos = fileName.length()-1;
  while ((pos >= 0) && fileName[pos]!='.') {
    bwExtn.push(fileName[pos]);
    pos--;
  }
  
  if (fileName[pos] == '.') 
    while (!bwExtn.empty()) {
      extn += bwExtn.top();
      bwExtn.pop();
    }
  else
    extn = "";
  return (extn);
}


/// This function takes a filename, determines it extension, creates a
/// new IOTreeASCIIClass or IOTreeHDF5Class based on the
/// extension, and calls OpenFile on the new object.
/// Extensions:  
/// .h5:            HDF5
/// .xml:           XML
/// .anything_else  ASCII
inline IOTreeClass *ReadTree (string fileName, 
			      string myName,
			      IOTreeClass *parent)
{
  IOTreeClass *newTree;
  string extn = Extension (fileName);
  if (extn == "h5")
    newTree = new IOTreeHDF5Class;
  //  else if (extn == "xml")
  //    newTree = newIOTreeXMLClass;
  /////else
  /////    newTree = new IOTreeASCIIClass;
  
  newTree->FileName = fileName;
  bool success = newTree->OpenFile (fileName, myName, parent);
  if (success)
    return (newTree);
  else{
    delete newTree;
    return (NULL);
  }
}


inline IOTreeClass *NewTree (string fileName,
			     string myName,
			     IOTreeClass *parent)
{
  IOTreeClass *newTree;
  string extn = Extension (fileName);
  if (extn == "h5")
    newTree = new IOTreeHDF5Class;
  //  else if (extn == "xml")
  //    newTree = newIOTreeXMLClass;
  ////// else
  /////    newTree = new IOTreeASCIIClass;

  bool success = newTree->NewFile (fileName, myName, parent);
  if (success)
    return (newTree);
  else{
    delete newTree;
    return (NULL);
  }
}


///  Wrapper class for IOTreeClass that gives a nearly identical
///  interface as the OutputSectionClass.
class IOSectionClass
{
private:
  IOTreeClass *CurrentSection;
  bool IsModified;
public:

  /// Opens the file reference by fileName and reads the contents into
  /// the tree in CurrentSection.  Creates a new object based on the
  /// extnesion of the filename.  For ".h5", it creates an
  /// IOTreeHDF5Class.  For ".xml" it creaes an IOTreeXMLClass.
  /// After creating the object, it calls the objects virtual OpenFile
  /// function, reading the contents of the file into the tree.
  bool OpenFile (string fileName)
  {
    CurrentSection = ReadTree (fileName, "Root", NULL);
    if (CurrentSection == NULL)
      return (false);
    else
      return (true);
  }

  /// Creates a file at the top level, choosing the appropriate type
  /// based on the file extension.
  bool NewFile (string fileName)
  {
    CurrentSection = NewTree (fileName, "Root", NULL);
    if (CurrentSection == NULL)
      return (false);
    else
      return (true);
  }

  /// Calls CurrentSections close file and then deletes the
  /// CurrentSection.  
  void CloseFile ()
  {
    while (CurrentSection->Parent != NULL)
      CloseSection();
    CurrentSection->CloseFile();
    delete (CurrentSection);
  }

  /// Opens the num'th section with the given name.  The default
  /// value for num is 0.
  inline bool OpenSection (string name, int num=0)
  {
    IOTreeClass *newSection;
    bool success;
    success = CurrentSection->FindSection(name, newSection, num);
    if (success)
      CurrentSection=newSection;
    return success;
  }

  /// Opens the num'th section below CurrentSection.
  inline bool OpenSection (int num)
  {
    IOTreeClass *newSection;
    list<IOTreeClass*>::iterator Iter=CurrentSection->SectionList.begin();
    int i = 0;
    while ((i<num) && 
	   (Iter != CurrentSection->SectionList.end())){
      i++;
      Iter++;
    }
    if (i<num)
      return false;
    else {
      CurrentSection = *Iter;
      return true;
    }
  }

  /// This mounts a file in the current tree under CurrentSection at
  /// the end of CurrentsSection's SectionList.  It does not change
  /// what CurrentSection points to, ie. it does not descend to the
  /// newly-opened section.
  inline bool IncludeSection (string name, string fileName)
  {
    IOTreeClass *newSection;
    newSection = ReadTree (fileName, name, CurrentSection);
    if (newSection == NULL)
      return false;
    else {
      CurrentSection->IncludeSection(newSection);
      return true;
    }
  }

  /// Creates a new section of the same type as currentSection under
  /// currentSection.  Pushes the new section to the end of the
  /// section list.
  inline void NewSection (string name)
  {
    CurrentSection = CurrentSection->NewSection(name);
  }

  /// This function creates a new file of the appropriate type as
  /// determined by the extension of fileName and mounts it at the end
  /// of the list under CurrentSection.  Returns false if the file
  /// couldn't be created.
  inline bool NewSection (string name, string fileName)
  {
    IOTreeClass *newSection;
    newSection = NewTree (fileName, name, CurrentSection);
     if (newSection == NULL)
      return false;
    else {
      CurrentSection->IncludeSection(newSection);
      CurrentSection = newSection;
      return true;
    }
  }

  /// Closes the current section.  That is, CurrentSection becomes
  /// CurrentSection's parent.
  inline void CloseSection ()
  {
    assert (CurrentSection->Parent != NULL);
    CurrentSection = CurrentSection->Parent;
  }

  /// Template function which reads a variable in the present section
  /// into the passed-by-reference T variable.
  template<class T>
  bool ReadVar(string name, T &var)
  {
    return (CurrentSection->ReadVar(name, var));
  }

  /// Writes a variable under the current section.
  void WriteVar(string name, double val)
  {
    CurrentSection->WriteVar(name, val);
  }
  void WriteVar(string name, Array<double,1> &val)
  {
    CurrentSection->WriteVar(name, val);
  }
  void WriteVar(string name, Array<double,2> &val)
  {
    CurrentSection->WriteVar(name, val);
  }
  void WriteVar(string name, Array<double,3> &val)
  {
    CurrentSection->WriteVar(name, val);
  }
  bool AppendVar(string name, double val)
  {
    return CurrentSection->AppendVar(name, val);
  }
  bool AppendVar(string name, Array<double,1> &val)
  {
    return CurrentSection->AppendVar(name, val);
  }



  void WriteVar(string name, int val)
  {
    CurrentSection->WriteVar(name, val);
  }
  void WriteVar(string name, Array<int,1> &val)
  {
    CurrentSection->WriteVar(name, val);
  }
  void WriteVar(string name, Array<int,2> &val)
  {
    CurrentSection->WriteVar(name, val);
  }
  void WriteVar(string name, Array<int,3> &val)
  {
    CurrentSection->WriteVar(name, val);
  }


//   void WriteVar(string name, bool val)
//   {
//     CurrentSection->WriteVar(name, val);
//   }
//   void WriteVar(string name, Array<bool,1> &val)
//   {
//     CurrentSection->WriteVar(name, val);
//   }
//   void WriteVar(string name, Array<bool,2> &val)
//   {
//     CurrentSection->WriteVar(name, val);
//   }
//   void WriteVar(string name, Array<bool,3> &val)
//   {
//     CurrentSection->WriteVar(name, val);
//   }


  void WriteVar(string name, string val)
  {
    CurrentSection->WriteVar(name, val);
  }
  void WriteVar(string name, Array<string,1> &val)
  {
    CurrentSection->WriteVar(name, val);
  }
  void WriteVar(string name, Array<string,2> &val)
  {
    CurrentSection->WriteVar(name, val);
  }
  void WriteVar(string name, Array<string,3> &val)
  {
    CurrentSection->WriteVar(name, val);
  }

  /// Returns the number of subsections within the present section
  /// which have the name name.  If called without a name, it returns
  /// the total number of sections.
  inline int CountSections(string name="")
  {
    return (CurrentSection->CountSections(name));
  }

  /// Calls CurrentSections virtual PrintTree() function.  This is for
  /// debugging purposes.  It spits out a hierarchy of the sections
  /// and variable names.
  void PrintTree()
  {
    CurrentSection->PrintTree();
  }

  IOSectionClass() 
  {
    CurrentSection=NULL;
    IsModified=false;
  }
	
};





#endif
