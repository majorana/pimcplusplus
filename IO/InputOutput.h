#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include "InputOutputBase.h"
#include "InputOutputHDF5.h"
//#include "InputOutputASCII.h"

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




///  Wrapper class for InputTreeClass that gives a nearly identical
///  interface as the OutputSectionClass.
class InputSectionClass
{
private:
  InputTreeClass *CurrentSection;
public:

  /// Opens the file reference by fileName and reads the contents into
  /// the tree in CurrentSection.  Creates a new object based on the
  /// extnesion of the filename.  For ".h5", it creates an
  /// InputTreeHDF5Class.  For ".xml" it creaes an InputTreeXMLClass.
  /// After creating the object, it calls the objects virtual OpenFile
  /// function, reading the contents of the file into the tree.
  bool OpenFile (string fileName)
  {
    bool success;
    string extn = Extension (fileName);
    if (extn == "h5") 
      CurrentSection = new InputTreeHDF5Class;
    else {
      cerr << "Unrecognized extension " << extn << endl;
      exit(1);
    }
    //else if (extn == "xml")
    //CurrentSection = new InputTreeXMLClass;
    //else
    //  CurrentSection = new InputTreeASCIIClass;
    success = CurrentSection->OpenFile(fileName, "Root", NULL);
    if (!success)
      delete (CurrentSection);
    
    return (success);
  }

  /// Calls CurrentSections virtual PrintTree() function.  This is for
  /// debugging purposes.  It spits out a hierarchy of the sections
  /// and variable names.
  void PrintTree()
  {
    CurrentSection->PrintTree();
  }

  /// Calls CurrentSections close file and then deletes the
  /// CurrentSection.  
  void CloseFile ()
  {
    CurrentSection->CloseFile();
    delete (CurrentSection);
  }

  /// Opens the num'th section with the given name.  The default
  /// value for num is 0.
  inline bool OpenSection (string name, int num=0)
  {
    InputTreeClass *newSection;
    CurrentSection->Rewind();
    bool success;
    for (int i=0; i<=num; i++)
      success = CurrentSection->FindSection(name, newSection, false);
    if (success)
      CurrentSection=newSection;
    return success;
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

  /// Returns the number of subsections within the present section
  /// which have the name name.
  inline int CountSections(string name)
  {
    return (CurrentSection->CountSections(name));
  }
};





#endif
