#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include "InputOutputBase.h"
#include "InputOutputHDF5.h"
//#include "InputOutputASCII.h"

#include <stack>
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
  bool OpenFile (string fileName)
  {
    bool success;
    string extn = Extension (fileName);
    if (extn == "h5") 
      CurrentSection = new InputTreeHDF5Class;
    //else if (extn == "xml")
    //CurrentSection = new InputTreeXMLClass;
    //else
    //  CurrentSection = new InputTreeASCIIClass;
    success = CurrentSection->OpenFile(fileName, "Root", NULL);
    if (!success)
      delete (CurrentSection);
    
    return (success);
  }

  void CloseFile ()
  {
    CurrentSection->CloseFile();
    delete (CurrentSection);
  }

  inline bool OpenSection (string name)
  {
    InputTreeClass *newSection;
    bool success = CurrentSection->FindSection(name, newSection);
    if (success)
      CurrentSection=newSection;
    return success;
  }
  
  inline void CloseSection ()
  {
    assert (CurrentSection->Parent != NULL);
    CurrentSection = CurrentSection->Parent;
  }

  template<class T>
  bool ReadVar(string name, T &var)
  {
    CurrentSection->ReadVar(name, var);
  }

  inline int CountSections(string fileName)
  {
    CurrentSection->CountSections(fileName);
  }
};





#endif
