#ifndef INPUT_OUTPUT_BASE_H
#define INPUT_OUTPUT_BASE_H

#include <string>
#include <list>
#include <stack>
#include "../Blitz.h"
#include <fstream>

using namespace std;

/// This typedef enumerates the atomic datatype that are supported in
/// inputfiles.   Note that 1,2, and 3-dimensional arrays of these
/// types are also supported
typedef enum{NOT_ATOMIC,INT_TYPE,DOUBLE_TYPE,STRING_TYPE,BOOL_TYPE} AtomicType; 

/// This is the abstract base class for the storage of variables in
/// the InputTreeClass lists.   It contains the name of the variable,
/// it's atomic type (see above), and it's dimensionality.
/// Specializations of this class will also have some reference to the
/// data in the variable itself.
class VarClass 
{  
public:
  string Name;
  AtomicType Type;
  int Dim;
  

  virtual bool ReadInto (double &val)              = 0;
  virtual bool ReadInto (Array<double,1> &val)     = 0;
  virtual bool ReadInto (Array<double,2> &val)     = 0;
  virtual bool ReadInto (Array<double,3> &val)     = 0;
  virtual bool ReadInto (int &val)                 = 0;
  virtual bool ReadInto (Array<int,1> &val)        = 0;
  virtual bool ReadInto (Array<int,2> &val)        = 0;
  virtual bool ReadInto (Array<int,3> &val)        = 0;
  virtual bool ReadInto (string &val)              = 0;
  virtual bool ReadInto (Array<string,1> &val)     = 0;
  virtual bool ReadInto (Array<string,2> &val)     = 0;
  virtual bool ReadInto (Array<string,3> &val)     = 0; 
  /*virtual bool ReadInto (bool &val)                = 0;
  virtual bool ReadInto (Array<bool,1> &val)       = 0;
  virtual bool ReadInto (Array<bool,2> &val)       = 0;
  virtual bool ReadInto (Array<bool,3> &val)       = 0;*/

};




/// This class stores a tree of input file sections.  Each member of
/// the tree contains a list of tree nodes below it and a list of
/// variables contained in the present node.
class InputTreeClass
{
 public:
  virtual void PrintTree(int numIndent)=0;
  virtual void PrintTree()=0;
  InputTreeClass* Parent;
  list<InputTreeClass*> SectionList;
  list<VarClass*> VarList;
  list<InputTreeClass*>::iterator Iter;
public:
  string Name;
  inline bool FindSection (string name, InputTreeClass * &sectionPtr, 
		    bool rewind=true);
  inline void Rewind(){
    Iter=SectionList.begin();
  }
  inline int CountSections(string name);

  template<class T>
    bool ReadVar(string name, T &var)
    {
      bool readVarSuccess;
      list<VarClass*>::iterator varIter=VarList.begin();
      while ((varIter!=VarList.end() && (*varIter)->Name!=name)){
	varIter++;
      }
      bool found = varIter != VarList.end();
      if (found){
	readVarSuccess=(*varIter)->ReadInto(var);
      }
      else if (Parent!=NULL){
	readVarSuccess=Parent->ReadVar(name,var);
      }
      else {
	cerr<<"Couldn't find variable "<<name;
	return false;
      }  
      return readVarSuccess;	 
    }
  
  virtual bool OpenFile (string fileName, 
			 string mySectionName, 
			 InputTreeClass *parent) = 0;
  virtual void CloseFile() = 0;
};



/// Returns the number of subsections with the given name within the
/// present section.
inline int InputTreeClass::CountSections(string name)
{
  list<InputTreeClass*>::iterator sectionIter;
  sectionIter=SectionList.begin();
  int numSections=0;
  while (sectionIter!=SectionList.end()){
    if (name==(*sectionIter)->Name){
      numSections++;
    }
    sectionIter++;
  }
  return numSections;
}





/// FindSection locates a subsection with the given name within the
/// section in contains and returns it in the pointer, sectionPtr,
/// which is passed a reference.  Returns true if the section is
/// found.  The final parameter, which default value "true",
/// optionally resets the section iterator to the beginning of the
/// section.  Thus, one may control whether or not order is
/// significant.  
inline bool InputTreeClass::FindSection (string name, 
					 InputTreeClass* &sectionPtr,
					 bool rewind)
{
  
  list<InputTreeClass*>::iterator tempIter=Iter;
  while ((Iter != SectionList.end()) && ((*Iter)->Name!=name)){
    Iter++;
  }
  bool found = Iter != SectionList.end(); 
  if (found){
    sectionPtr = *Iter;
    sectionPtr->Iter=sectionPtr->SectionList.begin();
    Iter++;
  }
  if (rewind)
    Iter = SectionList.begin();

  if (!found){
    Iter=tempIter;
  }
  return (found); 
}







class OutputSectionClass
{
public:
  virtual bool OpenFile(string fileName)                 = 0;
  virtual void CloseFile()                               = 0;
  virtual void OpenSection(string name)                  = 0;
  virtual void CloseSection()                            = 0;
  virtual void WriteVar(string name, double T)           = 0;
  virtual void WriteVar(string name, Array<double,1> &v) = 0;
  virtual void WriteVar(string name, Array<double,2> &v) = 0;
  virtual void WriteVar(string name, Array<double,3> &v) = 0;
  virtual void WriteVar(string name, int T)              = 0;
  virtual void WriteVar(string name, Array<int,1> &v)    = 0;
  virtual void WriteVar(string name, Array<int,2> &v)    = 0;
  virtual void WriteVar(string name, Array<int,3> &v)    = 0;
  virtual void WriteVar(string name, string str)         = 0;
  virtual void WriteVar(string name, Array<string,1> &v) = 0;
  virtual void WriteVar(string name, Array<string,2> &v) = 0;
  virtual void WriteVar(string name, Array<string,3> &v) = 0;

};




#endif
