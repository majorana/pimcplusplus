#ifndef INPUT_OUTPUT_H
#define INPUT_OUTPUT_H

#include <string>
#include <list>
#include "Blitz.h"
#include <fstream>

using namespace std;



class VarClass 
{  
public:
  string Name;

  virtual bool ReadInto (double &val) = 0;
  virtual bool ReadInto (int &val)                 = 0;
  virtual bool ReadInto (string &val)              = 0;
  /*  virtual bool ReadInto (Array<double,1> &val) = 0;
  virtual bool ReadInto (Array<double,2> &val)     = 0;
  virtual bool ReadInto (Array<double,3> &val)     = 0;
  virtual bool ReadInto (int &val)                 = 0;
  virtual bool ReadInto (Array<int,1> &val)        = 0;
  virtual bool ReadInto (Array<int,2> &val)        = 0;
  virtual bool ReadInto (Array<int,3> &val)        = 0;
  virtual bool ReadInto (bool &val)                = 0;
  virtual bool ReadInto (Array<bool,1> &val)       = 0;
  virtual bool ReadInto (Array<bool,2> &val)       = 0;
  virtual bool ReadInto (Array<bool,3> &val)       = 0;

  virtual bool ReadInto (Array<string,1> &val)     = 0;
  virtual bool ReadInto (Array<string,2> &val)     = 0;
  virtual bool ReadInto (Array<string,3> &val)     = 0; */
};





class InputSectionClass
{
 public:
  virtual void PrintTree(InputSectionClass *sec)=0;
  InputSectionClass* Parent;
  list<InputSectionClass*> SectionList;
  list<VarClass*> VarList;
  list<InputSectionClass*>::iterator Iter;
public:
  string Name;
  inline bool FindSection (string name, InputSectionClass * &sectionPtr, 
		    bool rewind=true);
  inline void Rewind(){
    Iter=SectionList.begin();
  }
  inline int CountSections(string name);
  template<class T> bool ReadVar(string name, T &var);
  virtual bool OpenFile (string fileName, InputSectionClass *parent) = 0;
  virtual void CloseFile() = 0;

};



inline int InputSectionClass::CountSections(string name)
{
  list<InputSectionClass*>::iterator sectionIter;
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






inline bool InputSectionClass::FindSection (string name, 
				     InputSectionClass* &sectionPtr,
				     bool rewind)
{
  
  list<InputSectionClass*>::iterator tempIter=Iter;
  while ((Iter != SectionList.end()) && ((*Iter)->Name!=name)){
    Iter++;
  }
  bool found = Iter != SectionList.end(); 
  if (found){
    sectionPtr = *Iter;
    Iter++;
  }
  if (rewind)
    Iter = SectionList.begin();
  sectionPtr->Iter=sectionPtr->SectionList.begin();
  if (!found){
    Iter=tempIter;
  }
  return (found);
}




template <class vartype>
class SpecialInputSectionClass : public InputSectionClass
{
  void ReadWithoutComments(string fileName,Array<char,1> &buffer);
  void printTree(InputSectionClass *sec);
  
public:
  bool OpenFile (string fileName);
  void CloseFile();
};




template<class T>
bool InputSectionClass::ReadVar(string name, T &var)
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





class OutputSectionClass
{
public:
  virtual bool OpenFile(string fileName)=0;
  virtual void OpenSection(string name)=0;
  virtual void CloseSection()=0;
  virtual void WriteVar(string name, double T)=0;
  virtual void WriteVar(string name, Array<double,1> &v) = 0;
  virtual void WriteVar(string name, Array<double,2> &v) = 0;
  virtual void WriteVar(string name, Array<double,3> &v) = 0;
};




//class Object
//{
//public:
//  virtual string Name = 0;
///  virtual Read (InputSectionClass &section) = 0;
// virtual CheckPoint(OutputSectionClass &outSection);
//};



/* private: */
/*   list<InputSectionClass*> MyChildren; */
/*   InputSection *MyParent; */
/* public: */
/*   virtual void Open (char *FileName) = 0; */
/*   void Open (string FileName) { Open (FileName.c_str());} */

/*   InputSectionClass& FindSection (string Name, bool rewind=true); */
/*   /// Atomics */
/*   virtual bool ReadVar (string name, int var)    = 0; */
/*   virtual bool ReadVar (string name, double var) = 0; */
/*   virtual bool ReadVar (string name, bool var)   = 0; */
/*   virtual bool ReadVar (string name, string var) = 0; */
/*   /// 1D Arrays */
/*   virtual bool ReadVar (string name, Array<int,1> var)    = 0; */
/*   virtual bool ReadVar (string name, Array<double,1> var) = 0; */
/*   virtual bool ReadVar (string name, Array<bool,1> var)   = 0; */
/*   virtual bool ReadVar (string name, Array<bool,1> var)   = 0; */
/*   /// 1D Arrays */
/*   virtual bool ReadVar (string name, Array<int,2> var)    = 0; */
/*   virtual bool ReadVar (string name, Array<double,2> var) = 0; */
/*   virtual bool ReadVar (string name, Array<bool,2> var)   = 0; */
/*   virtual bool ReadVar (string name, Array<bool,2> var)   = 0; */

/*   template<class T> bool ReadVar (string name, Array<T,1> var); */



#endif
