#include <fstream>
#include <string.h>
#include <iostream>
#include "Blitz.h"
#include <list>

using namespace std;




class VarClass 
{  
public:
  string Name;

  virtual bool ReadInto (double &val) = 0;
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
  virtual bool ReadInto (string &val)              = 0;
  virtual bool ReadInto (Array<string,1> &val)     = 0;
  virtual bool ReadInto (Array<string,2> &val)     = 0;
  virtual bool ReadInto (Array<string,3> &val)     = 0; */
};



class VarASCIIClass : public VarClass
{  
  Array<char,1> Value;
public:
  bool ReadInto (double &val)
  {
    int l = Value.size();
    char str[l+1];
    for (int i=0; i<l; i++)
      str[i] = Value(i);
    str[l] ='\0';
    val = atof(str);
    return true;
  }
  /*  bool ReadInto (Array<double,1> &val);
  bool ReadInto (Array<double,2> &val);
  bool ReadInto (Array<double,3> &val);
  bool ReadInto (int &val);
  bool ReadInto (Array<int,1> &val);
  bool ReadInto (Array<int,2> &val);
  bool ReadInto (Array<int,3> &val);
  bool ReadInto (bool &val);
  bool ReadInto (Array<bool,1> &val);
  bool ReadInto (Array<bool,2> &val);
  bool ReadInto (Array<bool,3> &val);
  bool ReadInto (string &val);
  bool ReadInto (Array<string,1> &val);
  bool ReadInto (Array<string,2> &val);
  bool ReadInto (Array<string,3> &val); */
};



class SectionTreeClass
{
protected:
  SectionTreeClass* parent;
  list<SectionTreeClass*> SectionList;
  list<VarClass*> VarList;
  list<SectionTreeClass*>::iterator Iter;
public:
  string Name;
  bool FindSection (string name, SectionTreeClass * &sectionPtr, 
		    bool rewind=true);
  //  virtual template<class T> bool ReadVar(string name, T &var) = 0;
  virtual bool ReadVar (string name, double &var);
  virtual bool ReadFile (string fileName) = 0;
};

bool SectionTreeClass::FindSection (string name, SectionTreeClass* &sectionPtr,
				    bool rewind)
{
  while ((Iter != SectionList.end()) && ((*Iter)->Name!=name))
    Iter++;
  bool found = Iter != SectionList.end(); 
  if (found)
    sectionPtr = *Iter;
  if (rewind)
    Iter = SectionList.begin();
  return (found);
}


template <class vartype>
class SpecialSectionClass : public SectionTreeClass
{
public:
  bool ReadFile (string fileName);
  template<class T> bool ReadVar(string name, T &var)
  {
    bool readVarSuccess;
    list<VarClass*>::iterator varIter=VarList.begin();
    while ((varIter!=VarList.end() && varIter->Name!=name)){
      varIter++;
    }
    bool found = varIter != VarList.end();
    if (found)
      readVarSuccess=varIter.ReadInto(var);
    else if (parent!=NULL){
      readVarSuccess=parent->ReadVar(name,var);
    }
    else {
      cerr<<"Couldn't find variable "<<name;
      return false;
    }  
    return readVarSuccess;	 
  }
};





bool SpecialSectionClass<VarASCIIClass>::ReadFile(string FileName)
{
  return true;

}




inline bool checkPair(Array<char,1> &buffer,int counter,char* toSee)
{
  if (counter+1>=buffer.size()){
    return false;
  }
  if (buffer(counter)==toSee[0] && buffer(counter+1)==toSee[1]){
    return true;
  }
  else return false;

}



// void ReadWithoutComments(string fileName,Array<char,1> &buffer)
// {
//   ifstream infile;
//   infile.open(fileName.c_str());
//   Array<char,1> tmpBuffer;
//   int counter=0;
//   bool inQuote=false;
//   char dummyChar;
//   while (!infile.eof()){
//     infile.get(dummyChar);    
//     counter++;
//   }

//   tmpBuffer.resize(counter);
//   buffer.resize(counter);
//   counter=-1;
//   infile.close();
//   ifstream infile2;
//   infile2.open(fileName.c_str());
//   while (!infile2.eof()){
//     counter++;
//     infile2.get(dummyChar);
//     tmpBuffer(counter)=dummyChar;    
//   }
//   cout<<tmpBuffer;
//   int bufferLoc=0;
//   for (int counter=0;counter<tmpBuffer.size();counter++){
//     if (inQuote){
//       while (tmpBuffer(counter) != '\"' && counter<tmpBuffer.size()){
// 	buffer(bufferLoc)=tmpBuffer(counter);
// 	counter++;
// 	bufferLoc++;
//       }
//       buffer(bufferLoc)=tmpBuffer(counter);
//       bufferLoc++;
//       inQuote=false;
//     }
//     else {
//       if (checkPair(tmpBuffer,counter,"//")){
// 	while (tmpBuffer(counter)!='\n' && counter<tmpBuffer.size()){
// 	  counter++;
// 	}
// 	buffer(bufferLoc)=tmpBuffer(counter); //copy the \n over
// 	bufferLoc++;
//       }
//       else if (checkPair(tmpBuffer,counter,"/*")){
// 	while (!checkPair(tmpBuffer,counter,"*/") && counter<tmpBuffer.size()){
// 	  counter++;
// 	}
// 	counter++; //end up in the / of comment
//       }
//       else if (tmpBuffer(counter)=='\"'){
// 	inQuote=true;
// 	buffer(bufferLoc)=tmpBuffer(counter);
// 	bufferLoc++;
//       }
//       else {
// 	buffer(bufferLoc)=tmpBuffer(counter);
// 	bufferLoc++;
//       }
//     }
//   }
//   buffer.resizeAndPreserve(bufferLoc+1);
// }



// int main()
// {
//   Array<char,1> myBuffer;
//   ReadWithoutComments("tryFile",myBuffer);
//   for (int counter=0;counter<myBuffer.size();counter++){
//     cout<<myBuffer(counter);
//   }

// }
 

main()
{
  exit(1);
}
