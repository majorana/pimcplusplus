#ifndef INPUT_OUTPUT_ASCII_H
#define INPUT_OUTPUT_ASCII_H


#include "InputOutputBase.h"



class InputSectionASCIIClass : public InputSectionClass
{
  void ReadWithoutComments(string fileName, Array<char,1> &buffer);
  void PrintTree();
  void PrintTree(int num);
  bool ReadSection(string sectionName,
					 InputSectionClass* parent,
					   bool findBrace=true)
    ;

 public:
  //mySectionName isn't actually implemented here at all
  bool OpenFile (string fileName, string mySectionName,
		 InputSectionClass *parent);

  void CloseFile();
  ///Notice this doesn't do anything yet
  void Close(){; }
  
 
};


class VarASCIIClass : public VarClass
{  
public:
  Array<char,1> Value;
  void* PtrValue;



  bool ReadInto (Array<double,1> &val) {return true;}
  bool ReadInto (Array<double,2> &val) {return true;}
  bool ReadInto (Array<double,3> &val) {return true;}
  bool ReadInto (Array<int,1> &val)    {return true;}
  bool ReadInto (Array<int,2> &val)    {return true;}    
  bool ReadInto (Array<int,3> &val)    {return true;}
  //  virtual bool ReadInto (string &val) {return true;}
  bool ReadInto (Array<string,1> &val){return true;}
  bool ReadInto (Array<string,2> &val){return true;}
  bool ReadInto (Array<string,3> &val) {return true;}



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
  bool ReadInto (int &val){
    int l = Value.size();
    char str[l+1];
    for (int i=0; i<l; i++)
      str[i] = Value(i);
    str[l] ='\0';
    cerr<<"Value is "<<Value<<endl;
    val = atoi(str);
    cerr<<"val is "<<val<<endl;
    return true;
  }
  bool ReadInto (string &val){
    val="";
    for (int counter=1;counter<Value.size()-1;counter++){
      val=val+Value(counter);
    }
    return true;
  }
};

#endif
