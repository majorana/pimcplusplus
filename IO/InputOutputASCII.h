#ifndef INPUT_OUTPUT_ASCII_H
#define INPUT_OUTPUT_ASCII_H
#include "InputOutputBase.h"
#include <iostream>
#include <stack>
#include <string>
#include <list>



class VarASCIIClass : public VarClass
{  
public:
  /// Hack to hold different types.  We cast to to whatever type is
  /// indicated by the Type and Dim variables in the parent class.
  void *Value;

  /// The ReadIntos take copy their contents into val if the type
  /// that's passed to them matches the type that's actually stored in
  /// it's void *Value, as indicated by Type and Dim.
  bool ReadInto (double &val);
  bool ReadInto (int &val);
  bool ReadInto (string &val);
  bool ReadInto (bool &val);
  bool ReadInto (Array<double,1> &v);
  bool ReadInto (Array<double,2> &v);
  bool ReadInto (Array<double,3> &v);
  bool ReadInto (Array<int,1> &v);
  bool ReadInto (Array<int,2> &v);
  bool ReadInto (Array<int,3> &v);
  bool ReadInto (Array<string,1> &v);
  bool ReadInto (Array<string,2> &v);
  bool ReadInto (Array<string,3> &v);
  bool ReadInto (Array<bool,1> &v);
  bool ReadInto (Array<bool,2> &v);
  bool ReadInto (Array<bool,3> &v);

  bool Append (double val);
  bool Append (Array<double,1> &val);
  bool Append (Array<double,2> &val);
  bool Append (int val);
  bool Append (Array<int,1> &val);
  bool Append (Array<int,2> &val);
  bool Append (string val);
  bool Append (Array<string,1> &val);
  bool Append (Array<string,2> &val);

  inline void Print(ofstream &outFile);
  ~VarASCIIClass();
};
 
inline void VarASCIIClass::Print(ofstream &outFile)
{

  if (Dim==0  && 





  if (Dim>0)
    outFile <<"Array<";
  if (Type==DOUBLE_TYPE){
    outFile<<"double";
  }
  else if (Type==INT_TYPE){
    outFile<<"int";
  }
  else if (Type==STRING_TYPE){
    outFile<<"string";
  }
  else if (Type==BOOL_TYPE){
    outFile<<"bool";
  }
  if (Dim>0){
    outFile<<", "; 
    outFile<<Dim<<">";
  }
  outFile<<" "<<Name;
  if (Dim>0){
    outFile<<"(";
    for (int counter=0;counter<Dim;counter++){
      outFile<<
      
    
  


}

/// This class holds an ASCII token, which is just a string and the
/// line number in which it appeared in the file.
class TokenClass
{
public:
  string Str;
  int LineNumber;
};


/// This is the ASCII specialization of IOTreeClass for ASCII text
/// files.  It's syntax is as follows:
/// Section (SectionName)
/// {
///   double x = 3;
///   Array<int,1> y(3) = [1, 2, 3];
///   Array<int,3> z(2,2,1) = [ 1, 2, 
///                             3, 4 ];
///   Section (Species, "species1.h5");
/// }
class IOTreeASCIIClass : public IOTreeClass
{
  /// Reads a text file into a buffer eliminating c++ and c-style
  /// comments.  
  void ReadWithoutComments(string fileName, Array<char,1> &buffer);
  /// Reads a section from a list of TokenClass objects.  iter should
  /// refer to the current place in the list that we should start
  /// reading at.  iter should point to a place just after the '{'.
  /// If wantEndBrace is true, it will look for an ending '}'.
  /// Otherwise it will read until the list of Tokens runs out.  
  bool ReadSection (IOTreeClass *parent, string name,
		    list<TokenClass>::iterator &iter,
		    list<TokenClass> &tokenList,
		    bool wantEndBrace);

 public:
  void FlushFile();
  void WriteSection(ofstream &outFile,int indent);
  /// Print an indented tree of section variable names.
  void PrintTree(int level);
  /// Same thing, just calls above with level 0;
  void PrintTree();
  /// Takes the name of a file to read, the name of my section and a
  /// pointer to my parent.  Reads the file into a tree of
  /// IOTreeClass's.
  bool OpenFile (string filename, string myName, 
		 IOTreeClass *parent);
  /// Do any file handling necessary and delete the whole tree of data.
  void CloseFile();
  IOTreeClass* NewSection(string name);
  void WriteVar(string name, double val);
  void WriteVar(string name, Array<double,1> &val);
  void WriteVar(string name, Array<double,2> &val);
  void WriteVar(string name, Array<double,3> &val);

  void WriteVar(string name, int val);
  void WriteVar(string name, Array<int,1> &val);
  void WriteVar(string name, Array<int,2> &val);
  void WriteVar(string name, Array<int,3> &val);

  void WriteVar(string name, bool val);
  void WriteVar(string name, Array<bool,1> &val);
  void WriteVar(string name, Array<bool,2> &val);
  void WriteVar(string name, Array<bool,3> &val);

  void WriteVar(string name, string val);
  void WriteVar(string name, Array<string,1> &val);
  void WriteVar(string name, Array<string,2> &val);
  void WriteVar(string name, Array<string,3> &val);
  bool NewFile (string fileName,
			string mySectionName,
			IOTreeClass *parent);
  void IncludeSection (IOTreeClass *);


};




#endif
