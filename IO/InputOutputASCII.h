#ifndef INPUT_OUTPUT_ASCII_H
#define INPUT_OUTPUT_ASCII_H
#include "InputOutputBase.h"
#include <iostream>
#include <stack>
#include <string>
#include <list>



class VarASCIIClass : public VarClass
{  
protected:
  inline void ComplainReadInto()
  { 
    cerr << "Trying to read into wrong type in ReadInto.\n";
    exit(1);
  }
  inline void ComplainAppend()
  {
    cerr << "Trying to Append to a variable of wrong type.\n";
  }

public:
  /// Hack to hold different types.  We cast to to whatever type is
  /// indicated by the Type and Dim variables in the parent class.

  /// The ReadIntos take copy their contents into val if the type
  /// that's passed to them matches the type that's actually stored in
  /// it's void *Value, as indicated by Type and Dim.
  virtual bool ReadInto (double &val);
  virtual bool ReadInto (int &val);
  virtual bool ReadInto (string &val);
  virtual bool ReadInto (bool &val);
  virtual bool ReadInto (Array<double,1> &v);
  virtual bool ReadInto (Array<double,2> &v);
  virtual bool ReadInto (Array<double,3> &v);
  virtual bool ReadInto (Array<int,1> &v);
  virtual bool ReadInto (Array<int,2> &v);
  virtual bool ReadInto (Array<int,3> &v);
  virtual bool ReadInto (Array<string,1> &v);
  virtual bool ReadInto (Array<string,2> &v);
  virtual bool ReadInto (Array<string,3> &v);
  virtual bool ReadInto (Array<bool,1> &v);
  virtual bool ReadInto (Array<bool,2> &v);
  virtual bool ReadInto (Array<bool,3> &v);

  virtual bool Append (double val);
  virtual bool Append (Array<double,1> &val);
  virtual bool Append (Array<double,2> &val);
  virtual bool Append (int val);
  virtual bool Append (Array<int,1> &val);
  virtual bool Append (Array<int,2> &val);
  virtual bool Append (string val);
  virtual bool Append (Array<string,1> &val);
  virtual bool Append (Array<string,2> &val);
  virtual bool Append (bool val);
  virtual bool Append (Array<bool,1> &val);
  virtual bool Append (Array<bool,2> &val);


  virtual void Print(ofstream &outFile) = 0;
};

class VarASCIIdouble0Class : public VarASCIIClass
{
public:
  double Value;
  bool ReadInto (double &val);
  void Print(ofstream &outFile);
};

class VarASCIIdouble1Class : public VarASCIIClass
{
public:
  Array<double,1> Value;
  bool ReadInto (Array<double,1> &val);
  bool Append (double val);
  void Print(ofstream &outFile);
};

class VarASCIIdouble2Class : public VarASCIIClass
{
public:
  Array<double,2> Value;
  bool ReadInto (Array<double,2> &val);
  bool Append (Array<double,1> &val);
  void Print(ofstream &outFile);
};

class VarASCIIdouble3Class : public VarASCIIClass
{
public:
  Array<double,3> Value;
  bool ReadInto (Array<double,3> &val);
  bool Append (Array<double,2> &val);
  void Print(ofstream &outFile);
};


class VarASCIIint0Class : public VarASCIIClass
{
public:
  int Value;
  bool ReadInto (int &val);
  void Print(ofstream &outFile);
};

class VarASCIIint1Class : public VarASCIIClass
{
public:
  Array<int,1> Value;
  bool ReadInto (Array<int,1> &val);
  bool Append (int val);
  void Print(ofstream &outFile);
};

class VarASCIIint2Class : public VarASCIIClass
{
public:
  Array<int,2> Value;
  bool ReadInto (Array<int,2> &val);
  bool Append (Array<int,1> &val);
  void Print(ofstream &outFile);
};

class VarASCIIint3Class : public VarASCIIClass
{
public:
  Array<int,3> Value;
  bool ReadInto (Array<int,3> &val);
  bool Append (Array<int,2> &val);
  void Print(ofstream &outFile);
};

 
class VarASCIIstring0Class : public VarASCIIClass
{
public:
  string Value;
  bool ReadInto (string &val);
  void Print(ofstream &outFile);
};

class VarASCIIstring1Class : public VarASCIIClass
{
public:
  Array<string,1> Value;
  bool ReadInto (Array<string,1> &val);
  bool Append (string val);
  void Print(ofstream &outFile);
};

class VarASCIIstring2Class : public VarASCIIClass
{
public:
  Array<string,2> Value;
  bool ReadInto (Array<string,2> &val);
  bool Append (Array<string,1> &val);
  void Print(ofstream &outFile);
};

class VarASCIIstring3Class : public VarASCIIClass
{
public:
  Array<string,3> Value;
  bool ReadInto (Array<string,3> &val);
  bool Append (Array<string,2> &val);
  void Print(ofstream &outFile);
};

 

class VarASCIIbool0Class : public VarASCIIClass
{
public:
  bool Value;
  bool ReadInto (bool &val);
  void Print(ofstream &outFile);
};

class VarASCIIbool1Class : public VarASCIIClass
{
public:
  Array<bool,1> Value;
  bool ReadInto (Array<bool,1> &val);
  bool Append (bool val);
  void Print(ofstream &outFile);
};

class VarASCIIbool2Class : public VarASCIIClass
{
public:
  Array<bool,2> Value;
  bool ReadInto (Array<bool,2> &val);
  bool Append (Array<bool,1> &val);
  void Print(ofstream &outFile);
};

class VarASCIIbool3Class : public VarASCIIClass
{
public:
  Array<bool,3> Value;
  bool ReadInto (Array<bool,3> &val);
  bool Append (Array<bool,2> &val);
  void Print(ofstream &outFile);
};

 


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
  void WriteSection(ofstream &outFile,int indent);
  /// Print an indented tree of section variable names.
  void PrintTree(int level);
  /// Same thing, just calls above with level 0;
  void PrintTree();

  IOTreeClass* NewSection(string name);
  void IncludeSection (IOTreeClass *);
  /// Takes the name of a file to read, the name of my section and a
  /// pointer to my parent.  Reads the file into a tree of
  /// IOTreeClass's.
  bool OpenFile (string filename, string myName, 
		 IOTreeClass *parent);
  bool NewFile (string fileName, string mySectionName,
		IOTreeClass *parent);
  /// Do any file handling necessary and delete the whole tree of data.
  void CloseFile();
  void FlushFile();

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
  IOTreeASCIIClass()
  { IsModified = false; }
};




#endif
