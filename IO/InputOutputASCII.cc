#include "InputOutput.h"


/// Simply prints 3*num spaces
inline void ASCIIPrintIndent(int num)
{
  for (int counter=0;counter<num*3;counter++)
    cout<<' ';
}


/// Prints an indented hierarchy of sections and variable names to
/// cout. 
void InputTreeASCIIClass::PrintTree(int indentNum)
{
  ASCIIPrintIndent(indentNum);
  cout<<"Section: "<<Name<<endl;
  list<VarClass*>::iterator varIter=VarList.begin();
  while (varIter!=VarList.end()){
    ASCIIPrintIndent(indentNum+1);
    cout<<"Variable: "<<(*varIter)->Name<<" "<<endl;
    varIter++;
  }
  list<InputTreeClass*>::iterator secIter=SectionList.begin();
  while (secIter!=SectionList.end()){
    //    cout<<"Section: "<<(*secIter)->Name<<endl;
    (*secIter)->PrintTree(indentNum+1);
    secIter++;
  }
}

/// Calls PrintTree(0)
void InputTreeASCIIClass::PrintTree()
{
  PrintTree(0);
}




/// Returns true if theChar is a special character that should be
/// parsed into its own token.
bool isSpecial(char theChar)
{
  return ( (theChar=='(') ||
	   (theChar==')') ||
	   (theChar=='{') ||
	   (theChar=='}') ||
	   (theChar=='[') ||
	   (theChar==']') ||
	   (theChar=='<') ||
	   (theChar=='>') ||	
	   (theChar=='=') ||
	   (theChar==';') ||
	   (theChar==','));
}
	   
/// Returns true if theChar is a space, newline, tab, or carriage return.      
bool isWhiteSpace(char theChar)
{
  return ( (theChar=='\n') ||
	   (theChar==' ' ) ||
	   (theChar=='\t') ||
	   (theChar=='\r'));
}
      

/// Returns true if theChar is a letter or underscore		      
bool isAlpha(char theChar)
{
  return ((theChar>='a' && theChar<='z') || (theChar>='A' && theChar<='Z')
	  ||theChar=='_');
}

/// Returns true if theChar is a digit
bool isDigit(char theChar)
{
  return (theChar>='0' && theChar<='9');
}

/// Returns true if theChar is the a valid character for starting a
/// number.  Includes a digit, a '.' or a '-'
bool isNumStart(char theChar)
{
  return ((isDigit(theChar)) || (theChar=='.') || (theChar=='-'));
}

/// Returns true if ch is a valid character comprising a number.
bool isNumChar (char ch)
{
  return (isDigit(ch) || (ch =='.') || (ch=='e') || (ch=='-'));
}


/// Tokenize takes an array of characters and constructs a list of
/// TokenClass objects.  Each token has a string and a line number.
/// Valid tokens are special characters: "(){}[]<>,", quoted strings,
/// words, or numbers.  White space is not significant, except in
/// separating tokens.
void Tokenize(Array<char,1> buffer, list<TokenClass>& tokenList)
{
  int pos=0;
  int lineNum=1;
  while (pos<buffer.size()){
    if (isSpecial(buffer(pos))){
      TokenClass tempToken;
      tempToken.Str+=buffer(pos);
      tempToken.LineNumber = lineNum;
      tokenList.push_back(tempToken);
      pos++;
    }
    else if (buffer(pos)=='\"'){
      TokenClass tempToken;
      tempToken.Str="\"";
      pos++;
      while (buffer(pos)!='\"'){
	tempToken.Str+=buffer(pos);
	pos++;
      }
      pos++;
      tempToken.Str+='\"';
      tempToken.LineNumber = lineNum;
      tokenList.push_back(tempToken);
    }
    else if (isAlpha(buffer(pos))){
      TokenClass tempToken;
      while ((isAlpha(buffer(pos)) || isDigit(buffer(pos)))){
	tempToken.Str+=buffer(pos);
	pos++;
      }
      tempToken.LineNumber = lineNum;
      tokenList.push_back(tempToken);
    }
    else if (isWhiteSpace(buffer(pos))){
      if (buffer(pos)=='\n')
	lineNum++;
      pos++;
    }
    else if (isNumStart(buffer(pos))){
      TokenClass tempToken;
      while (isNumChar(buffer(pos))){
	tempToken.Str+=buffer(pos);
	pos++;
      }
      tempToken.LineNumber = lineNum;
      tokenList.push_back(tempToken);
    }
    else if (buffer(pos)=='\0')
      break;
    else {
      cerr<<"There was a token we do not recognize in line "<<lineNum<<endl;
      cerr <<"The rest of the file is as follows:\n";
      while (pos<buffer.size()) {
	cerr << (int)buffer(pos);
	pos++;
      }
      exit(1);
    }
  }
}	     
	  



/// Just a shortcut to look at two characters at a time.
bool checkPair(Array<char,1> &buffer,int counter,char* toSee)
{
  if (counter+1>=buffer.size()){
    return false;
  }
  if (buffer(counter)==toSee[0] && buffer(counter+1)==toSee[1]){
    return true;
  }
  else return false;

}

/// Reads a file into a character array, removing C and C++ style
/// comments. 
void InputTreeASCIIClass::ReadWithoutComments(string fileName,
					      Array<char,1> 
					      &buffer)
{
  ifstream infile;
  infile.open(fileName.c_str());
  if (!infile.is_open()) {
    cerr << "Cannot open file " << fileName 
	 << " for reading.  Exiting.\n";
    exit(1);
  }
  Array<char,1> tmpBuffer;
  int counter=0;
  bool inQuote=false;
  char dummyChar;
  while (!infile.eof()){
    infile.get(dummyChar);    
    counter++;
  }

  tmpBuffer.resize(counter);
  buffer.resize(counter);
  counter=-1;
  infile.close();
  ifstream infile2;
  infile2.open(fileName.c_str());
  while (!infile2.eof()){
    counter++;
    infile2.get(dummyChar);
    tmpBuffer(counter)=dummyChar;    
  }
  //  cout<<tmpBuffer;
  int bufferLoc=0;
  for (int counter=0;counter<tmpBuffer.size();counter++){
    if (inQuote){
      while ( counter<tmpBuffer.size() && tmpBuffer(counter) != '\"'){
	buffer(bufferLoc)=tmpBuffer(counter);
	counter++;
	bufferLoc++;
      }
      buffer(bufferLoc)=tmpBuffer(counter);
      bufferLoc++;
      inQuote=false;
    }
    else {
      if (checkPair(tmpBuffer,counter,"//")){
	while (tmpBuffer(counter)!='\n' && counter<tmpBuffer.size()){
	  counter++;
	}
	buffer(bufferLoc)=tmpBuffer(counter); //copy the \n over
	bufferLoc++;
      }
      else if (checkPair(tmpBuffer,counter,"/*")){
	while (!checkPair(tmpBuffer,counter,"*/") && counter<tmpBuffer.size()){
	  counter++;
	  if (tmpBuffer(counter)=='\n'){
	    buffer(bufferLoc)=tmpBuffer(counter);
	    bufferLoc++;
	  }		   
	}
	counter++; //end up in the / of comment
      }
      else if (tmpBuffer(counter)=='\"'){
	inQuote=true;
	buffer(bufferLoc)=tmpBuffer(counter);
	bufferLoc++;
      }
      else {
	buffer(bufferLoc)=tmpBuffer(counter);
	bufferLoc++;
      }
    }
  }
  buffer.resizeAndPreserve(bufferLoc);
  infile2.close();
}



/// If isError is true, then we print out an error message giving the
/// line number and the string passed to us in ErrorStr.
inline void ReadAbort (bool isError, int lineNumber, string ErrorStr)
{
  if (isError) {
    cerr << "Error in input file at line number " << lineNumber 
	 << ":\n";
    cerr << ErrorStr;
    exit(1);
  }
}

/// Removes all double quotes from the input string and return it.
string StripQuote(string str)
{
  string newString;
  int i=0;
  assert (str[0] == '\"');
  assert (str[str.length()-1] == '\"');
  while (i<str.length())
    {
      if (str[i] != '\"')
	newString += str[i];
      i++;
    }
  return newString;
}


/// Looks at the string passed to it and returns the corresponding
/// enumerated type.  If the type is not recognized, it returns
/// NOT_ATOMIC.  
AtomicType GetType (string typeString)
{
  if (typeString=="double")
    return DOUBLE_TYPE;
  else if (typeString=="int")
    return INT_TYPE;
  else if (typeString=="string")
    return STRING_TYPE;
  else if (typeString=="bool")
    return BOOL_TYPE;
  else return NOT_ATOMIC;
}


/// Takes a token and reads its value into a double, aborting if there
/// is a problem.
void ReadAtomicVar(TokenClass token,double &d)
{

  char* endPtr;
  d=strtod(token.Str.c_str(),&endPtr);
  ReadAbort(*endPtr!='\0',token.LineNumber,"Expected Double\n");
}

/// Takes a token and reads its value into an int, aborting if there
/// is a problem.
void ReadAtomicVar(TokenClass token,int &d)
{

  char* endPtr;
  d=strtol(token.Str.c_str(),&endPtr,10);
  ReadAbort(*endPtr!='\0',token.LineNumber,"Expected Int\n");
}

/// Takes a token and reads its value into a string, aborting if there
/// is a problem.
void ReadAtomicVar(TokenClass token,string &d)
{

  ReadAbort (token.Str[0] != '\"', token.LineNumber, 
	     "Expected '\"'.");
  ReadAbort (token.Str[token.Str.length()-1] != '\"', token.LineNumber, 
	     "Expected '\"'.");
  d=StripQuote(token.Str);
}

/// Takes a token and reads its value into a bool, aborting if there
/// is a problem.
void ReadAtomicVar(TokenClass token, bool &b)
{
      if (token.Str=="true"){
	b=true;
      }
      else if (token.Str=="false"){
	b=false;
      }
      else ReadAbort(true,token.LineNumber,"Expected true or false\n");
}

/// This template function reads a 1-D array from a token list into
/// the array.  The syntax requires an opening '[' the a
/// comma-separated list of values, then a closing ']'.
template <class T>
void ReadArrayData(list<TokenClass>::iterator &iter,
		   list<TokenClass> &tokenList,
		   Array<T,1> valArray)
{
  ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
  iter++;
  for (int counter=0;counter<valArray.extent(0)-1;counter++){
    ReadAtomicVar(*iter,valArray(counter));
    iter++;
    ReadAbort(iter->Str != ",", iter->LineNumber, "Expected , not found\n");
    iter++;
  }
  //Read last value
  ReadAtomicVar(*iter,valArray(valArray.extent(0)-1));
  iter++;
  ReadAbort(iter->Str != "]", iter->LineNumber, "Expected ] not found\n");
  iter++;
  ReadAbort(iter->Str != ";", iter->LineNumber, "Expected ; not found\n");
  iter++;
}


/// This template function reads a 2-D array from a token list into
/// the array.  The syntax requires an opening '[' the a
/// comma-separated list of values, then a closing ']'.  The data is
/// read row-ordered, i.e. the first index changes fastests as we read
/// in the values.
template <class T>
void ReadArrayData(list<TokenClass>::iterator &iter,
		   list<TokenClass> &tokenList,
		   Array<T,2> valArray)
{
  ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
  iter++;
  for (int i=0;i<valArray.extent(1);i++)
    for (int j=0; j<valArray.extent(0); j++)
      {
	ReadAtomicVar(*iter,valArray(j,i));
	iter++;
	// Read comma if this isn't the last value.
	if ((i!=valArray.extent(1)-1) || 
	    (j!=valArray.extent(0)-1)) {
	  ReadAbort(iter->Str != ",", iter->LineNumber, 
		    "Expected , not found\n");
	  iter++;
	}
      }
  ReadAbort(iter->Str != "]", iter->LineNumber, "Expected ] not found\n");
  iter++;
  ReadAbort(iter->Str != ";", iter->LineNumber, "Expected ; not found\n");
  iter++;
}


/// This template function reads a 3-D array from a token list into
/// the array.  The syntax requires an opening '[' the a
/// comma-separated list of values, then a closing ']'.  The data is
/// read row-ordered, i.e. the first index changes fastests as we read
/// in the values.
template <class T>
void ReadArrayData(list<TokenClass>::iterator &iter,
		   list<TokenClass> &tokenList,
		   Array<T,3> valArray)
{
  ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
  iter++;
  for (int i=0;i<valArray.extent(2);i++)
    for (int j=0; j<valArray.extent(1); j++)
      for (int k=0; k<valArray.extent(0); k++)
      {
	ReadAtomicVar(*iter,valArray(k,j,i));
	iter++;
	// Read comma if this isn't the last value.
	if ((i!=valArray.extent(2)-1) || 
	    (j!=valArray.extent(1)-1) ||
	    (k!=valArray.extent(0)-1)) {
	  ReadAbort(iter->Str != ",", iter->LineNumber, 
		    "Expected , not found\n");
	  iter++;
	}
      }
  ReadAbort(iter->Str != "]", iter->LineNumber, "Expected ] not found\n");
  iter++;
  ReadAbort(iter->Str != ";", iter->LineNumber, "Expected ; not found\n");
  iter++;
}



/// Reads an array from a list of tokens, starting at the token
/// pointed to by iter.  It places the array into the newVar object.
/// It expects to begin reading after the word "Array".
void ReadArray(list<TokenClass>::iterator &iter,
	       list<TokenClass> &tokenList,
	       VarASCIIClass *newVar)
{
  ReadAbort(iter->Str != "<", iter->LineNumber, "Expected < not found\n");
  iter++;
  AtomicType myType=GetType(iter->Str);
  ReadAbort(myType==NOT_ATOMIC,iter->LineNumber,
	    "Array does not have atomic type\n");
  iter++;
  ReadAbort(iter->Str != ",", iter->LineNumber, "Expected , not found\n");
  iter++;
  int numDim;
  ReadAtomicVar(*iter,numDim);
  iter++;
  ReadAbort(iter->Str != ">", iter->LineNumber, "Expected , not found\n");
  iter++;
  
  Array<int,1> dimSize(numDim);
  
  string myName=iter->Str;
  newVar->Name = myName;
  iter++;
  ReadAbort(iter->Str != "(", iter->LineNumber, "Expected ( not found\n");
  iter++;
  for (int counter=0;counter<numDim-1;counter++){
    ReadAtomicVar(*iter,dimSize(counter));
    iter++;
    ReadAbort(iter->Str != ",", iter->LineNumber, "Expected , not found\n");
    iter++;
  }
  //Read the last dimension
  ReadAtomicVar(*iter,dimSize(numDim-1));
  iter++;
  ReadAbort(iter->Str != ")", iter->LineNumber, "Expected ) not found\n");
  iter++;
  ReadAbort(iter->Str!="=",iter->LineNumber,"Expected = not found\n");
  iter++;
  newVar->Dim=numDim;
  newVar->Type=myType;
  if (numDim==1){
    if (myType==INT_TYPE){
      Array<int,1> *valArray=new Array<int,1>(dimSize(0));
      ReadArrayData(iter,tokenList,*valArray);
      newVar->Value=valArray;
    }
    else if (myType==DOUBLE_TYPE){
      Array<double,1> *valArray =new Array<double,1>(dimSize(0));
      ReadArrayData(iter,tokenList,*valArray);
      newVar->Value=valArray;
    }
    else if (myType==BOOL_TYPE){
      Array<bool,1> *valArray=new Array<bool,1>(dimSize(0));
      ReadArrayData(iter,tokenList,*valArray);
      newVar->Value=valArray;
    }
    else if (myType==STRING_TYPE){
      Array<string,1> *valArray=new Array<string,1>(dimSize(0));
      ReadArrayData(iter,tokenList,*valArray);
	newVar->Value=valArray;
    }
  }
  else if (numDim==2){
    if (myType==INT_TYPE){
      Array<int,2> *valArray=new Array<int,2>(dimSize(0),dimSize(1));
      ReadArrayData(iter,tokenList,*valArray);
      newVar->Value=valArray;
    }
    else if (myType==DOUBLE_TYPE){
      Array<double,2> *valArray =new Array<double,2>(dimSize(0),
						     dimSize(1));
      ReadArrayData(iter,tokenList,*valArray);
      newVar->Value=valArray;
    }
    else if (myType==BOOL_TYPE){
      Array<bool,2> *valArray=new Array<bool,2>(dimSize(0),
						dimSize(1));
      ReadArrayData(iter,tokenList,*valArray);
      newVar->Value=valArray;
    }
    else if (myType==STRING_TYPE){
      Array<string,2> *valArray=new Array<string,2>(dimSize(0),
						    dimSize(1));
      ReadArrayData(iter,tokenList,*valArray);
      newVar->Value=valArray;
    }
  }
  else if (numDim==3){
    if (myType==INT_TYPE){
      Array<int,3> *valArray=new Array<int,3>(dimSize(0),
					      dimSize(1),
					      dimSize(2));
      ReadArrayData(iter,tokenList,*valArray);
      newVar->Value=valArray;
    }
    else if (myType==DOUBLE_TYPE){
      Array<double,3> *valArray =new Array<double,3>(dimSize(0),
						     dimSize(1),
						     dimSize(2));
      ReadArrayData(iter,tokenList,*valArray);
      newVar->Value=valArray;
    }
    else if (myType==BOOL_TYPE){
      Array<bool,3> *valArray=new Array<bool,3>(dimSize(0),
						dimSize(1),
						dimSize(2));
      ReadArrayData(iter,tokenList,*valArray);
      newVar->Value=valArray;
      
    }
    else if (myType==STRING_TYPE){
      Array<string,3> *valArray=new Array<string,3>(dimSize(0),
						    dimSize(1),
						    dimSize(2));
      ReadArrayData(iter,tokenList,*valArray);
      newVar->Value=valArray;
    }
  }
  else if (numDim>1){
    cerr<<"We haven't implemented this yet\n";
  }
}



/// This function parses a variable assigment from the list of tokens,
/// creates a new VarASCIIClass object and puts the appropriate value
/// in that object.  It recognizes any of the atomic types or arrays
/// of theose atomic types.
VarASCIIClass* ReadASCIIVar (list<TokenClass>::iterator &iter,
			     list<TokenClass> &tokenList)
{
  VarASCIIClass *newVar = new VarASCIIClass;  
  AtomicType myType=GetType(iter->Str);
  if (myType==NOT_ATOMIC){
    ReadAbort(iter->Str!="Array",iter->LineNumber,
	      "Invalid Type: "+iter->Str+"\n");
    iter++;
    ReadArray(iter,tokenList,newVar);
  }
  else {
    iter++;
    string myName=iter->Str;
    newVar->Name = myName;
    iter++;
    ReadAbort(iter->Str!="=",iter->LineNumber,"Expected equals sign\n");
    iter++;
    TokenClass valToken=*iter;
    iter++;
    ReadAbort(iter->Str!=";",iter->LineNumber,"Expected semicolon\n");
    iter++;
    newVar->Type=myType;
    newVar->Dim=0;
    if (myType==INT_TYPE){
      newVar->Value=new int;
      ReadAtomicVar(valToken,*((int*)newVar->Value));
    }
    else if (myType==DOUBLE_TYPE){
      newVar->Value=new double;
      ReadAtomicVar(valToken,*((double*)newVar->Value));
    }
    else if (myType==STRING_TYPE){
      newVar->Value=new string;
      ReadAtomicVar(valToken,*((string*)newVar->Value));
    }
    else if (myType==BOOL_TYPE){
      newVar->Value=new bool();
      ReadAtomicVar(valToken,*((bool*)newVar->Value));
    }
  }
  return(newVar);
}


/// ReadSection parses a section in the input file.  It takes as
/// arguments this sections parent, its name, a tokenlist iterator,
/// the tokenlist, and a bool telling us if we want to look for a "}"
/// at the end of the input.  If we don't, we keep parsing until the
/// buffer runs out.  Calls itself recursively as necessary, builing
/// up a tree of sections and variables.
bool InputTreeASCIIClass::ReadSection (InputTreeClass *parent,
				       string myName,
				       list<TokenClass>::iterator &iter,
				       list<TokenClass> &tokenList,
				       bool wantEndBrace)
{
  Parent = parent;
  Name = myName;
  while ((iter != tokenList.end()) && (iter->Str != "}")) {
    if (iter->Str == "Section") {
      InputTreeClass *newTree;
      iter++;
      ReadAbort(iter->Str != "(", iter->LineNumber, "Expected ( not found\n");
      iter++;
      string newName = iter->Str;
      iter++;
      // Check for included section
      if (iter->Str == ",") {
	// Get filename
	iter++;
	string fileName = StripQuote(iter->Str);
	iter++;
	ReadAbort (iter->Str!=")", iter->LineNumber, "Expected ) not found\n");
	iter++;
	ReadAbort (iter->Str!=";", iter->LineNumber, "Expected ; not found\n");
	iter++;	
	newTree = ReadTree (fileName, newName, this);
      }
      else {
	ReadAbort(iter->Str != ")", iter->LineNumber, 
		  "Expected ) not found\n");
	iter++;
	ReadAbort(iter->Str != "{", iter->LineNumber, 
		  "Expected { not found\n");
	iter++;
	newTree = new InputTreeASCIIClass();
	((InputTreeASCIIClass*)newTree)->ReadSection((InputTreeClass*)this,
						     newName,iter,tokenList,
						     true);         
      }
      SectionList.push_back(newTree);
    }
    else {
      VarClass *newVar =  ReadASCIIVar(iter, tokenList);
      VarList.push_back(newVar);
    }
  }
  if ((iter==tokenList.end()) && wantEndBrace) {
    cerr << "Unexpected end of file before } \n";
    exit (1);
  }
	    
  if (iter!=tokenList.end())  
    iter++;
  return (true);
}


/// OpenFile takes a filename to open, the name of this section and
/// the parent of this section.  It reads the file into a buffer,
/// converts it to a list of tokens, then parses the tokens,
/// constructing a tree of sections containing variables lists.  
bool InputTreeASCIIClass::OpenFile(string fileName, 
				   string myName, 
				   InputTreeClass *parent)
{
  Array<char,1> buffer;
  ReadWithoutComments(fileName,buffer);
  list<TokenClass> tokenList;
  Tokenize(buffer,tokenList);
  list<TokenClass>::iterator iter=tokenList.begin();
  ReadSection(parent,myName,iter,tokenList, false);
  return true;
}


/// CloseFile recursively destroys the tree of data we constructed.
void InputTreeASCIIClass::CloseFile()
{  
  // First, free all the variables in the list
  while (!VarList.empty()) {
    delete(VarList.front());
    VarList.pop_front();
  }
  
  // Now, call all closes recursively and delete all sections
  while (!SectionList.empty())
    {
      SectionList.front()->CloseFile();
      delete SectionList.front();
      SectionList.pop_front();
    }
}



//////////////////////////////////////////////////////////////////////
//                             ReadInto's                           //
//////////////////////////////////////////////////////////////////////

bool VarASCIIClass::ReadInto (double &val)
{
  assert (Type == DOUBLE_TYPE);
  assert (Dim == 0);
  val = *((double *)Value);
}

bool VarASCIIClass::ReadInto (int &val)
{
  assert (Type == INT_TYPE);
  assert (Dim == 0);
  val = *((int *)Value);
}

bool VarASCIIClass::ReadInto (string &val)
{
  assert (Type == STRING_TYPE);
  assert (Dim == 0);
  val = *((string *)Value);
}

bool VarASCIIClass::ReadInto (bool &val)
{
  assert (Type == BOOL_TYPE);
  assert (Dim == 0);
  val = *((bool *)Value);
}



bool VarASCIIClass::ReadInto (Array<double,1> &val)
{
  assert (Type == DOUBLE_TYPE);
  assert (Dim == 1);
  Array<double,1> &myVal = *((Array<double,1>*)Value);
  val.resize(myVal.extent(0));
  val = myVal;
}
bool VarASCIIClass::ReadInto (Array<double,2> &val)
{
  assert (Type == DOUBLE_TYPE);
  assert (Dim == 2);
  Array<double,2> &myVal = *((Array<double,2>*)Value);
  val.resize(myVal.extent(0), myVal.extent(1));
  val = myVal;
}
bool VarASCIIClass::ReadInto (Array<double,3> &val)
{
  assert (Type == DOUBLE_TYPE);
  assert (Dim == 3);
  Array<double,3> &myVal = *((Array<double,3>*)Value);
  val.resize(myVal.extent(0), myVal.extent(1), myVal.extent(2));
  val = myVal;
}


bool VarASCIIClass::ReadInto (Array<int,1> &val)
{
  assert (Type == INT_TYPE);
  assert (Dim == 1);
  Array<int,1> &myVal = *((Array<int,1>*)Value);
  val.resize(myVal.extent(0));
  val = myVal;
}
bool VarASCIIClass::ReadInto (Array<int,2> &val)
{
  assert (Type == INT_TYPE);
  assert (Dim == 2);
  Array<int,2> &myVal = *((Array<int,2>*)Value);
  val.resize(myVal.extent(0), myVal.extent(1));
  val = myVal;
}
bool VarASCIIClass::ReadInto (Array<int,3> &val)
{
  assert (Type == INT_TYPE);
  assert (Dim == 3);
  Array<int,3> &myVal = *((Array<int,3>*)Value);
  val.resize(myVal.extent(0), myVal.extent(1), myVal.extent(2));
  val = myVal;
}



bool VarASCIIClass::ReadInto (Array<string,1> &val)
{
  assert (Type == STRING_TYPE);
  assert (Dim == 1);
  Array<string,1> &myVal = *((Array<string,1>*)Value);
  val.resize(myVal.extent(0));
  val = myVal;
}
bool VarASCIIClass::ReadInto (Array<string,2> &val)
{
  assert (Type == STRING_TYPE);
  assert (Dim == 2);
  Array<string,2> &myVal = *((Array<string,2>*)Value);
  val.resize(myVal.extent(0), myVal.extent(1));
  val = myVal;
}
bool VarASCIIClass::ReadInto (Array<string,3> &val)
{
  assert (Type == STRING_TYPE);
  assert (Dim == 3);
  Array<string,3> &myVal = *((Array<string,3>*)Value);
  val.resize(myVal.extent(0), myVal.extent(1), myVal.extent(2));
  val = myVal;
}



bool VarASCIIClass::ReadInto (Array<bool,1> &val)
{
  assert (Type == BOOL_TYPE);
  assert (Dim == 1);
  Array<bool,1> &myVal = *((Array<bool,1>*)Value);
  val.resize(myVal.extent(0));
  val = myVal;
}
bool VarASCIIClass::ReadInto (Array<bool,2> &val)
{
  assert (Type == BOOL_TYPE);
  assert (Dim == 2);
  Array<bool,2> &myVal = *((Array<bool,2>*)Value);
  val.resize(myVal.extent(0), myVal.extent(1));
  val = myVal;
}
bool VarASCIIClass::ReadInto (Array<bool,3> &val)
{
  assert (Type == BOOL_TYPE);
  assert (Dim == 3);
  Array<bool,3> &myVal = *((Array<bool,3>*)Value);
  val.resize(myVal.extent(0), myVal.extent(1), myVal.extent(2));
  val = myVal;
}


VarASCIIClass::~VarASCIIClass()
{
  if (Dim==0){
    if (Type==INT_TYPE)
      delete ((int *)Value);   
    else if (Type==DOUBLE_TYPE)
      delete ((double *)Value);
    else if (Type==BOOL_TYPE)
      delete ((bool *)Value);
    else if (Type==STRING_TYPE)
      delete ((string *)Value);
  }
  else if (Dim==1){
    if (Type==INT_TYPE)
      delete ((Array<int,1> *)Value);   
    else if (Type==DOUBLE_TYPE)
      delete ((Array<double,1> *)Value);
    else if (Type==BOOL_TYPE)
      delete ((Array<bool,1> *)Value);
    else if (Type==STRING_TYPE)
      delete ((Array<string,1> *)Value);
  }
  else if (Dim==2){
    if (Type==INT_TYPE)
      delete ((Array<int,2> *)Value);   
    else if (Type==DOUBLE_TYPE)
      delete ((Array<double,2> *)Value);
    else if (Type==BOOL_TYPE)
      delete ((Array<bool,2> *)Value);
    else if (Type==STRING_TYPE)
      delete ((Array<string,2> *)Value);
  }
  else if (Dim==3){
    if (Type==INT_TYPE)
      delete ((Array<int,3> *)Value);   
    else if (Type==DOUBLE_TYPE)
      delete ((Array<double,3> *)Value);
    else if (Type==BOOL_TYPE)
      delete ((Array<bool,3> *)Value);
    else if (Type==STRING_TYPE)
      delete ((Array<string,3> *)Value);
  }
}
