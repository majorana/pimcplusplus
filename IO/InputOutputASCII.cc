#include "InputOutput.h"

// void InputTreeASCIIClass::PrintTree()
// {
//   cout<<"Section: "<<sec->Name<<endl;
//    list<VarClass*>::iterator varIter=sec->VarList.begin();
//    while (varIter!=sec->VarList.end()){
//      cout<<"Variable: "<<(*varIter)->Name<<" ";
//      VarASCIIClass *tempVar=(VarASCIIClass*)*varIter;
//      printCharArray(tempVar->Value);
//      //    cout<<tempVar->Value;
//      varIter++;
//    }
//    list<InputTreeClass*>::iterator secIter=sec->SectionList.begin();
//    while (secIter!=sec->SectionList.end()){
//      //    cout<<"Section: "<<(*secIter)->Name<<endl;
//      (*secIter). PrintTree();
//      secIter++;
//    }
// }


inline void ASCIIPrintIndent(int num)
{
  for (int counter=0;counter<num*3;counter++){
    cout<<' ';
  }
}


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

void InputTreeASCIIClass::PrintTree()
{
  PrintTree(0);
}





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
	   
      
bool isWhiteSpace(char theChar)
{
  return ( (theChar=='\n') ||
	   (theChar==' ' ) ||
	   (theChar=='\t') ||
	   (theChar=='\r'));
}
      

		      
bool isAlpha(char theChar)
{
  return ((theChar>='a' && theChar<='z') || (theChar>='A' && theChar<='Z')
	  ||theChar=='_');
}

bool isDigit(char theChar)
{
  return (theChar>='0' && theChar<='9');
}

bool isNumStart(char theChar)
{
  return ((isDigit(theChar)) || (theChar=='.') || (theChar=='-'));
}

bool isNumChar (char ch)
{
  return (isDigit(ch) || (ch =='.') || (ch=='e') || (ch=='-'));
}



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


void 
InputTreeASCIIClass::ReadWithoutComments(string fileName,
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
}


inline void ReadAbort (bool isError, int lineNumber, string ErrorStr)
{
  if (isError) {
    cerr << "Error in input file at line number " << lineNumber 
	 << ":\n";
    cerr << ErrorStr;
    exit(1);
  }
}


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


void ReadAtomicVar(TokenClass token,double &d)
{

  char* endPtr;
  d=strtod(token.Str.c_str(),&endPtr);
  ReadAbort(*endPtr!='\0',token.LineNumber,"Expected Double\n");
}

void ReadAtomicVar(TokenClass token,int &d)
{

  char* endPtr;
  d=strtol(token.Str.c_str(),&endPtr,10);
  ReadAbort(*endPtr!='\0',token.LineNumber,"Expected Int\n");
}

void ReadAtomicVar(TokenClass token,string &d)
{
  d=StripQuote(token.Str);

}


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


template <class T>
void ReadArrayData(list<TokenClass>::iterator &iter,
		   list<TokenClass> &tokenList,
		   Array<T,2> valArray)
{
  ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
  iter++;
  for (int i=0;i<valArray.extent(0);i++)
    for (int j=0; j<valArray.extent(1); j++)
      {
	ReadAtomicVar(*iter,valArray(i,j));
	iter++;
	// Read comma if this isn't the last value.
	if ((i!=valArray.extent(0)-1) || 
	    (j!=valArray.extent(1)-1)) {
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


template <class T>
void ReadArrayData(list<TokenClass>::iterator &iter,
		   list<TokenClass> &tokenList,
		   Array<T,3> valArray)
{
  ReadAbort(iter->Str != "[", iter->LineNumber, "Expected [ not found\n");
  iter++;
  for (int i=0;i<valArray.extent(0);i++)
    for (int j=0; j<valArray.extent(1); j++)
      for (int k=0; k<valArray.extent(2); k++)
      {
	ReadAtomicVar(*iter,valArray(i,j,k));
	iter++;
	// Read comma if this isn't the last value.
	if ((i!=valArray.extent(0)-1) || 
	    (j!=valArray.extent(1)-1) ||
	    (k!=valArray.extent(2)-1)) {
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



bool InputTreeASCIIClass::OpenFile(string fileName, 
				   string myName, 
				   InputTreeClass *parent)
{
  //  Name = myName;
  //  Parent = parent;
  Array<char,1> buffer;
  ReadWithoutComments(fileName,buffer);
  list<TokenClass> tokenList;
  Tokenize(buffer,tokenList);
  list<TokenClass>::iterator iter=tokenList.begin();
  ReadSection(parent,myName,iter,tokenList, false);
//   list<TokenClass>::iterator listIt;
//   listIt=tokenList.begin();
//   while (listIt!=tokenList.end()){
//     cerr<<listIt->Str << " Line Number: " << listIt->LineNumber << endl;
//     listIt++;
//   }
  

  return true;
}



void InputTreeASCIIClass::CloseFile()
{
  return;

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
