#ifndef DavidPAClass_H
#define DavidPAClass_H

#include "PAFitBase.h"
#include "../Blitz.h"
#include <fstream>
#include <iostream>
//#include "../../PathDataClass.h"


///This is the pair action class. It uses the following formula in
///order to calculate the pair action
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
class DavidPAClass : public PairActionFitClass 
{

 private:
  
  ///Holds the Ukj coefficients for a given q
  Array<double,1> TempukjArray;
  Array<double,1> TempdukjArray;

  //  DistanceTableClass *DistanceTable;
  /// Skips to the next string in the file whose substring matches skipToString
  inline string SkipTo(ifstream &infile, string skipToString);
  /// Reads a Fortran 3 tensor
  inline void ReadFORTRAN3Tensor(ifstream &infile, Array<double,3> &tempUkj);
 public:
  Array<double,1> Potential; 
  string type1,type2;
  inline bool Read(IOSectionClass &IOSection,double x, int y);
  inline void Print();
  /// This stores the interaction potential
  CubicSpline V;
  /// This stores the coefficients of the expansion specified above.
  /// The array index is over the levels.  You call it with the q
  /// value and a temporary array to get all of the values in that
  /// column. 
  Array<MultiCubicSpline,1> ukj; ///<(level )
  ///Same as ukj but stores the beta derivatives.
  Array<MultiCubicSpline,1> dukj; ///<(level )
  /// Calculate the U(s,q,z) value when given s,q,z and the level 
  inline void calcUsqz(double s,double q,double z,int level,
		       double &U, double &dU, double &V);
  /// This is the order of the fit to use. 
  int n;
  /// This is the temperature 
  double tau;
  /// Function to read David's squarer file input.
  inline void ReadDavidSquarerFile(string DMFile);
  double U (double q, double z, double s2, int level);
  double dU(double q, double z, double s2, int level);
#ifdef MAKE_FIT
  void ReadParams  (IOSectionClass &inSection);
  void WriteBetaIndependentInfo (IOSectionClass &outSection);
  /// Returns weighter RMS error
  void Error (Rho &rho, double &Uerror, double &dUerror);
  void AddFit (Rho &rho);
  void WriteFits(IOSectionClass &outSection);
#endif
};


/// Calculate the U(s,q,z) value when given s,q,z and the level 
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
inline void DavidPAClass::calcUsqz(double s,double q,double z,int level,
				      double &U, double &dU, double &V)
{
  U=0.0;
  dU=0.0;
  double r=q+0.5*z;
  double rprime=q-0.5*z;
  level=level+4;
  // Check to make sure we're inside the grid.
  if (q > ukj(level).grid->End) {
    U = 0.0; dU=0.0; V = 0.0;
    return;
  }
  if (r > ukj(level).grid->End) {
    U = 0.0; dU=0.0; V = 0.0;
    return;
  }
  if (rprime > ukj(level).grid->End) {
    U = 0.0; dU=0.0; V = 0.0;
    return;
  }

  // This is the endpoint action 
  V = 0.5*(ukj(level)(0,r) + ukj(level)(0,rprime));
  U+=0.5*((ukj(level))(1,r)+(ukj(level))(1,rprime)); 
  dU+=0.5*((dukj(level))(1,r)+(dukj(level))(1,rprime));  

  if (s > 0.0)
    {
      double zsquared=z*z;
      double ssquared=s*s; 
      double ssquaredinverse=1.0/ssquared;
      double Sto2k=ssquared;
      (ukj(level))(q,TempukjArray); 
      (dukj(level))(q,TempdukjArray); 
      for (int k=1;k<=n;k++){  
	
	double Zto2j=1;
	double currS=Sto2k;
	
	for (int j=0;j<=k;j++){
	  // indexing into the 2darray
	  double Ucof  = TempukjArray(k*(k+1)/2+j+1); 
				
	  double dUcof = TempdukjArray(k*(k+1)/2+j+1);
	  U+=(Ucof)*Zto2j*currS;
	  dU+=(dUcof+V)*Zto2j*currS; //+V = HACK!
	  Zto2j*=zsquared;
	  currS=currS*ssquaredinverse;				
	}				
	Sto2k=Sto2k*ssquared;
      }
    }
  //  cerr<<dU<<" "<<V<<" "<<r<<" "<<rprime<<" "<<s<<" "<<z<<" "<<q<<endl;
}

inline bool DavidPAClass::Read(IOSectionClass &in,double x, int y)
{
  string fileName;
  // Read Particles;
  assert(in.OpenSection("Fits"));
  assert(in.OpenSection("Particle1"));
  Particle1.Read(in);
  in.CloseSection();
  assert(in.OpenSection("Particle2"));
  Particle2.Read(in);
  in.CloseSection();
  lambda = Particle1.lambda + Particle2.lambda;
  



  assert(in.ReadVar("Daviddmfile",fileName));
  ReadDavidSquarerFile(fileName.c_str());
  //  assert(in.ReadVar("type1",type1));
  //  assert(in.ReadVar("type2",type2));
  in.CloseSection();
  return true;
}
     

inline string DavidPAClass::SkipTo(ifstream &infile,string skipToString)
{
  string lineString;
  do{
    getline(infile,lineString);

  } while (lineString.find(skipToString,0)==string::npos && !infile.eof());

  if (infile.eof()){
    cerr<<"ERROR!!!  No "<<skipToString<<" found in Davids squarer file\n";
  }
  return lineString;

}


inline void DavidPAClass::ReadFORTRAN3Tensor(ifstream &infile,Array<double,3> &tempUkj)
{

   
  for (int counterTau=0;counterTau<tempUkj.extent(2);counterTau++){
    for (int counterUkj=0;counterUkj<tempUkj.extent(1);counterUkj++){
      for (int counterGridPoints=0;counterGridPoints<tempUkj.extent(0);counterGridPoints++){

	infile>>tempUkj(counterGridPoints,counterUkj,counterTau);
	//	tempUkj(counterGridPoints,counterUkj,counterTau)=0;
      }
    }
  }

}


inline bool IsDigit(char c)
{
  return (((c>='0') && (c<='9')) || (c == '-'));
}

inline bool IsNumberChar(char c)
{
  bool IsNumber = false;
  
  IsNumber = IsNumber || ((c>='0')&&(c<='9'));
  IsNumber = IsNumber || (c=='-');
  IsNumber = IsNumber || (c=='.');
  return IsNumber;
}



inline int GetNextInt(string &s)
{
  int i=0;
  int size = s.size();
  while ((i<size) && (!IsDigit(s[i])))
    i++;
  s = s.substr (i, size-i);

  istringstream sstream(s);

  int num;
  sstream >> num;
  int pos = sstream.tellg();
  if (pos == -1)
    s = "";
  else
    s = s.substr(pos, s.size()-pos);

  return (num);
}


inline string GetNextWord(string &s)
{
  int counter=0;
  while (s[counter] == ' '){
    counter++;
  }
  int startSpot=counter;
  while (s[counter] != ' '){
    counter++;
  }
  string toReturn=s.substr(startSpot,counter-startSpot);
  s=s.substr(counter,s.size()-counter);
  return toReturn;
}

inline double GetNextDouble(string &s)
{
  int i=0;
  int size = s.size();
  while ((i<size) && (!IsNumberChar(s[i])))
    i++;
  s = s.substr (i, size-i);

  istringstream sstream(s);

  double num;
  sstream >> num;
  int pos = sstream.tellg();
  if (pos == -1)
    s = "";
  else
    s = s.substr(pos, s.size()-pos);

  return (num);
}

inline void DavidPAClass::Print()
{ cerr<<"hi\n";}
            
inline void DavidPAClass::ReadDavidSquarerFile(string DMFile)
{
  ifstream infile;
  //cout <<DMFile<<endl;
  infile.open(DMFile.c_str());  
  if (infile.fail()){
    cerr<<"CAN'T OPEN THE FILE!!!!!!!!!!";
  }
  
  string numOfFitsString=SkipTo(infile,"SQUARER");
  cerr<<GetNextWord(numOfFitsString)<<endl;
  cerr<<GetNextWord(numOfFitsString)<<endl;
  cerr<<GetNextWord(numOfFitsString)<<endl;
  cerr<<GetNextWord(numOfFitsString)<<endl;

  int numOfFits=GetNextInt(numOfFitsString);
  n = numOfFits;
  cerr<<endl;
  // Read in  the potential
  Array<double,1> potential;
  string potGridString = SkipTo(infile, "RANK");
  GetNextWord(potGridString);
  GetNextWord(potGridString);//HACK?
  int numPotPoints = GetNextInt(potGridString);
  potential.resize(numPotPoints);
  SkipTo(infile, "potential");
  for (int i=0; i<numPotPoints; i++)
    infile >> potential(i);

  //  string NDERIVString = SkipTo(infile,"NDERIV");


  //  NDERIVString.erase(NDERIVString.find("NDERIV"),strlen("NDERIV"));

  ///  2*(NDERIV+1);
  Grid *theGrid;

  for (int counter=0;counter<=numOfFits;counter++){ //Get the U's 
    string RankString =SkipTo(infile,"RANK");
    int theRank=GetNextInt(RankString);
    //cout<<theRank<<endl;

    if (theRank!=3){
      //cerr<<"ERROR! ERROR! Rank was not 3" << endl;
      counter--;
    }
    else {
      int NumGridPoints=GetNextInt(RankString);
      int NumUKJ=GetNextInt(RankString);
      int NumTau=GetNextInt(RankString);
      
      
      string RGridString =SkipTo(infile,"GRID 1");
      string GridType=GetNextWord(RGridString);
      GridType=GetNextWord(RGridString);
      GridType=GetNextWord(RGridString);
      double startGrid = GetNextDouble(RGridString);
      double endGrid = GetNextDouble(RGridString);
    
      if (GridType=="LINEAR"){
	theGrid=new LinearGrid(startGrid,endGrid,NumGridPoints);
      }
      else if (GridType=="LOG"){
	//cout<<"We're really in log grid here\n";
	double delta=pow((endGrid/startGrid),1.0/(NumGridPoints-1.0));
	//cerr << "delta = " << delta << endl;
	theGrid = new LogGrid(startGrid,delta,NumGridPoints);
      }
      else {
	cerr << "Unrecognized grid type in ReadDavidSquarerFile.\n";
	cerr << "GridType = \"" << GridType << "\"\n";
      }
	  
      
      string TauGridString = SkipTo(infile,"GRID   3"); //We hope this is a log grid
      GetNextWord(TauGridString);
      GetNextWord(TauGridString); /// takes out the Grid  3
      string shouldBeLog;
      if  ((shouldBeLog=GetNextWord(TauGridString))!="LOG"){
	cerr<<"ERROR!!! ERROR!!! The tau grid is not a LOG Grid\n";
	cerr<<shouldBeLog<<endl;
      }
      double smallestTau=GetNextDouble(TauGridString);
      double largestTau=GetNextDouble(TauGridString);
      int numTauCalc=(int)floor(log(largestTau/smallestTau)/log(2.0)+0.5+1.0); ///I think this -1 is correct but who knows
      if (NumTau!=numTauCalc){
	
	cerr<<"ERROR!!! ERROR!!! num tau inconsistency \n";
	cerr<<NumTau<< " "<<numTauCalc<<"  "<<log(largestTau/smallestTau)/log(2.0) + 1.0<< endl;
      }
      string beginString=SkipTo(infile,"BEGIN density matrix table");
      int NMax=GetNextInt(beginString); //This is magically the most accurate fit i.e. NDERIV-1
      if (GetNextInt(beginString)!=1){ //i.e. if it's not U
	cerr<<"ERROR!!! ERROR!!! We got the beta derivative and not U\n";
      }
      Array<double,3> tempUkj(NumGridPoints,NumUKJ,NumTau);
      
      ukj.resize(NumTau);
      ReadFORTRAN3Tensor(infile,tempUkj);
      Array<double,3> tempUkj2(NumGridPoints,NumUKJ+1,NumTau);
      for(int i=0; i<NumTau; i++)
	tempUkj2(Range::all(),0,i) = potential;
      tempUkj2(Range::all(),Range(1,NumUKJ),Range::all()) = tempUkj;
      
      tau=largestTau; //HACK!
      for (int levelCounter=0;levelCounter<NumTau;levelCounter++){//the -3 here is a HACK!
	ukj(levelCounter).Init(theGrid,tempUkj2(Range::all(),Range::all(),levelCounter));
	tau=tau/2; //HACK!
      }
      //      tau=smallestTau; HACK REMOVAL!
      n=NMax;
      
    }
  }



  for (int counter=0;counter<=numOfFits;counter++){ //Get the beta derivative of U's 
    string RankString =SkipTo(infile,"RANK");
    int theRank=GetNextInt(RankString);
    //cout<<theRank<<endl;

    if (theRank!=3){
      //cerr<<"ERROR! ERROR! Rank was not 3" << endl;
      counter--;
    }
    else {
      int NumGridPoints=GetNextInt(RankString);
      int NumUKJ=GetNextInt(RankString);
      int NumTau=GetNextInt(RankString);
      
      
      string RGridString =SkipTo(infile,"GRID 1");
      string GridType=GetNextWord(RGridString);
      GridType=GetNextWord(RGridString);
      GridType=GetNextWord(RGridString);
      double startGrid = GetNextDouble(RGridString);
      double endGrid = GetNextDouble(RGridString);
    
      if (GridType=="LINEAR"){
	theGrid=new LinearGrid(startGrid,endGrid,NumGridPoints);
      }
      else if (GridType=="LOG"){
	//cout<<"We're really in log grid here\n";
	double delta=pow((endGrid/startGrid),1.0/(NumGridPoints-1.0));
	//cerr << "delta = " << delta << endl;
	theGrid = new LogGrid(startGrid,delta,NumGridPoints);
      }
      else {
	cerr << "Unrecognized grid type in ReadDavidSquarerFile.\n";
	cerr << "GridType = \"" << GridType << "\"\n";
      }
	  
      
      string TauGridString = SkipTo(infile,"GRID   3"); //We hope this is a log grid
      GetNextWord(TauGridString);
      GetNextWord(TauGridString); /// takes out the Grid  3
      string shouldBeLog;
      if  ((shouldBeLog=GetNextWord(TauGridString))!="LOG"){
	cerr<<"ERROR!!! ERROR!!! The tau grid is not a LOG Grid\n";
	cerr<<shouldBeLog<<endl;
      }
      double smallestTau=GetNextDouble(TauGridString);
      double largestTau=GetNextDouble(TauGridString);
      int numTauCalc=(int)floor(log(largestTau/smallestTau)/log(2.0)+0.5+1.0); ///I think this -1 is correct but who knows
      if (NumTau!=numTauCalc){
	
	cerr<<"ERROR!!! ERROR!!! num tau inconsistency \n";
	cerr<<NumTau<< " "<<numTauCalc<<"  "<<log(largestTau/smallestTau)/log(2.0) + 1.0<< endl;
      }
      string beginString=SkipTo(infile,"BEGIN density matrix table");
      int NMax=GetNextInt(beginString); //This is magically the most accurate fit i.e. NDERIV-1
      if (GetNextInt(beginString)!=2){ //i.e. if it's not U
	cerr<<"ERROR!!! ERROR!!! We didn't get the beta derivative.\n";
      }
      Array<double,3> tempdUkj(NumGridPoints,NumUKJ,NumTau);
      Array<double,3> tempdUkj2(NumGridPoints,NumUKJ+1,NumTau);
      TempukjArray.resize(NumUKJ+1);      
      TempdukjArray.resize(NumUKJ+1);      
      dukj.resize(NumTau);
      ReadFORTRAN3Tensor(infile,tempdUkj);
      tau=largestTau; //HACK
      for(int i=0; i<NumTau; i++){ //HACK!
	tempdUkj2(Range::all(),0,i) = potential;
	tau=tau/2; //HACK!
      }
      tempdUkj2(Range::all(),Range(1,NumUKJ),Range::all()) = tempdUkj;
      for (int levelCounter=0;levelCounter<NumTau;levelCounter++){
	dukj(levelCounter).Init(theGrid,tempdUkj2(Range::all(),Range::all(),levelCounter));
      }
      //      tau=smallestTau; HACK REMOVAL!
      cerr<<"My tau is "<<tau<<endl;
      n=NMax;
      
    }
  }
  Potential.resize(potential.size());
  for (int counter=0;counter<potential.size();counter++){
    Potential(counter)=potential(counter);
  }
}
#endif
