#include "ActionClass.h"
#include <sstream>
#include <iostream>
#include <string>
#include <fstream>

double ActionClass::calcTotalAction(int startSlice, int endSlice, 
				    Array<int,1> changedParticles,
				    int level)
{
  // First, sum the pair actions
  double TotalU = 0.0;
  double TotalK = 0.0;
  int numChangedPtcls = changedParticles.size();
  int species1;
  //  int NumSpecies;
  //  SpeciesArrayClass &IdentPtcls = mySpeciesArray;
  //  NumSpecies = IdentPtcls.Size();
  //  cerr<<"What the action class thinks the size is in calctotalaction: ";
  //  cerr<<  mySpeciesArray->size()<<endl;
  int skip = 1<<level;
  double levelTau = tau* (1<<level);
  for (int ptcl1=0; ptcl1<numChangedPtcls; ptcl1++){
    for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++){
      int species1=Path.ParticleSpeciesNum(ptcl1);
      double notMyself=(double)(ptcl1!=ptcl2);
      for (int slice=startSlice;slice<endSlice;slice+=skip){
	///I think we're double counting particles here
	dVec r1=Path(slice,ptcl1);
	dVec r2=Path(slice,ptcl2);
	dVec rp1=Path(slice+skip,ptcl1);
	dVec rp2=Path(slice+skip,ptcl2);
	dVec r=r1-r2;
	dVec rp=(rp1-rp2);
	double rmag=sqrt(dot(r,r));
	double rpmag=sqrt(dot(rp,rp));
	double s = sqrt(dot (r-rp, r-rp));
	double q = 0.5 * (rmag + rpmag);
	double z = (rmag - rpmag);
	int PairIndex = PairMatrix(species1,
				   Path.ParticleSpeciesNum(ptcl2));
	//I bet it's cheaper to have an if statement that avoids all this nonsense, instead of hte NotMyself thing
	TotalU += notMyself*
	  PairActionVector(PairIndex).calcUsqz(s,q,z, level);
      
      }
    }
    double FourLambdaTauInv=1.0/(4.0*Path.Species(species1).lambda*levelTau);
    for (int slice=startSlice; slice < endSlice;slice+=skip) {
      dVec r1 = Path(slice,ptcl1);
      dVec r2 = Path(slice+skip,ptcl1);
      //This function has to be written and possibly memoized or something?
      double LinkDistSqrd=distSqrd(r1,r2);  
      //We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
      TotalK += LinkDistSqrd*FourLambdaTauInv; 
    }
  }
   


//       Species1 = changedParticles(i)[0];
//       Ptcl1 = changedParticles(i)[1];
//       for (int Species2=0; Species2<NumSpecies; Species2++) {
// 	int NumPtcls2 = IdentPtcls(Species2).NumParticles();
// 	for (int Ptcl2=0; Ptcl2<NumPtcls2; Ptcl2++) {
// 	  for (int Slice=StartSlice; Slice < EndSlice; Slice+=skip) {	    

/////Aren't we calculating everybody in the list 
//// against themselves multiple times here?
// 	    double NotMyself = (double)((Ptcl1!=Ptcl2)||(Species1!=Species2));
// 	    dVec r1 = IdentPtcls(Species1,Ptcl1,Slice);
// 	    dVec r2 = IdentPtcls(Species2,Ptcl2,Slice);
// 	    dVec rp1 = IdentPtcls(Species1,Ptcl1,Slice+skip);
// 	    dVec rp2 = IdentPtcls(Species2,Ptcl2,Slice+skip);
// 	    dVec r = r1 - r2;
// 	    dVec rp = (rp1 -rp2);
// 	    double rmag = sqrt(dot(r,r));
// 	    double rpmag = sqrt(dot(rp,rp));
	    
// 	    double s = sqrt(dot (r-rp, r-rp));
// 	    double q = 0.5 * (rmag + rpmag);
// 	    double z = (rmag - rpmag);
// 	    int PairIndex = PairMatrix(Species1, Species2);
// 	    TotalU += NotMyself*
// 	      PairActionVector(PairIndex).calcUsqz(s,q,z, level);
// 	  }
// 	}

      
      // Now, sum up the kinetic action
//       double FourLambdaTauInv=1.0/(4.0*IdentPtcls(Species1).lambda*levelTau);
//       for (int Slice=StartSlice; Slice < EndSlice; Slice+=skip) {
// 	dVec r1 = IdentPtcls(Species1,Ptcl1,Slice);
// 	dVec r2 = IdentPtcls(Species1,Ptcl1,Slice+skip);
// 	//This function has to be written and possibly memoized or something?
// 	double LinkDistSqrd=distSqrd(r1,r2);  
// 	//We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
// 	TotalK += LinkDistSqrd*FourLambdaTauInv; 
//       }
//     }
    

  

  
  return (TotalK + TotalU);


}


string PairActionClass::SkipTo(ifstream &infile,string skipToString)
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


void PairActionClass::ReadFORTRAN3Tensor(ifstream &infile,Array<double,3> &tempUkj)
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



int GetNextInt(string &s)
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


string GetNextWord(string &s)
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

double GetNextDouble(string &s)
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


void PairActionClass::ReadDavidSquarerFile(string DMFile)
{
  ifstream infile;
  //cout <<DMFile<<endl;
  infile.open(DMFile.c_str());  
  if (infile.fail()){
    cerr<<"CAN'T OPEN THE FILE!!!!!!!!!!";
  }
  
  string numOfFitsString=SkipTo(infile,"SQUARER");
  GetNextWord(numOfFitsString);
  GetNextWord(numOfFitsString);
  GetNextWord(numOfFitsString);
  GetNextWord(numOfFitsString);

  int numOfFits=GetNextInt(numOfFitsString);
  n = numOfFits;
  string NDERIVString = SkipTo(infile,"NDERIV");


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
      for (int levelCounter=0;levelCounter<NumTau;levelCounter++){
	ukj(levelCounter).Init(theGrid,tempUkj(Range::all(),Range::all(),levelCounter));
      }
      tau=smallestTau;
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
      TempukjArray.resize(NumUKJ);      
      dukj.resize(NumTau);
      ReadFORTRAN3Tensor(infile,tempdUkj);
      for (int levelCounter=0;levelCounter<NumTau;levelCounter++){
	dukj(levelCounter).Init(theGrid,tempdUkj(Range::all(),Range::all(),levelCounter));
      }
      tau=smallestTau;
      n=NMax;
      
    }
  }







  /* for (int i=0; i<ukj(6).grid->NumPoints; i++)
    {
      double r = (*ukj(6).grid)(i);
      double u = ukj(6).Params(i,0);
      cerr << r << " " << u  << "\n";
    }


  for (int i=0; i<ukj(6).grid->NumPoints; i++)
    {
      double r = (*dukj(6).grid)(i);
      double du = dukj(6).Params(i,0);
      cerr << r << " " << du  << "\n";
    }

  */

//   for (counter=0;counter<=nderiv;counter++){ // Get the Beta derivatives

//   }

  
  
  

}
