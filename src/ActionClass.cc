#include "ActionClass.h"

double ActionClass::calcTotalAction(Array<ParticleID,1> changedParticles,
				  int StartSlice, int EndSlice, int level)
{
  // First, sum the pair actions
  double TotalU = 0.0;
  double TotalK = 0.0;
  int NumChangedPtcls = changedParticles.size();
  int Species1, Species2, Ptcl1, Ptcl2;
  int NumSpecies;

  ArrayOfIdenticalParticlesClass &IdentPtcls = *myIdenticalParticleArray;
  NumSpecies = IdentPtcls.size();

  int skip = 1<<level;
  int levelTau = tau*1<<level;
  for (int i=0; i<NumChangedPtcls; i++)
    {
      Species1 = changedParticles(i)[0];
      Ptcl1 = changedParticles(i)[1];
      for (int Species2=0; Species2<NumSpecies; Species2++) {
	int NumPtcls2 = IdentPtcls(Species2).NumParticles;
	for (int Ptcl2=0; Ptcl2<NumPtcls2; Ptcl2++) {
	  for (int Slice=StartSlice; Slice < EndSlice; Slice+=skip) {	    
	    double NotMyself = (double)((Ptcl1!=Ptcl2)||(Species1!=Species2));
	    dVec r1 = IdentPtcls(Species1).Path(Ptcl1,Slice);
	    dVec r2 = IdentPtcls(Species2).Path(Ptcl2,Slice);
	    dVec rp1 = IdentPtcls(Species1).Path(Ptcl1,Slice+skip);
	    dVec rp2 = IdentPtcls(Species2).Path(Ptcl2,Slice+skip);
	    dVec r = r1 - r2;
	    dVec rp = rp1 - rp2;
	    double rmag = sqrt(dot(r,r));
	    double rpmag = sqrt(dot(rp,rp));
	    
	    double s = sqrt(dot (r-rp, r-rp));
	    double q = 0.5 * (rmag + rpmag);
	    double z = (rmag - rpmag);
	    int PairIndex = PairMatrix(Species1, Species2);
	    TotalU += NotMyself*
	      PairActionVector(PairIndex).calcUsqz(s,q,z, level);
	  }
	}
      }
      
      // Now, sum up the kinetic action
      double FourLambdaTauInv=1/(4*IdentPtcls(Species1).lambda*tau);
      for (int Slice=StartSlice; Slice < EndSlice; Slice+=skip) {
	dVec r1 = IdentPtcls(Species1).Path(Ptcl1,Slice);
	dVec r2 = IdentPtcls(Species1).Path(Ptcl1,Slice+skip);
	//This function has to be written and possibly memoized or something?
	double LinkDistSqrd=distSqrd(r1,r2);  
	//We are ignoring the \$\frac{3N}{2}*\log{4*\Pi*\lambda*\tau}
	TotalK += LinkDistSqrd*FourLambdaTauInv; 
      }
    }
    
      return (TotalK+TotalU);
}


string PairActionClass::SkipTo(ifsteam &infile,string skipToString)
{
  string lineString;
  do{
    getline(infile,lineString);
  } while (lineString.find(skipToString,0)==string::npos && !infile.eof());

  if (infile.eof()){
    cout<<"ERROR!!! ERROR!!! ERROR! No NDERIV found in Davids squarer file\n";
  }
  return lineString;


}


void PairActionClass::ReadFORTRAN3Tensor(ifstream &infile,Array<double,3> &tempUkj)
{

  
  for (int counterTau=0;counterTau<tempUkj.size(3);counterTau++){
    for (int counterUkj=0;counterUkj<tempUkj.size(2);counterUkj++){
      for (int counterGridPoints;counterGridPoints<tempUkj.size(1);counterGridPoints++){
	infile>>tempUkj(counterGridPoints,counterUkj,counterTau);
      }
    }
  }

}

void PairActionClass::ReadDavidSquarerFile(string DMFile)
{
  ifstream infile;
  infile.open(DMFile);  
  

  string NDERIVString = SkipTo(infile,"NDERIV");
  int numOfFits=GetNextInt(NDERIVString);
  //  NDERIVString.erase(NDERIVString.find("NDERIV"),strlen("NDERIV"));

  ///  2*(NDERIV+1);
  Grid *theGrid;
  for (counter=0;counter<=numOfFits;counter++){ //Get the U's 
    string RankString =SkipTo("RANK");
    int theRank=GetNextInt(RankString);
    if (theRank!=3){
      cerr<<"ERROR! ERROR! Rank was not 3";
    }
    int numGridPoints=GetNextInt(RankString);
    int numUKJ=GetNextInt(RankString);
    int numTau=GetNextInt(RankString);
    
    string RGridString =SkipTo("GRID 1");
    string GridType=GetNextWord(RGridString);
    double startGrid = GetNextDouble(RGridString);
    double endGrid = GetNextDouble(RGridString);
    
    if (GridType=="LINEAR"){
      theGrid=new LinearGrid(startGrid,endGrid,numGridPoints);
    }
    else if (GridType=="LOG"){
      double delta=pow((endGrid/startGrid),1.0/(numGridPoints-1.0));
      theGrid = new LogGrid(startGrid,delta,numGridPoints);
    }
    
    String TauGridString = SkipTo("GRID  3"); //We hope this is a log grid
    if  (getNextWord(TauGridString)!="LOG"){
      cerr<<"ERROR!!! ERROR!!! The tau grid is not a LOG Grid\n";
    }
    double smallestTau=getNextDouble(TauGridString);
    double largestTau=getNextDouble(TauGridString);
    int numTauCalc=floor(log(largestTau/smallestTau)/log(2)+0.5);
    if (numTau!=numTauCalc){
      cerr<<"ERROR!!! ERROR!!! num tau inconsistency \n";
    }
    String beginString=SkipTo("BEGIN density matrix table");
    NMax=GetNextInt(beginString); //This is magically the most accurate fit i.e. NDERIV-1
    if (GetNextInt(beginString)!=1){ //i.e. if it's not U
      cerr<<"ERROR!!! ERROR!!! We got the beta derivative and not U";
    }
    Array<double,3> tempUkj(NumGridPoints,NumUkj,NumTau);

    ReadFORTRAN3Tensor(infile,tempUkj);
    for (levelCounter=0;levelCounter<numTau;levelCounter++){
      ukj(levelCounter).Init(theGrid,tempUkj(Range::all(),Range::all(),levelCounter));
    }
    tau=smallestTau;
    n=NMax;
	
    
    
  }

  for (counter=0;counter<=nderiv;counter++){ // Get the Beta derivatives

  }

  
  
  
  
  




}
