#include "DavidPAClass.h"

#ifdef MAKE_FIT

void DavidPAClass::ReadParams(IOSectionClass &in)
{

}

void DavidPAClass::WriteBetaIndependentInfo (IOSectionClass &out)
{

}

void DavidPAClass::Error (Rho &rho, double &Uerror, double &dUerror)
{

}

void DavidPAClass::AddFit (Rho &rho)
{

}

void DavidPAClass::WriteFits(IOSectionClass &out)
{


}
#endif

double DavidPAClass::U (double q, double z, double s2, int level)
{
  double s=sqrt(s2);
  double uTemp;
  double duTemp;
  double vTemp;
  calcUsqz(s,q,z,level,uTemp,duTemp,vTemp);
  return uTemp;

}

double DavidPAClass::V(double r)
{
  if (r>ukj(0).grid->End){
    r=ukj(0).grid->End;
  }
  if (r<ukj(0).grid->Start){
    r=ukj(0).grid->Start;
  }
  return ukj(0)(0,r);
}

double DavidPAClass::dU(double q, double z, double s2, int level)
{
  double s=sqrt(s2);
  double uTemp;
  double duTemp;
  double vTemp;
  calcUsqz(s,q,z,level,uTemp,duTemp,vTemp);
  //  cerr<<"my vtemp is "<<vTemp<<endl;
  return duTemp;


}
/// Calculate the U(s,q,z) value when given s,q,z and the level 
/*! \f[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n 
  \sum_{j=1}^k u_{kj}(q;\tau)z^{2j}s^{2(k-j)}\f]   */
void DavidPAClass::calcUsqz(double s,double q,double z,int level,
				      double &U, double &dU, double &V)
{

  level=level+(NumTau-(TauPos+1));
  //  cerr<<"My level is "<<level<<endl;
  double rmin = ukj(level).grid->Start;
  
  U=0.0;
  dU=0.0;
  //  level=level+4;
  // Check to make sure we're inside the grid.
  //  //  if (q > ukj(level).grid->End) {
  //  //    U = 0.0; dU=0.0; V = 0.0;
  //  //    return;
  //  //  }
  
  double r=q+0.5*z;
  double rprime=q-0.5*z;
  if (r > ukj(level).grid->End) {
    //    U =0.0 ; dU=0.0; V = 0.0;
    //    return;
    r=ukj(level).grid->End;

  }
  if (rprime > ukj(level).grid->End) {
    //    U = 0.0; dU=0.0; V = 0.0;
    //    return;
    rprime=ukj(level).grid->End;

  }
//   if (q>ukj(level).grid->End){
//     q=ukj(level).grid->End;
//   }
//   q=0.5*(r+rprime);
//   z=z*0.5;
//   s=s*0.5;
 //   if (fabs(z)>q){
//     z=q;
//   }
//   if (s>q){
//     s=q;
//   }
  
  // This is the endpoint action 
  
  
  if ((rprime < rmin) || (r < rmin)){
    U = 5000.0; dU = 0.0;
    return;
  }
  
//   if (q<ukj(level).grid->Start){
//     q=ukj(level).grid->Start;
//   }

  // Compensate for potential, which is subtracted from diaganal action in
  // dm file.
  V = 0.5*(ukj(level)(0,r) + ukj(level)(0,rprime));
  U+= 0.5*(ukj(level)(1,r)+ukj(level)(1,rprime));
  dU+=0.5*((dukj(level))(1,r)+(dukj(level))(1,rprime));
  dU+=V;

//   cerr<<"printing u_10 as function of q<<endl";
//   for (int qInt=0;qInt<200;qInt++){
//     q=(double)qInt/10.0;
//     (ukj(level))(q,TempukjArray); 
//     int k=1;
//     int j=0;
//     cerr<<TempukjArray(k*(k+1)/2+j+1)<< " " <<endl;
//     //  cerr<<dukj(level)(1,q)<<endl;
//   }
//   cerr<<"done printing as a functino of q"<<endl;
//   cerr<<"printing u_11 as function of q<<endl";
//   for (int qInt=0;qInt<200;qInt++){
//     q=(double)qInt/10.0;
//     (dukj(level))(q,TempukjArray); 
//     int k=1;
//     int j=1;
//     cerr<<TempukjArray(k*(k+1)/2+j+1)<<endl;
//   }
//   cerr<<"done printing as a functino of q"<<endl;
//   cerr<<"printing u_20 as function of q<<endl";
//   for (int qInt=0;qInt<200;qInt++){
//     q=(double)qInt/10.0;
//     (dukj(level))(q,TempukjArray); 
//     int k=3;
//     int j=0;
//     cerr<<TempukjArray(k*(k+1)/2+j+1)<<endl;
//   }
//   cerr<<"done printing as a functino of q"<<endl;
//   cerr<<"printing u_21 as function of q<<endl";
//   for (int qInt=0;qInt<200;qInt++){
//     q=(double)qInt/10.0;
//     (dukj(level))(q,TempukjArray); 
//     int k=3;
//     int j=1;
//     cerr<<TempukjArray(k*(k+1)/2+j+1)<<endl;
//   }
//   cerr<<"done printing as a functino of q"<<endl;
//   cerr<<"printing u_22 as function of q<<endl";
//   for (int qInt=0;qInt<200;qInt++){
//     q=(double)qInt/10.0;
//     (dukj(level))(q,TempukjArray); 
//     int k=3;
//     int j=3;
//     cerr<<TempukjArray(k*(k+1)/2+j+1)<<endl;
//   }
//   cerr<<"done printing as a functino of q"<<endl;
//  q=0.5*(r+rprime);

  double UDiag=U;
  //assert(1==2);
  //  return;
  //  if (q>ukj(level).grid->End){
  //    q=ukj(level).grid->End;
  //  }
  if (s > 0.0 && q<ukj(level).grid->End) { // && q<ukj(level).grid->End){
//     if (fabs(z)>2*q){
//       z=2*q*fabs(z)/z;
//     }
//     if (s>2*q){
//       s=2*q;
//     }


    double zsquared=z*z;
    double ssquared=s*s; 
    double ssquaredinverse=1.0/ssquared;
    double Sto2k=ssquared;
    
    (ukj(level))(q,TempukjArray); 
    (dukj(level))(q,TempdukjArray); 
    ////HACK! 
    n=2;
    for (int k=1;k<=n;k++){  
     
      double Zto2j=1;
      double currS=Sto2k;
      
      for (int j=0;j<=k;j++){
	// indexing into the 2darray
	double Ucof  = TempukjArray(k*(k+1)/2+j+1);
	
	double dUcof = TempdukjArray(k*(k+1)/2+j+1);
	U+=(Ucof)*Zto2j*currS;
	dU+=(dUcof)*Zto2j*currS; //+V = HACK!
	Zto2j*=zsquared;
	currS=currS*ssquaredinverse;				
      }				
      Sto2k=Sto2k*ssquared;
    } 
  } 
//   cerr<<U<<" "<<UDiag<<" "<<U-UDiag<<endl;
//   if ((U-UDiag)>1e-1){
//     cerr<<"ERROR! ERROR!"<<endl;
//     cerr<<dU<<" "<<V<<" "<<r<<" "<<rprime<<" "<<s<<" "<<z<<" "<<q<<endl;
//   }
}

void DavidPAClass::ReadDavidSquarerFile(string DMFile)
{


  ifstream infile;
  //cout <<DMFile<<endl;
  infile.open(DMFile.c_str());  
  if (infile.fail()){
    cerr<<"CAN'T OPEN THE FILE!!!!!!!!!!";
  }
  
  string numOfFitsString=SkipTo(infile,"SQUARER");
  cerr<<GetNextWord(numOfFitsString)<<endl;
  //  double LowTau=atof(GetNextWord(numOfFitsString));
  cerr<<GetNextWord(numOfFitsString)<<endl;
  //  cerr<<LowTau<<endl;
  cerr<<GetNextWord(numOfFitsString)<<endl;
  cerr<<GetNextWord(numOfFitsString)<<endl;

  int numOfFits=GetNextInt(numOfFitsString);
  n = numOfFits;
  cerr<<"N is "<<n<<endl;
  // Read in  the potential
  Array<double,1> potential;
  string potGridString = SkipTo(infile, "RANK");
  GetNextWord(potGridString);
  GetNextWord(potGridString);//HACK?
  int numPotPoints = GetNextInt(potGridString);
  potential.resize(numPotPoints);
  Array<double,1> startDeriv(numPotPoints+1); //I think this is the right number of grid poins
  startDeriv=5.0e30;
  Array<double,1> endDeriv(numPotPoints+1);
  endDeriv=0.0;


  SkipTo(infile, "potential");
  cerr<<"I'm reading the potential\n";
  for (int i=0; i<numPotPoints; i++){
    infile >> potential(i);
    cerr<<potential(i)<<endl;
  }
  cerr<<"done\n";
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
      NumTau=GetNextInt(RankString);
      
      
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
      cerr<<"NumTau is"<<NumTau<<endl;
      ukj.resize(NumTau);
      ReadFORTRAN3Tensor(infile,tempUkj);
      Array<double,3> tempUkj2(NumGridPoints,NumUKJ+1,NumTau);
      for(int i=0; i<NumTau; i++){
	tempUkj2(Range::all(),0,i) = potential;
      }
      tempUkj2(Range::all(),Range(1,NumUKJ),Range::all()) = tempUkj;
      tempUkj2(NumGridPoints-1,Range::all(),Range::all())=0.0;
	
      tau=largestTau; //HACK!
      for (int levelCounter=0;levelCounter<NumTau;levelCounter++){//the -3 here is a HACK!
	if (NMax==2){ //MORE HACK!
	  ukj(levelCounter).Init(theGrid,tempUkj2(Range::all(),Range::all(),levelCounter),startDeriv,endDeriv);
	}
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
      cerr<<"NMAX is "<<NMax<<endl;
      if (GetNextInt(beginString)!=2){ //i.e. if it's not U
	cerr<<"ERROR!!! ERROR!!! We didn't get the beta derivative.\n";
      }
      Array<double,3> tempdUkj(NumGridPoints,NumUKJ,NumTau);
      Array<double,3> tempdUkj2(NumGridPoints,NumUKJ+1,NumTau);
      if (NMax==2){
	TempukjArray.resize(NumUKJ+1);      
	TempdukjArray.resize(NumUKJ+1);      
      }
      dukj.resize(NumTau);
      ReadFORTRAN3Tensor(infile,tempdUkj);
      tau=largestTau; //HACK
      for(int i=0; i<NumTau; i++){ //HACK!
	tempdUkj2(Range::all(),0,i) = potential;
	cerr<<"Current tau is "<<tau<<" "<<i<<endl;
	if (fabs(tau-DesiredTau)<1e-10){
	  cerr<<"The tau I've chosen is "<<tau;
	  TauPos=i;
	}
	tau=tau/2; //HACK!
      }
      cerr<<"I'm about ot actually initialize dukj now!"<<endl;
      tempdUkj2(Range::all(),Range(1,NumUKJ),Range::all()) = tempdUkj;
      tempdUkj2(NumGridPoints-1,Range::all(),Range::all())=0.0; ///NOT SURE ABOUT THIS!!!
      for (int levelCounter=0;levelCounter<NumTau;levelCounter++){
	if (NMax==2){ //NMax again
	  dukj(levelCounter).Init(theGrid,tempdUkj2(Range::all(),Range::all(),levelCounter),startDeriv,endDeriv);
	}
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

