#include "GlobalJosephsonMove.h"


double GlobalJosephsonMove::g(int j)
{
  int N=PathData.Path.NumTimeSlices();
  double total=0;
  if (abs(j)==1 || abs(j)==N-2)
    total=total+1.0/(32*Ec*PathData.Path.tau);
  double PiOverN=(M_PI/(double)N);
  double num=PiOverN*PiOverN;
  double denom=sin(PiOverN*j);
  denom=denom*denom;
  total+=PathData.Actions.Josephson.alpha/(8*M_PI*M_PI)*num/denom;
  return total;
}
  
void GlobalJosephsonMove::Read(IOSectionClass &moveInput)
{
  ///do nothing for now

  BuildA();

  
}



void GlobalJosephsonMove::BuildA()
{
  A.resize(PathData.Path.NumTimeSlices()-1);
  A(0)=1.0;
  for (int slice=1;slice<PathData.Path.NumTimeSlices()-1;slice++){
    A(slice)=A(slice-1)*exp(-8*g(slice));

    cerr<<A(slice)<<endl;
  }
  
}
void GlobalJosephsonMove::MakeMove()
{
  PathClass &Path=PathData.Path;
  vector<int> flipList;
  Array<bool,1> flippedAlready(PathData.Path.NumTimeSlices()-1);
  for (int counter=0;counter<flippedAlready.size();counter++)
    flippedAlready(counter)=false;
  //Choose axis
  int n=Path.Random.LocalInt(2*nMax+1)-nMax;
  //Choose j_root
  //  int jroot=Path.Random.LocalInt(Path.NumTimeSlices()-1);
  int jroot=0;
  flipList.push_back(jroot);
  int flipCounter=0;
  flippedAlready(jroot)=true;
  //  cerr<<"Beginning the move: "<<n<<endl;
  double maxDist=abs(Path(0,0)[0]-n*M_PI);
  for (int slice=0;slice<Path.NumTimeSlices()-1;slice++)
    if (maxDist<abs(Path(slice,0)[0]-n*M_PI))
      maxDist=abs(Path(slice,0)[0]-n*M_PI);
  
  while (flipCounter<flipList.size()){
    int flipVal=flipList[flipCounter]; 
    int s=flipVal+1;
    double rMult=1.0; //A(s-flipVal);
    //    cerr<<"Flipval: "<<flipVal<<endl;
    //    cerr<<"Popped: ";
    while (s!=flipVal+A.size()){
      double r=Path.Random.Local();
      //      cerr<<"r: "<<r<<" "<<rMult*r<<" "<<s-flipVal<<" "<<A(s-flipVal)<<endl;
      while ((s!=flipVal+A.size()) && (A(s-flipVal)*exp(-maxDist*maxDist*(s-flipVal))>rMult*r))
	s++;
      if (s==flipVal+A.size())
	cerr<<"we are off the edge!!"<<flipVal<<" "<<endl;
      else{
	double iDiff=8*g(s-flipVal);
	double sDiff=8*g(s-flipVal)*(Path(s % A.size(),0)[0]-n*M_PI)*
	  (Path(flipVal,0)[0]-n*M_PI);
	if (s % A.size() !=flipVal && flippedAlready(s % A.size())==false &&
	    ((1-exp(-sDiff))/(1-exp(-iDiff)*exp(-maxDist*maxDist))>Path.Random.Local())){ //flip!
	  //	  cerr<<"W"<<s%A.size()<<" ";
	  flipList.push_back(s % A.size());
	  flippedAlready(s % A.size())=true;
	  //	  cerr<<s % A.size()<<", ";
	}
	rMult=A(s-flipVal);
	s++;
      }
    }
    //    cerr<<endl;
    flipCounter++;
  }
  for (int slice=0;slice<Path.NumTimeSlices()-1;slice++)
    if (flippedAlready(slice))
      Path(slice,0)[0]=-Path(slice,0)[0]+2*n*M_PI;
  Path(Path.NumTimeSlices()-1,0)=Path(0,0);
  //Recenter
  double phiBar=0;
  for (int slice=0;slice<Path.NumTimeSlices();slice++){
    phiBar+=Path(slice,0)[0];
  }
  phiBar=phiBar/Path.NumTimeSlices();
  int toShift=0;
  while (phiBar>2*M_PI){
    toShift=toShift-1;
    phiBar=phiBar-2*M_PI;
  }
  while (phiBar<-2*M_PI){
    toShift=toShift+1;
    phiBar=phiBar+2*M_PI;
  }
  for (int slice=0;slice<Path.NumTimeSlices();slice++)
    Path(slice,0)[0]=Path(slice,0)[0]+toShift*2*M_PI;
  //Cluster move always gets accepted
  int slice1=0;
  int slice2=PathData.NumTimeSlices()-1;
  ActiveParticles.resize(1);
  ActiveParticles(0)=0;

  PathData.AcceptMove(slice1,slice2,ActiveParticles);
  
  




}
//Move by Werner, Troyer, et al in PRL 95 060201 (2005)
void GlobalJosephsonMove::MakeMoveSlow()
{
  PathClass &Path=PathData.Path;
  vector<int> flipList;
  Array<bool,1> flippedAlready(PathData.Path.NumTimeSlices()-1);
  for (int counter=0;counter<flippedAlready.size();counter++)
    flippedAlready(counter)=false;
  //Choose axis
  int n=Path.Random.LocalInt(2*nMax+1)-nMax;
  //Choose j_root
  int jroot=Path.Random.LocalInt(Path.NumTimeSlices()-1);
  //  double phi_jAxis=Path(jroot,0)[0]-n*M_PI;
  flipList.push_back(jroot);
  int flipCounter=0;
  flippedAlready(jroot)=true;
  while (flipCounter<flipList.size()){
    for (int slice=0;slice<Path.NumTimeSlices()-1;slice++){
      int flipVal=flipList[flipCounter];
      double sDiff=8*g(slice-flipVal)*(Path(slice,0)[0]-n*M_PI)*
	(Path(flipVal,0)[0]-n*M_PI);
      //  if (slice !=jroot && log(1-Path.Random.Local())>-sDiff){ //flip!
      if (slice !=flipVal && flippedAlready(slice)==false &&
	  (1-exp(-sDiff)>Path.Random.Local())){ //flip!
	flipList.push_back(slice);
	flippedAlready(slice)=true;
      }
    }
    flipCounter++;
  }
  for (int slice=0;slice<Path.NumTimeSlices()-1;slice++)
    if (flippedAlready(slice))
      Path(slice,0)[0]=-Path(slice,0)[0]+2*n*M_PI;
  Path(Path.NumTimeSlices()-1,0)=Path(0,0);
  //Recenter
  double phiBar=0;
  for (int slice=0;slice<Path.NumTimeSlices();slice++){
    phiBar+=Path(slice,0)[0];
  }
  phiBar=phiBar/Path.NumTimeSlices();
  int toShift=0;
  while (phiBar>2*M_PI){
    toShift=toShift-1;
    phiBar=phiBar-2*M_PI;
  }
  while (phiBar<-2*M_PI){
    toShift=toShift+1;
    phiBar=phiBar+2*M_PI;
  }
  for (int slice=0;slice<Path.NumTimeSlices();slice++)
    Path(slice,0)[0]=Path(slice,0)[0]+toShift*2*M_PI;
  //Cluster move always gets accepted
  int slice1=0;
  int slice2=PathData.NumTimeSlices()-1;
  ActiveParticles.resize(1);
  ActiveParticles(0)=0;
//   SetMode(NEWMODE);
//   double newA=PathData.Actions.Josephson.SingleAction(slice1,
// 						      slice2,
// 						      ActiveParticles,
// 						      0);
//   SetMode(OLDMODE);
//   double oldA=PathData.Actions.Josephson.SingleAction(slice1,
// 						      slice2,
// 						      ActiveParticles,
// 						      0);
//   SetMode(NEWMODE);
//   cerr<<"Actions: "<<newA<<" "<<oldA<<" "<<newA-oldA<<" "<<n<<" "<<jroot<<endl;

  PathData.AcceptMove(slice1,slice2,ActiveParticles);
}
