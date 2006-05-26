

// #include <gsl/gsl_linalg.h>
// #include <gsl/gsl_blas.h>
#include <algorithm>
#include "TruncatedInverse.h"
#include "../PathDataClass.h"
#include <Common/MatrixOps/MatrixOps.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "fstream.h"
// #ifdef ORDER_N_FERMIONS
extern "C"{
#include "../det_calc_uekt.h"
}
/// #endif
// int testme();


//We really want to know what the changed particle is here. Currently
//I'm going to try to just store it when we are doing the singleactionbelow
void 
TruncatedInverseClass::AcceptCopy(int slice1, int slice2)
{
//   cout<<"I'm being accepted! YEA!"<<ChangedColumn<<endl;
//   SetMode(OLDMODE);
//   cerr<<"Pre-checking"<<endl;
//   CheckDeterminantMatrix();
//   SetMode(NEWMODE);
//   cerr<<"post checking"<<endl;
//   for (int theRow=0;theRow<DetMatrix.extent(0);theRow++){
//     DetMatrix(ChangedColumn,theRow)=newCol(theRow);
//   }
//   for (int col=0;col<DetMatrix.extent(0);col++){
//     if (col!=ChangedColumn)
//       DetMatrix(col,ChangedColumn)=newCol(col);
//   }
  
  //  CheckDeterminantMatrix();
  //  cerr<<DetMatrix<<endl;
  //  cerr<<"My new determinant is "<<Determinant(DetMatrix)<<endl;
  //  BuildDeterminantMatrix();
}

struct doubleint
{
  double distance;
  int ptcl;
};
struct intdouble
{
  double distance;
  int ptcl;
};

bool operator<(const doubleint& a, const doubleint& b) {
  return a.distance < b.distance;
}

bool operator<(const intdouble& a, const intdouble& b) {
  return a.ptcl < b.ptcl;
}


void 
TruncatedInverseClass::RejectCopy(int slice1, int slice2)
{
  //  cerr<<"I'm being rejected! YEA!"<<endl;
  
}


void 
TruncatedInverseClass::BuildSmallDeterminantMatrix()
{
  cerr<<"Into it"<<endl;
  double T=1.0/(Path.TotalNumSlices*Path.tau);
  double lambda=PathData.Path.Species(0).lambda;
  Array<bool,1> doneAlready(PathData.Path.NumParticles());
  doneAlready=false;
  vector<doubleint> determinantPtcl;
  vector<intdouble> tempdeterminantPtcl;
  int xBox,yBox,zBox;
  int slice=0;
  int ptcl1=ChangedColumn;
  intdouble toPushBack;
  toPushBack.distance=0.0;
  toPushBack.ptcl=ptcl1;
  tempdeterminantPtcl.push_back(toPushBack);
  doneAlready(ptcl1)=true;
  Path.Cell.FindBox(olddvec,xBox,yBox,zBox);
  //      cerr<<"Beginning"<<endl;
  for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
    int rxbox,rybox,rzbox;
    rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
    rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
    rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
    list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
    for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
      int ptcl2=*i;
      //      if (!doneAlready(ptcl2))

      toPushBack.ptcl=ptcl2;
      dVec disp=PathData.Path(0,ptcl2)-olddvec;
      PathData.Path.PutInBox(disp);
      toPushBack.distance=sqrt(dot(disp,disp));
      tempdeterminantPtcl.push_back(toPushBack);
      doneAlready(ptcl2)=true;
    }
  }



  Path.Cell.FindBox(newdvec,xBox,yBox,zBox);
  //      cerr<<"Beginning"<<endl;
  for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
    int rxbox,rybox,rzbox;
    rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
    rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
    rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
    list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
    for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
      int ptcl2=*i;
      //      if (!doneAlready(ptcl2))
      //	tempdeterminantPtcl.push_back(ptcl2);
      toPushBack.ptcl=ptcl2;
      dVec disp=PathData.Path(0,ptcl2)-newdvec;
      PathData.Path.PutInBox(disp);
      toPushBack.distance=sqrt(dot(disp,disp));
      tempdeterminantPtcl.push_back(toPushBack);
      doneAlready(ptcl2)=true;
    }
  }

  sort(tempdeterminantPtcl.begin(),tempdeterminantPtcl.end());

  for (int counter=0;counter<tempdeterminantPtcl.size();counter++){
    doubleint temp;
    temp.ptcl=tempdeterminantPtcl[counter].ptcl;
    temp.distance=tempdeterminantPtcl[counter].distance;
    determinantPtcl.push_back(temp);
    if (counter+1!=tempdeterminantPtcl.size()-1 && 
	counter!=tempdeterminantPtcl.size() &&
	tempdeterminantPtcl[counter].ptcl==tempdeterminantPtcl[counter+1].ptcl &&
	tempdeterminantPtcl[counter].ptcl==tempdeterminantPtcl[counter+2].ptcl)
      counter=counter+2;
    else if (counter!=tempdeterminantPtcl.size()-1 && 
	tempdeterminantPtcl[counter].ptcl==tempdeterminantPtcl[counter+1].ptcl)
      counter=counter+1;
  }
  sort(determinantPtcl.begin(),determinantPtcl.end());
  for (int counter=0;counter<determinantPtcl.size();counter++)
    cerr<<"Distance: "<<counter<<" "<<determinantPtcl[counter].distance<<endl;
  //  for (int counter=0;counter<determinantPtcl.size();counter++)
  //    cerr<<determinantPtcl[counter].ptcl<<" "<<determinantPtcl[counter].distance<<endl;


//   int N=determinantPtcl.size();
//   cerr<<"Size is "<<N<<endl;
//   SmallDetMatrix.resize(N,N);
//   int counter1=-1;

//   vector<doubleint>::const_iterator iter1;
//   vector<doubleint>::const_iterator iter2;
//   for (iter1=determinantPtcl.begin();iter1!=determinantPtcl.end();iter1++){
//     int ptcl1=(*iter1).ptcl;
//     counter1++;
//     int counter2=-1;
//     for (iter2=determinantPtcl.begin();iter2!=determinantPtcl.end();iter2++){
//       counter2++;
//       int ptcl2=(*iter2).ptcl;
//       dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
//       PathData.Path.PutInBox(disp);
//       double dist2=dot(disp,disp);
//       //      cerr<<ptcl1<<" "<<ptcl2<<endl;
//       SmallDetMatrix(counter1,counter2)=exp(-T*dist2/(4.0*lambda));
//     }
//   }      



//  int N=determinantPtcl.size();
//  cerr<<"Size is "<<N<<endl;
  
  for (int N=1;N<=determinantPtcl.size();N++){
    SmallDetMatrix.resize(N,N);
    for (int counter1=0;counter1<N;counter1++){
      int ptcl1=determinantPtcl[counter1].ptcl;
      for (int counter2=0;counter2<N;counter2++){
	int ptcl2=determinantPtcl[counter2].ptcl;
	dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
	PathData.Path.PutInBox(disp);
	double dist2=dot(disp,disp);
	SmallDetMatrix(counter1,counter2)=exp(-T*dist2/(4.0*lambda));
      }
    }

    DeterminantList(N-1)=Determinant(SmallDetMatrix);
  }
  
  //  cerr<<SmallDetMatrix<<endl;



  //  cerr<<SmallDetMatrix<<endl;
  cerr<<"out of it"<<endl;
}


// void 
// TruncatedInverseClass::BuildSmallDeterminantMatrix()
// {
//   cerr<<"Into it"<<endl;
//   double T=1.0/(Path.TotalNumSlices*Path.tau);
//   double lambda=PathData.Path.Species(0).lambda;

//   list<int> determinantPtcl;
//   int xBox,yBox,zBox;
//   int slice=0;
//   int ptcl1=ChangedColumn;
//   //  determinantPtcl.push_back(ptcl1);
//   Path.Cell.FindBox(Path(slice,ptcl1),xBox,yBox,zBox);
//   //      cerr<<"Beginning"<<endl;
//   for (int cellVal=0;cellVal<Path.Cell.AffectedCells.size();cellVal++){
//     int rxbox,rybox,rzbox;
//     rxbox=(xBox+Path.Cell.AffectedCells(cellVal)[0] +2 * Path.Cell.GridsArray.extent(0)) % Path.Cell.GridsArray.extent(0);
//     rybox=(yBox+Path.Cell.AffectedCells(cellVal)[1] + 2 * Path.Cell.GridsArray.extent(1)) % Path.Cell.GridsArray.extent(1);
//     rzbox=(zBox+Path.Cell.AffectedCells(cellVal) [2] +2 * Path.Cell.GridsArray.extent(2)) % Path.Cell.GridsArray.extent(2);
//     list<int> &ptclList=Path.Cell.GridsArray(rxbox,rybox,rzbox).Particles(slice);
//     for (list<int>::iterator i=ptclList.begin();i!=ptclList.end();i++) {
//       int ptcl2=*i;
//       determinantPtcl.push_back(ptcl2);
//     }
//   }
//   int N=determinantPtcl.size();
//   cerr<<"Size is "<<N<<endl;
//   SmallDetMatrix.resize(N,N);
//   int counter1=-1;

//   list<int>::const_iterator iter1;
//   list<int>::const_iterator iter2;
//   for (iter1=determinantPtcl.begin();iter1!=determinantPtcl.end();iter1++){
//     int ptcl1=*iter1;
//     counter1++;
//     int counter2=-1;
//     for (iter2=determinantPtcl.begin();iter2!=determinantPtcl.end();iter2++){
//       counter2++;
//       int ptcl2=*iter2;
//       //      dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
//       //      PathData.Path.PutInBox(disp);
//       //      double dist2=dot(disp,disp);
//       //      cerr<<ptcl1<<" "<<ptcl2<<endl;
//       SmallDetMatrix(counter1,counter2)=DetMatrix(ptcl1,ptcl2); // exp(-T*dist2/(4.0*lambda));
//     }
//   }      
//   //  cerr<<SmallDetMatrix<<endl;
//   cerr<<"out of it"<<endl;
// }

void 
TruncatedInverseClass::BuildDeterminantMatrix()
{
  cerr<<"Original buildling"<<endl;
  double T=1.0/(Path.TotalNumSlices*Path.tau);
  double lambda=PathData.Path.Species(0).lambda;
  for (int ptcl1=0;ptcl1<Path.NumParticles();ptcl1++){
    for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++){
      dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
      PathData.Path.PutInBox(disp);
      double dist2=dot(disp,disp);
      DetMatrix(ptcl1,ptcl2)=exp(-T*dist2/(4.0*lambda));
    }
  }
  //  CheckDeterminantMatrix();
}


void 
TruncatedInverseClass::CheckDeterminantMatrix()
{

  double T=1.0/(Path.TotalNumSlices*Path.tau);
  double lambda=PathData.Path.Species(0).lambda;
  for (int ptcl1=0;ptcl1<Path.NumParticles();ptcl1++){
    for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++){
      dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
      PathData.Path.PutInBox(disp);
      double dist2=dot(disp,disp);
      if (DetMatrix(ptcl1,ptcl2)-exp(-T*dist2/(4.0*lambda))>=1e-3){
	cerr<<ptcl1<<" "<<ptcl2<<endl;
	cerr<<DetMatrix(ptcl1,ptcl2)<<endl;
	cerr<<exp(-T*dist2/(4.0*lambda))<<endl;
      }
      assert(DetMatrix(ptcl1,ptcl2)-exp(-T*dist2/(4.0*lambda))<=1e-3);
    }
  }
}




TruncatedInverseClass::TruncatedInverseClass (PathDataClass &pathData) :
					    
  NodalActionClass (pathData), 
  Path (pathData.Path)

{
}

void 
TruncatedInverseClass::Read (IOSectionClass &in)
{
  // Do nothing for now
  //  SetupFreeActions();
  cerr<<"Initializing"<<endl;
  int N=Path.NumParticles();
  DetMatrix.resize(N,N);
  SmallDetMatrix.resize(N,N);
  DeterminantList.resize(N);

}


double 
TruncatedInverseClass::SingleAction (int startSlice, int endSlice,
				  const Array<int,1> &changePtcls,
				  int level)
{
#ifdef ORDER_N_FERMIONS
  ChangedColumn=changePtcls(0);
  SetMode(OLDMODE);
  olddvec=Path(0,ChangedColumn);
  SetMode(NEWMODE);
  newdvec=Path(0,ChangedColumn);
  SetMode(OLDMODE);
  //  BuildDeterminantMatrix();
  BuildSmallDeterminantMatrix();
  cerr<<"Matrix size is "<<SmallDetMatrix.extent(0)<<endl;
  ofstream infile;
//   infile.open("bigMatrixA");
//   for (int i=0;i<DetMatrix.extent(0);i++){
//     for (int j=0;j<DetMatrix.extent(1);j++){
//       infile<<DetMatrix(i,j);
//       infile<<" ";
//     }
//     infile<<endl;
//   }
//   infile.close();

//   infile.open("smallMatrixA");
//   for (int i=0;i<SmallDetMatrix.extent(0);i++){
//     for (int j=0;j<SmallDetMatrix.extent(1);j++){
//       infile<<SmallDetMatrix(i,j);
//       infile<<" ";
//     }
//     infile<<endl;
//   }
//   infile.close();

  //  cerr<<"Det matrix is "<<DetMatrix<<endl;
  //  cerr<<"Small det matrix is "<<SmallDetMatrix<<endl;
  //  //  //  //  Array<double,2> smallInverse=Inverse(SmallDetMatrix);
  //  double det_old=Determinant(DetMatrix);
  double det_old=1.0;
  double small_det_old=Determinant(SmallDetMatrix);
  int oldSize=SmallDetMatrix.extent(0);
  SetMode(NEWMODE);
  //  BuildDeterminantMatrix();
  BuildSmallDeterminantMatrix();
  cerr<<"Matrix size is "<<SmallDetMatrix.extent(0)<<endl;
  //  //  //  //  Array<double,2> MultToDet;
  //  //  //  //  MultToDet.resize(SmallDetMatrix.extent(0),SmallDetMatrix.extent(1));
  //  //  //  //  MatMult(smallInverse,SmallDetMatrix,MultToDet);
  //  MultToDet*=10;
  //  cerr<<smallInverse<<endl;

//   infile.open("bigMatrixA2");
//   for (int i=0;i<DetMatrix.extent(0);i++){
//     for (int j=0;j<DetMatrix.extent(1);j++){
//       infile<<DetMatrix(i,j);
//       infile<<" ";
//     }
//     infile<<endl;
//   }
//   infile.close();

//   infile.open("smallMatrixA2");
//   for (int i=0;i<SmallDetMatrix.extent(0);i++){
//     for (int j=0;j<SmallDetMatrix.extent(1);j++){
//       infile<<SmallDetMatrix(i,j);
//       infile<<" ";
//     }
//     infile<<endl;
//   }
//   infile.close();

  //  sleep(100);
  //  cerr<<"Determinant: "<<Determinant(MultToDet)<<endl;
  //  double det_new=Determinant(DetMatrix);
  double det_new=1.0;
  double small_det_new=Determinant(SmallDetMatrix);
  int newSize=SmallDetMatrix.extent(0);
  cerr<<det_old<<" "<<det_new<<" "<<small_det_old<<" "<<small_det_new<<endl;
  cerr<<det_new/det_old<<" "<<small_det_new/small_det_old<<endl;
  cerr<<oldSize<<" "<<newSize<<endl;
  
  for (int i=0;i<newSize;i++){
    SetMode(OLDMODE);
    small_det_old=DeterminantList(i);
    SetMode(NEWMODE);
    small_det_new=DeterminantList(i);
    cerr<<"Using "<<i<<" "<<small_det_new/small_det_old<<endl;
  }

  return 1.0;

#else
  cerr << "PIMC++ was not configured with --enable-on-fermions.\n"
       << "This function doesn't work.  Please reconfigure.\n";
  abort();
  return 0.0;
#endif

}


//just here sot hatt it will compile
bool
TruncatedInverseClass::IsPositive(int x)
{
  return true;
}

double 
TruncatedInverseClass::d_dBeta (int slice1, int slice2, int level)
{ 
  return 0.0;
}

NodeType TruncatedInverseClass::Type()
{
  return FREE_PARTICLE;
}


bool
TruncatedInverseClass::IsGroundState()
{
  return (false);
}

void
TruncatedInverseClass::WriteInfo (IOSectionClass &out)
{
  out.WriteVar ("Type", "FREE_PARTICLE");
}


