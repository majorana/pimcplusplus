
#include <boost/random.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
// #include <boost/numeric/ublas/iterator.hpp>
//  #include <boost/numeric/ublas/detail/iterator.hpp>
#include <boost/iterator.hpp>
#include <vector>
#include <Common/MatrixOps/MatrixOps.h>
#include <Common/IO/IO.h>
#include <Common/Random/Random.h>
#include <blitz/array.h>
#include "Common.h"
#include "Observables/ObservableBase.h"
#include <list>

#include <algorithm>
using namespace boost::numeric::ublas;

#include <blitz/mstruct.h>
#include <fstream.h>


double k=1.0;
dVec Box;

class distTuple
{
public:
  double dist;
  int ptcl;
};

// // class mePtrClass
// // {
// //   int nextPtcl;
// //   int prevPtcl
// //   int orbital;
// //   int ptcl;
////}


// class MInfoClass
// {
//   Array<mePtrClass,2> MInfo;
//   Array<double,2> &M;
//   MInfoClass(Array<double,2> &tempM) : M(tempM)
//   {
//   }
//   double Cutoff;
//   void Initialize(double cutoff)
//   {
//     Cutoff=cutoff;
//     for (int orbital=0;orbital<M.extent(1);orbital++){
//       int currLarge=-1;
//       for (int ptcl=0;ptcl<M.extent(0);ptcl++){
// 	MInfo(ptcl,orbital).ptcl=ptcl;
// 	MInfo(ptcl,orbital).orbital=orbital;
// 	if (abs(M(ptcl,orbital))>cutoff){
// 	  MInfo(ptcl,orbital).prevPtcl=currLarge;
// 	  if (currLarge!=-1)
// 	    MInfo(currLarge,orbital).nextPtcl=ptcl;
// 	  currLarge=ptcl;
// 	}
//       }
//     }
//   }
//   void Next(int ptcl,int orbital)
//   {
//     assert(abs(M(ptcl,orbital))>Cutoff);
//     return MInfo(ptcl,orbital).nextPtcl;
//   }
//   void Prev(int ptcl,int orbital)
//   {
//     assert(abs(M(ptcl,orbital))>Cutoff);
//     return MInfo(ptcl,orbital).prevPtcl;
//   }
  
  

// };

// class MInfoClass
// {
//   Array<double,2> &M;
//   Array<list<int>,1> orbitals;
//   void Initialize()
//   {
//     for (int ptcl=0;ptcl<M.extent(0);ptcl++)
//       for (int orbital=0;orbital<M.extent(1);orbital++)
// 	M(ptcl,orbital)=

//   }
//   MInfoClass(Array<double,2> &tempM) : M(tempM){}
  



// }
void 
PutInBox (dVec &v)
{
  for (int i=0; i<NDIM; i++) {
    double n = -floor(v(i)*(1.0/Box(i))+0.5);
    v(i) += n*Box(i);
  }
  
}

double dist2(dVec x,dVec y)
{
  dVec diff=x-y;
  PutInBox(diff);
  return dot(diff,diff);

}

bool operator<(const distTuple& a, const distTuple& b) {
  return a.dist < b.dist;
}


void SortParticles(Array<dVec,1> &config, dVec sortFrom,
		   std::vector<distTuple> &distArray)
{
  distArray.resize(config.size());
  for (int i=0;i<distArray.size();i++){
    distArray[i].ptcl=i;
    distArray[i].dist=sqrt(dist2(sortFrom,config(i)));
  }
  sort(distArray.begin(),distArray.end());

}

		   
	     

void NearbyParticle(Array<dVec,1> &lattice, Array<dVec,1> &config,
		    int orbital, double cutoff, list<int> &ptclList)
{
  double cutoff2=cutoff*cutoff;
  dVec orbitalLocation=lattice(orbital);
  for (int ptcl=0;ptcl<config.size();ptcl++)
    if (dist2(orbitalLocation,config(ptcl))<cutoff2){
      ptclList.push_back(ptcl);
    }
	 
}

void SaveConfiguration(ofstream &outFile,Array<dVec,1> &config)
{
  //  cerr<<"Beginning config write"<<endl;
  for (int counter=0;counter<config.size();counter++){
    for (int dim=0;dim<3;dim++){
      outFile<<config(counter)[dim]<<" ";
      //      cerr<<config(counter)[dim]<<" ";
    }
    //    cerr<<endl;
    outFile<<endl;
  }
  outFile<<endl;
  //  cerr<<endl;
  //  cerr<<"Ending config write"<<endl;
}

void NearbyParticleToParticle(Array<dVec,1> &lattice, Array<dVec,1> &config,
			      int ptcl, double cutoff, list<int> &ptclList)
{
  double cutoff2=(cutoff)*(cutoff);
  dVec ptclLocation=config(ptcl);
  for (int ptclA=0;ptclA<config.size();ptclA++)
    if (dist2(ptclLocation,config(ptclA))<cutoff2){
      ptclList.push_back(ptclA);
    }

}


void NearbyOrbital(Array<dVec,1> &lattice, Array<dVec,1> &config,
		    int particle, double cutoff, list<int> &orbitalList)
{
  double cutoff2=cutoff*cutoff;
  dVec particleLocation=config(particle);
  for (int orbital=0;orbital<config.size();orbital++)
    if (dist2(particleLocation,lattice(orbital))<cutoff2){
      orbitalList.push_back(orbital);
    }
	 
}

void BuildRelevant2(Array<dVec,1> &lattice, Array<dVec,1> &config,
		    int ptcl,
		    list<int> &orbitalList, list<int> &ptclList,
		    Array<double,2> &M)
{
  double cutoff=2.0;
  double MCutoff=0.01;
  //First choose relevant orbitals. 
  NearbyOrbital(lattice, config,ptcl, cutoff, orbitalList);
  //Now, you want to choose the particles that are relevant to the orbitals
  //you've chosen. There must be more particles then orbitals.  Ideally
  //you want to choose particles that have non-zero values at these orbitals
  list<int> tempPtclList;
  list<int>::iterator orbitalIter;
  for (int ptcl=0;ptcl<config.size();ptcl++)
    for (orbitalIter=orbitalList.begin();orbitalIter!=orbitalList.end();
	 orbitalIter++)
      if (abs(M(ptcl,*orbitalIter))>MCutoff)
	tempPtclList.push_back(ptcl);


  
  tempPtclList.sort();
  int currPtcl=-1;
  list<int>::iterator ptclIter;
  for (ptclIter=tempPtclList.begin();ptclIter!=tempPtclList.end();
       ptclIter++){
    if (currPtcl!=*ptclIter)
      ptclList.push_back(*ptclIter);
    currPtcl=*ptclIter;
  }
  //  cerr<<"The total number of SPAI orbitals is "<<orbitalList.size()<<endl;
  //  cerr<<"The total number of SPAI particles is "<<ptclList.size()<<endl;
}


void BuildRelevant(Array<dVec,1> &lattice, Array<dVec,1> &config,
		   int ptcl,
		   list<int> &orbitalList,list<int>  &ptclList)
{
  //  double cutoff=4.5;
  double cutoff=160.0;
  list<int> tempPtclList;
  NearbyOrbital(lattice, config,ptcl, cutoff, orbitalList);
  NearbyParticleToParticle(lattice,config,ptcl,cutoff,ptclList);
  double tempCutoff=cutoff;
  ///  while (ptclList.size()<=orbitalList.size()+10 &&
  while (ptclList.size()<config.size()){
    ptclList.clear();
      tempCutoff=tempCutoff+0.05*Box[0];
    NearbyParticleToParticle(lattice,config,ptcl,tempCutoff,ptclList);
  }

  return;
  list<int>::iterator iter;
  for (iter=orbitalList.begin();iter!=orbitalList.end();iter++)
    NearbyParticle(lattice,config,*iter,cutoff,tempPtclList);
  tempPtclList.sort();
  int currPtcl=-1;
  list<int>::iterator ptclIter;
  for (ptclIter=tempPtclList.begin();ptclIter!=tempPtclList.end();
       ptclIter++){
    if (currPtcl!=*ptclIter)
      ptclList.push_back(*ptclIter);
    currPtcl=*ptclIter;
  }
}
     
    
 


void MakeBlitzIdentity(Array<double,2> &Id)
{
  for (int i=0;i<Id.extent(0);i++)
    for (int j=0;j<Id.extent(1);j++)
      if (i==j)
	Id(i,j)=1.0;
      else
	Id(i,j)=0.0;
}

void x()
{
  cerr<<"Hi"<<endl;
  mapped_vector<double> x(3,3);
}

void sparsify_outer_prod_add(mapped_vector<double> &a,matrix_row<coordinate_matrix<double> > &br,coordinate_matrix<double> &C,double cutoff,double detRatio, int theSize)
{
  matrix_row<coordinate_matrix<double> >::const_iterator iterr;
  mapped_vector<double> b(theSize,theSize);
  for (iterr=br.begin();iterr!=br.end();iterr++){
   b(iterr.index())=*iterr;
  }
  mapped_vector<double>::const_iterator itera;
  mapped_vector<double>::const_iterator iterb;
  for (itera=a.begin();itera!=a.end();itera++)
    for (iterb=b.begin();iterb!=b.end();iterb++){
      int ia=itera.index();
      int ib=iterb.index();
      double val= *itera * *iterb;
      double newVal=C(ia,ib)-1.0/detRatio*val;
      if (abs(newVal)>=cutoff)
	C(ia,ib) =newVal;
      else
	C.erase_element(ia,ib);  
    }
}


void BuildBCC(dVec &box,int numParticles,Array<dVec,1> &config)
{
  int num = numParticles;
  config.resize(num);
  bool isCubic = (box[0]==box[1]) && (box[1]==box[2]);
  if (!isCubic) {
    perr << "A cubic box is current required for cubic initilization\n";
    abort();
  }
  int numPerDim = (int) ceil (pow(0.5*(double)num, 1.0/3.0)-1.0e-6);
  double delta = box[0] / numPerDim;
  for (int ptcl=0; ptcl<config.size(); ptcl++) {
    int ip = (ptcl)/2;
    int ix, iy, iz;
    ix = ip/(numPerDim*numPerDim);
    iy = (ip-(ix*numPerDim*numPerDim))/numPerDim;
    iz = ip - ix*numPerDim*numPerDim - iy*numPerDim;
    dVec r;
    r[0] = ix*delta-0.5*box[0];
    r[1] = iy*delta-0.5*box[1];
    r[2] = iz*delta-0.5*box[2];
    if (ptcl % 2) 
      r += 0.5*delta;
    config(ptcl)=r;
  }
}

void Read(IOSectionClass &in, Array<dVec,1> &lattice,Array<dVec,1> &config,
	  double &cutoff)
{
  int numParticles;
  assert(in.ReadVar("NumParticles",numParticles));
  assert(in.ReadVar("k",k));
  assert(in.ReadVar("cutoff",cutoff));
  Box(0)=pow(numParticles*4.0/3.0*M_PI,1.0/3.0);
  Box(1)=pow(numParticles*4.0/3.0*M_PI,1.0/3.0);
  Box(2)=pow(numParticles*4.0/3.0*M_PI,1.0/3.0);
  BuildBCC(Box,numParticles,config);
  BuildBCC(Box,numParticles,lattice);
  dVec toAdd;
  toAdd(0)=0.2;
  toAdd(1)=0.3;
  toAdd(2)=0.1;
  for (int i=0;i<lattice.size();i++)
    config(i)=config(i)+toAdd;
  return;


  Array<double,2> tempLattice;
  assert(in.ReadVar("Lattice",tempLattice));
  lattice.resize(tempLattice.extent(0));
  for (int dim=0;dim<NDIM;dim++)
    for (int i=0;i<lattice.size();i++)
      lattice(i)[dim]=tempLattice(i,dim);

  Array<double,2> tempConfig;
  
  assert(in.ReadVar("InitialConfig",tempConfig));
  config.resize(tempConfig.extent(0));
  for (int dim=0;dim<NDIM;dim++)
    for (int i=0;i<config.size();i++)
      config(i)[dim]=tempConfig(i,dim);
  assert(lattice.size()==config.size());
}

void CalcMtInverse(Array<double,2> &M,Array<double,2> &MtInverse)
{
    //M^{-1t}
    MtInverse=M;
    Transpose(MtInverse);
    MtInverse=Inverse(MtInverse);
}

void SetBox(double x,double y, double z)
{
  Box(0)=x;
  Box(1)=y;
  Box(2)=z;
}
double calcOrbital(int ptclNum, int orbitalNum,Array<dVec,1> &lattice, Array<dVec,1> &config)
{
  dVec diff=config(ptclNum)-lattice(orbitalNum);
  PutInBox(diff);
  double dist2=dot(diff,diff);
  return exp(-k*dist2);
  

}
void BuildInitialM(Array<double,2> &M,Array<dVec,1> &lattice, Array<dVec,1> &config)
{
   //Build initial matrix
   for (int ptclNum=0;ptclNum<config.size();ptclNum++)
     for (int orbitalNum=0;orbitalNum<lattice.size();orbitalNum++)
       M(ptclNum,orbitalNum)=calcOrbital(ptclNum, orbitalNum,lattice,config);
}

class TruncatedClass
{
public:
  double detRatio;
  Array<double,2> &M;
  Array<double,2> MSmall;
  Array<double,1> &u;
  std::vector<distTuple> ptclDistArray;
  std::vector<distTuple> orbitalDistArray;
  void Initialize(int size)
  {
    ///do nothing for now
  }
  void BuildSmallOld(Array<dVec,1> &config, Array<dVec,1> &lattice,
		  int movedPtcl, int numParticles)
  {
    SortParticles(config, config(movedPtcl),ptclDistArray);
    SortParticles(lattice, config(movedPtcl),orbitalDistArray);
    MSmall.resize(numParticles,numParticles);
    for (int i=0;i<MSmall.extent(0);i++)
      for (int j=0;j<MSmall.extent(1);j++){
	MSmall(i,j)=M(ptclDistArray[i].ptcl,orbitalDistArray[j].ptcl);
      }
  }
  
  void BuildSmallNew(Array<dVec,1> &config, Array<dVec,1> &lattice,
		     int movedPtcl)
  {
    ///We assumed the moved particle always corresponds to the 0'th particle
    for (int j=0;j<MSmall.extent(1);j++){
      MSmall(0,j)=M(movedPtcl,orbitalDistArray[j].ptcl)+u(orbitalDistArray[j].ptcl);
    }
  }
  void CalcDeterminantRatio(Array<dVec,1> &config, Array<dVec,1> &lattice,
			       int movedPtcl, int numParticles)
  {
    BuildSmallOld(config,lattice,movedPtcl,numParticles);
    double oldDet=Determinant(MSmall);
    BuildSmallNew(config,lattice,movedPtcl);
    double newDet=Determinant(MSmall);
    detRatio=newDet/oldDet;
  }
  TruncatedClass(Array<double,2> &tempM, Array<double,1> &tempu) :
    M(tempM), u(tempu)
  {
  }

};

class DirectSingleParticleUpdateClass
{
public:
  double detRatio;
  Array<double,2> &M;
  Array<double,2> MtInverse;
  Array<double,1> &u;
  Array<double,1> MtInverseu;
  LinearGrid grid;
  Array<double,1> realHist;
  Array<double,1> inverseHist;
  int timesAccumulated;

  DirectSingleParticleUpdateClass(Array<double,2> &tempM, Array<double,1> &tempu) :
    M(tempM), u(tempu)
  {

  }
  void Initialize(int size)
  {
    MtInverse.resize(size,size);
    MtInverseu.resize(size);
    cerr<<"My box is "<<Box<<endl;
    grid.Init(0.0,Box[0]*2.0,size*2);
    timesAccumulated=0;
    realHist.resize(2*size);
    realHist=0.0;
    inverseHist.resize(2*size);
    inverseHist=0.0;
  }
  void CalcDeterminantRatio(int movePtcl)
  {
    CalcMtInverse(M,MtInverse);
    MatVecProd(MtInverse,u,MtInverseu);
    detRatio=1+MtInverseu(movePtcl); //1+M^t(-1)e_k^T
  }



  // (Ignore rest of class for time being)
  void MapM(Array<dVec,1> &config, Array<dVec,1> &lattice)
  {
    Array<double,2> MSquared(config.size(),config.size());
    Array<double,2> MInverseSquared(config.size(),config.size());
    //    MatMult(M,M,MSquared);
    //    MatMult(MtInverse,MtInverse,MInverseSquared);
    for (int i=0;i<config.size();i++){
      for (int j=0;j<config.size();j++){
	dVec diff=config(i)-config(j);
	PutInBox(diff);
	double dist2=dot(diff,diff);
	MSquared(i,j)=exp(-k*dist2);
      }
    }
    MInverseSquared=Inverse(MSquared);
    for (int i=0;i<config.size();i++)
      for (int j=0;j<config.size();j++){
	dVec ptcl=config(i);
	//	dVec orbital=lattice(j);
	dVec orbital=config(j);
	double dist=sqrt(dist2(ptcl,orbital));
	timesAccumulated++;
	if (dist<grid.End) {
	  int index=grid.ReverseMap(dist);
	  //	  cerr<<"The index is "<<index<<endl;
	  //	  cerr<<"The realHist size is "<<realHist.size()<<endl;
	  realHist(index)=realHist(index)+MSquared(i,j);
	  inverseHist(index)=inverseHist(index)+MInverseSquared(i,j);
	  //	  cerr<<"R: "<<realHist(index)<<endl;
	} 
      }
  }


  
};

class SPAIClass
{
public:
  double detRatio;
  int sparseElements;
  Array<double,2> &M;
  Array<double,2> MtInverse;
  Array<double,1> &u;
  Array<double,1> MtInverseu;
  Array<double,2> OuterProductTerm;
  Array<double,2> PseudoIdent;
  double cutoff;
  SPAIClass(Array<double,2> &tempM, Array<double,1> &tempu) :
    M(tempM), u(tempu)
  {
    cutoff=0.0;
  }
  void LinearSolve(int ptcl,Array<dVec,1> &lattice, Array<dVec,1> &config)
  {
    list<int> orbitalList;
    list<int> ptclList;
    list<int> ptclUpdateList;
    NearbyParticleToParticle(lattice, config,
    			     ptcl, 0.0001, ptclUpdateList);
    list<int>::iterator ptclUpdateIter;
    /////    cerr<<"Check: The number of columns I'm updating simultaneously is "<<ptclUpdateList.size()<<endl;
    //    cerr<<"into linear solve"<<endl;
    //    cerr<<"Box is "<<Box<<endl;
    for (ptclUpdateIter=ptclUpdateList.begin();ptclUpdateIter!=ptclUpdateList.end();ptclUpdateIter++){
      ptclList.clear();
      orbitalList.clear();
       
      BuildRelevant2(lattice, config,*ptclUpdateIter,orbitalList,ptclList,M);
      /////      cerr<<"SizesB: "<<ptclList.size()<<" "<<orbitalList.size()<<endl;
      Array<double,2> MSmall(ptclList.size(),orbitalList.size());
      //    Array<double,2> NSmall(orbitalList.size(),ptclList.size());
      Array<double,2> NSmall(orbitalList.size(),1);
      list<int>::iterator ptclIter;
      list<int>::iterator orbitalIter;
      

    int ptclCounter=-1;

    for (ptclIter=ptclList.begin();ptclIter!=ptclList.end();ptclIter++){
      ptclCounter++;
      int orbitalCounter=-1;
      for (orbitalIter=orbitalList.begin();orbitalIter!=orbitalList.end();orbitalIter++){
	orbitalCounter++;
	MSmall(ptclCounter,orbitalCounter)=M(*ptclIter,*orbitalIter);
	if (*ptclIter==*ptclUpdateIter)
	  NSmall(orbitalCounter,0)=MtInverse(*ptclIter,*orbitalIter);
	//	  NSmall(orbitalCounter,0)=MtInverse(*orbitalIter,*ptclIter);
      }
    }

    Array<double,1> NSmallVec(MSmall.extent(1));
    NSmallVec=NSmall(Range::all(),0);
    Array<double,1> AlmostIdentityVec(MSmall.extent(0));
    MatVecProd(MSmall,NSmallVec,AlmostIdentityVec);
    ptclCounter=-1;
    for (ptclIter=ptclList.begin();ptclIter!=ptclList.end();ptclIter++){
      ptclCounter++;
      if (*ptclIter==*ptclUpdateIter)
	AlmostIdentityVec(ptclCounter)=1.0;
      else
	AlmostIdentityVec(ptclCounter)=0.0;
    }
    if (MSmall.extent(0)<MSmall.extent(1)) //Not large enough to hold answer
      AlmostIdentityVec.resizeAndPreserve(MSmall.extent(1));
    LinearLeastSquares(MSmall,NSmallVec,AlmostIdentityVec);
    if (MSmall.extent(0)<MSmall.extent(1))
            AlmostIdentityVec.resizeAndPreserve(MSmall.extent(0));
    NSmallVec.resizeAndPreserve(MSmall.extent(0));
    AlmostIdentityVec.resizeAndPreserve(MSmall.extent(1));
    MatVecProd(MSmall,AlmostIdentityVec,NSmallVec);
    int orbitalCounter=-1;
    for (orbitalIter=orbitalList.begin();orbitalIter!=orbitalList.end();
	 orbitalIter++){
      orbitalCounter++;
      MtInverse(*ptclUpdateIter,*orbitalIter)=AlmostIdentityVec(orbitalCounter);
      //      MtInverse(*orbitalIter,*ptclIter)=AlmostIdentityVec(orbitalCounter);
    }
     }
  }

  void Initialize(int size,
		  Array<dVec,1> &lattice,Array<dVec,1> &config)
  {
    MtInverse.resize(size,size);
    MtInverseu.resize(size);
    OuterProductTerm.resize(size,size);
    //    CalcMtInverse(M,MtInverse);
    PseudoIdent.resize(size,size);
    //    MakeBlitzIdentity(Ident);
    for (int i=0;i<config.size();i++)
      LinearSolve(i,lattice,config);
  }
  ///for checking code
  double CheckMe()
  {
    Transpose(M);
    MatMult(M,MtInverse,PseudoIdent);
    Transpose(M);
    double frobNorm2=0.0;
    for (int i=0;i<PseudoIdent.extent(0);i++)
      for (int j=0;j<PseudoIdent.extent(1);j++)
	if (i!=j)
	  frobNorm2+=PseudoIdent(i,j)*PseudoIdent(i,j);
	else
	  frobNorm2+=(1-PseudoIdent(i,j))*(1-PseudoIdent(i,j));
    return sqrt(frobNorm2);
    
  }
   void CalcDeterminantRatio(int movePtcl)
  {
    MatVecProd(MtInverse,u,MtInverseu);
    detRatio=1+MtInverseu(movePtcl); //1+M^t(-1)e_k^T
  }

  //Must call calcDeterminantRatio before this is called
  //so that the detRatio and MtInverseu is correct!!!!
  ////  void UpdateInverse(int movePtcl)
  ////  {
  ////    OuterProduct(MtInverseu,MtInverse(movePtcl,Range::all()),OuterProductTerm);
  ////    MtInverse=MtInverse-1.0/detRatio*OuterProductTerm;
  ////  }
  void SetCutoff(double tempcutoff)
  {
    cutoff=tempcutoff;
  }
  ///Sparsify calculates how sparse things are
  void Sparsify()
  {
    sparseElements=0;
    for (int i=0;i<MtInverse.extent(0);i++)
      for (int j=0;j<MtInverse.extent(1);j++)
	if (abs(MtInverse(i,j))<cutoff){
	  MtInverse(i,j)=0.0;
	  sparseElements++;
	}

  }
  void UpdateInverse(Array<dVec,1> &lattice, Array<dVec,1> &config)
  {
    //currently we are updating particles and we are updating all of them.
    for (int i=0;i<config.size();i++)
      LinearSolve(i,lattice,config);
  }

};



class DenseUpdateInverseClass
{
public:
  double detRatio;
  int sparseElements;
  Array<double,2> &M;
  Array<double,2> MtInverse;
  Array<double,1> &u;
  Array<double,1> MtInverseu;
  Array<double,2> OuterProductTerm;
  Array<double,2> PseudoIdent;
  double cutoff;
  DenseUpdateInverseClass(Array<double,2> &tempM, Array<double,1> &tempu) :
    M(tempM), u(tempu)
  {
    cutoff=0.0;
  }
  void Initialize(int size)
  {
    MtInverse.resize(size,size);
    MtInverseu.resize(size);
    OuterProductTerm.resize(size,size);
    CalcMtInverse(M,MtInverse);
    PseudoIdent.resize(size,size);
    //    MakeBlitzIdentity(Ident);
  }
  ///for checking code
  double CheckMe()
  {
    Transpose(M);
    MatMult(M,MtInverse,PseudoIdent);
    Transpose(M);
    double frobNorm2=0.0;
    for (int i=0;i<PseudoIdent.extent(0);i++)
      for (int j=0;j<PseudoIdent.extent(1);j++)
	if (i!=j)
	  frobNorm2+=PseudoIdent(i,j)*PseudoIdent(i,j);
	else
	  frobNorm2+=(1-PseudoIdent(i,j))*(1-PseudoIdent(i,j));
    return sqrt(frobNorm2);
    
  }
   void CalcDeterminantRatio(int movePtcl)
  {
    MatVecProd(MtInverse,u,MtInverseu);
    detRatio=1+MtInverseu(movePtcl); //1+M^t(-1)e_k^T
  }

  //Must call calcDeterminantRatio before this is called
  //so that the detRatio and MtInverseu is correct!!!!
  void UpdateInverse(int movePtcl)
  {
    OuterProduct(MtInverseu,MtInverse(movePtcl,Range::all()),OuterProductTerm);
    MtInverse=MtInverse-1.0/detRatio*OuterProductTerm;
  }
  void SetCutoff(double tempcutoff)
  {
    cutoff=tempcutoff;
  }
  ///Sparsify calculates how sparse things are
  void Sparsify()
  {
    sparseElements=0;
    for (int i=0;i<MtInverse.extent(0);i++)
      for (int j=0;j<MtInverse.extent(1);j++)
	if (abs(MtInverse(i,j))<cutoff){
	  MtInverse(i,j)=0.0;
	  sparseElements++;
	}

  }
  ///For SPAI (ignore rest of class for time being)
  void LinearSolve(int ptcl,Array<dVec,1> &lattice, Array<dVec,1> &config)
  {
    list<int> orbitalList;
    list<int> ptclList;
    list<int> ptclUpdateList;
    NearbyParticleToParticle(lattice, config,
    			     ptcl, 0.0001, ptclUpdateList);
    list<int>::iterator ptclUpdateIter;
    //    cerr<<"Check: The number of columns I'm updating simultaneously is "<<ptclUpdateList.size()<<endl;
    //    cerr<<"into linear solve"<<endl;
    //    cerr<<"Box is "<<Box<<endl;
    for (ptclUpdateIter=ptclUpdateList.begin();ptclUpdateIter!=ptclUpdateList.end();ptclUpdateIter++){
      ptclList.clear();
      orbitalList.clear();
      
      BuildRelevant(lattice, config,*ptclUpdateIter,orbitalList,ptclList);
      //    cerr<<"SizesB: "<<ptclList.size()<<" "<<orbitalList.size()<<endl;
      Array<double,2> MSmall(ptclList.size(),orbitalList.size());
      //    Array<double,2> NSmall(orbitalList.size(),ptclList.size());
      Array<double,2> NSmall(orbitalList.size(),1);
      list<int>::iterator ptclIter;
      list<int>::iterator orbitalIter;
      
      
      int ptclCounter=-1;
      
      for (ptclIter=ptclList.begin();ptclIter!=ptclList.end();ptclIter++){
	ptclCounter++;
	int orbitalCounter=-1;
	for (orbitalIter=orbitalList.begin();orbitalIter!=orbitalList.end();orbitalIter++){
	  orbitalCounter++;
	  MSmall(ptclCounter,orbitalCounter)=M(*ptclIter,*orbitalIter);
	  if (*ptclIter==*ptclUpdateIter)
	    NSmall(orbitalCounter,0)=MtInverse(*ptclIter,*orbitalIter);
	  //	  NSmall(orbitalCounter,0)=MtInverse(*orbitalIter,*ptclIter);
	}
      }
      
      Array<double,1> NSmallVec(MSmall.extent(1));
      NSmallVec=NSmall(Range::all(),0);
      Array<double,1> AlmostIdentityVec(MSmall.extent(0));
      MatVecProd(MSmall,NSmallVec,AlmostIdentityVec);
      ptclCounter=-1;
      for (ptclIter=ptclList.begin();ptclIter!=ptclList.end();ptclIter++){
	ptclCounter++;
	if (*ptclIter==*ptclUpdateIter)
	  AlmostIdentityVec(ptclCounter)=1.0;
	else
	  AlmostIdentityVec(ptclCounter)=0.0;
      }
      if (MSmall.extent(0)<MSmall.extent(1)) //Not large enough to hold answer
	AlmostIdentityVec.resizeAndPreserve(MSmall.extent(1));
      LinearLeastSquares(MSmall,NSmallVec,AlmostIdentityVec);
      if (MSmall.extent(0)<MSmall.extent(1))
	AlmostIdentityVec.resizeAndPreserve(MSmall.extent(0));
      NSmallVec.resizeAndPreserve(MSmall.extent(0));
      AlmostIdentityVec.resizeAndPreserve(MSmall.extent(1));
      MatVecProd(MSmall,AlmostIdentityVec,NSmallVec);
      int orbitalCounter=-1;
      for (orbitalIter=orbitalList.begin();orbitalIter!=orbitalList.end();
	   orbitalIter++){
	orbitalCounter++;
	MtInverse(*ptclUpdateIter,*orbitalIter)=AlmostIdentityVec(orbitalCounter);
	//      MtInverse(*orbitalIter,*ptclIter)=AlmostIdentityVec(orbitalCounter);
      }
    }
  }
  
};

class SparseUpdateInverseClass
{
public:
  double detRatio;
  Array<double,2> &M;
  coordinate_matrix<double> MtInverse;
  mapped_vector<double> u;
  mapped_vector<double> MtInverseu;
  coordinate_matrix<double> OuterProductTerm;
  double cutoff;
  SparseUpdateInverseClass(Array<double,2> &tempM) :
    M(tempM)
  {
    cutoff=0.0;
  }
  void Initialize(int size)
  {
    MtInverse.resize(size,size,false);
    MtInverseu.resize(size);
    u.resize(size);
    OuterProductTerm.resize(size,size,false);
    Array<double,2> tempMtInverse(size,size);
    CalcMtInverse(M,tempMtInverse);
    for (int i=0;i<tempMtInverse.extent(0);i++)
      for (int j=0;j<tempMtInverse.extent(1);j++)
	if (abs(tempMtInverse(i,j))>cutoff)
	  MtInverse(i,j)=tempMtInverse(i,j);
	else
	  MtInverse(i,j)=0.0;
  }
  ///These could be done faster if you just wanted the determinant
  ///because you really only need to do one of the rows to calculate
  ///this. If you need to update the inverse you will need the
  ///information later though
  void CalcDeterminantRatio(int movePtcl)
  {
    MtInverseu.clear();
    MtInverseu= prec_prod(MtInverse,u);
    detRatio=1+MtInverseu(movePtcl); //1+M^t(-1)e_k^T
  }
  void calcU(int movePtcl,Array<dVec,1> &lattice, Array<dVec,1> &config)
  {
    for (int orbital=0;orbital<config.size();orbital++){
      u(orbital)=calcOrbital(movePtcl,orbital,lattice,config)-M(movePtcl,orbital);
    }
  }
  //Must call calcDeterminantRatio before this is called
  //so that the detRatio and MtInverseu is correct!!!!
  void UpdateInverse(int movePtcl)
  {
    matrix_row<coordinate_matrix<double> > mr (MtInverse, movePtcl);
    sparsify_outer_prod_add(MtInverseu,mr,MtInverse,cutoff,detRatio,MtInverse.size1());
  }
  void SetCutoff(double tempcutoff)
  {
    cutoff=tempcutoff;
  }
  ///Could be done much faster because you wouldn't have to be
  //iterating over all the elements, just the elements in the matrix
  void Sparsify()
  {
    for (int i=0;i<MtInverse.size1();i++)
      for (int j=0;j<MtInverse.size2();j++){
	double newVal=MtInverse(i,j);
	if (abs(newVal)>=cutoff)
	  MtInverse(i,j) =newVal;
	else
	  MtInverse.erase_element(i,j);
      }
  }
};


int main(int argc, char **argv)
{
  cerr<<"Starting code "<<endl;
  CommunicatorClass comm;
  RandomClass Random(comm);
  int Seed=Random.InitWithRandomSeed(1);
  SetBox(1.0,1.0,1.0);
  Array<dVec,1> config; Array<dVec,1> lattice;
  IOSectionClass in;
  assert (in.OpenFile(argv[1]));
  double cutoff;
  Read(in,lattice,config,cutoff); 
  double sigma=Box(0)*0.02;
  ofstream outFile;
  outFile.open("configs.out");
  Array<double,2> M(config.size(),config.size());
  Array<double,1> u(config.size());
  cerr<<"Building up initial M"<<endl;
  BuildInitialM(M,config,lattice);

  //Initialize truncated method
  TruncatedClass Truncated(M,u);
  Truncated.Initialize(config.size());

  //Initialize spai method
  SPAIClass Spai(M,u);
  Spai.Initialize(config.size(),lattice,config);

  //Initialize direct single particle update method
  /////  DirectSingleParticleUpdateClass Direct(M,u);
  ////  Direct.Initialize(config.size());
  
  //initialize dense update method
  DenseUpdateInverseClass Dense(M,u);
  Dense.SetCutoff(cutoff);
  Dense.Initialize(config.size());
  cerr<<"Dense is initialized"<<endl;
  ///  for (int i=0;i<config.size();i++)
  ///   Dense.LinearSolve(i,lattice,config);
  
  ///initialize sparse update method
  SparseUpdateInverseClass Sparse(M);
  Sparse.SetCutoff(cutoff);
  Sparse.Initialize(config.size());
  cerr<<"Sparse is initialized"<<endl;
  int accept=0;
  int totalCount=0;
  ///Run monte carlo
  int totalSteps=config.size()*config.size()*config.size()+100;
  cerr<<"CHECK Total size is "<<config.size()<<endl;
  cerr<<"CHECK Total size is "<<totalSteps<<endl;

  for (int mcSteps=0;mcSteps<totalSteps;mcSteps++){
    cerr<<"Starting mcStep "<<mcSteps<<endl;
    int movePtcl=Random.LocalInt(config.size());
    double ranNumber=Random.Local();
    dVec changeGaussian;
    Random.LocalGaussianVec(sigma,changeGaussian);
    config(movePtcl)=config(movePtcl)+changeGaussian;
    PutInBox(config(movePtcl));
    ///Calculate u
    for (int orbital=0;orbital<config.size();orbital++){
      u(orbital)=calcOrbital(movePtcl,orbital,lattice,config)-M(movePtcl,orbital);
    }
    ///The sparse matrix class has to calculate it separately because it wants it in a sparse vector
    //SP    Sparse.calcU(movePtcl,lattice,config);

    //Calculate the determinant ratio using all the methods
    //////    Direct.CalcDeterminantRatio(movePtcl);
    Truncated.CalcDeterminantRatio(config,lattice,movePtcl,40);
    Dense.CalcDeterminantRatio(movePtcl);
    //SP    Sparse.CalcDeterminantRatio(movePtcl);
    Spai.CalcDeterminantRatio(movePtcl);
    ///print out determinant ratio
    cerr<<"ratios are "
      //	<<Direct.detRatio<<" "
	<<Dense.detRatio<<" "
	<<Sparse.detRatio<<" "
	<<Truncated.detRatio<<" "
	<<Spai.detRatio<<endl;
    
    if (Dense.detRatio>ranNumber){//accept
      cerr<<"Accept"<<endl;
      accept++;
      totalCount++;
      //      cerr<<"ACCEPT"<<endl;
      //update M which is being shared
      for (int orbital=0;orbital<config.size();orbital++){
	M(movePtcl,orbital)=calcOrbital(movePtcl, orbital,lattice, config);
      }
      //      cerr<<"Accept mid"<<endl;
      Dense.UpdateInverse(movePtcl);
      Dense.Sparsify();
      //      cerr<<"Grr"<<endl;
      //SP      Sparse.UpdateInverse(movePtcl); 
      Spai.UpdateInverse(lattice,config);
      //      cerr<<"Accept end"<<endl;
      //     cerr<<"CHECK IS "<<Spai.CheckMe()<<endl;
      //      cerr<<"Accept ending"<<endl;
    }
    else { //reject
      cerr<<"Reject"<<endl;
      totalCount++;
      //      cerr<<"Reject"<<endl;
      config(movePtcl)=config(movePtcl)-changeGaussian;
      PutInBox(config(movePtcl));
      //      cerr<<"Reject ending"<<endl;
    }  
    cerr<<"CHECK RATIO IS "<<(double)accept/(double)totalCount<<" "
	<<accept<<" "<<totalCount<<endl;

    if (mcSteps % config.size()==0) {
      SaveConfiguration(outFile,config);
      // cerr<<"Hi"<<endl;
    }
      
  }


  
}
