#include <boost/random.hpp>
#include <boost/numeric/ublas/matrix_sparse.hpp>
#include <boost/numeric/ublas/vector_sparse.hpp>
#include <boost/numeric/ublas/vector.hpp>
#include <boost/numeric/ublas/matrix_expression.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/iterator.hpp>
#include <boost/iterator.hpp>

#include <Common/MatrixOps/MatrixOps.h>
#include <Common/IO/IO.h>
#include <Common/Random/Random.h>
#include <blitz/array.h>
#include "Common.h"
#include "Observables/ObservableBase.h"
#include <list>
using namespace boost::numeric::ublas;
#include <blitz/mstruct.h>
#include <fstream.h>


double k=1.0;
dVec Box;

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

void SaveConfiguration(ofstream outFile,Array<dVec,1> &config)
{
  for (int counter=0;counter<config.size();counter++){
    outFile<<config<<" ";
  }
  outFile<<endl;

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
  //  cerr<<"THe config.size is "<<config.size()<<endl;
  double cutoff2=cutoff*cutoff;
  dVec particleLocation=config(particle);
  for (int orbital=0;orbital<config.size();orbital++)
    if (dist2(particleLocation,lattice(orbital))<cutoff2){
      orbitalList.push_back(orbital);
    }
	 
}

void BuildRelevant(Array<dVec,1> &lattice, Array<dVec,1> &config,
		   int ptcl,
		   list<int> &orbitalList,list<int>  &ptclList)
{
  double cutoff=4.5;
  // double cutoff=160.0;
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
  //  cerr<<"check: "<<ptclList.size()<<" "<<orbitalList.size()<<endl;
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

void sparsify_outer_prod_add(sparse_vector<double> &a,matrix_row<sparse_matrix<double> > &br,sparse_matrix<double> &C,double cutoff,double detRatio, int theSize)
{
  matrix_row<sparse_matrix<double> >::const_iterator iterr;
  sparse_vector<double> b(theSize,theSize);
  for (iterr=br.begin();iterr!=br.end();iterr++){
    //    cerr<<"iterr is "<<iterr.index()<<endl;
   b(iterr.index())=*iterr;
  }
  sparse_vector<double>::const_iterator itera;
  sparse_vector<double>::const_iterator iterb;
  for (itera=a.begin();itera!=a.end();itera++)
    for (iterb=b.begin();iterb!=b.end();iterb++){
      //      if (itera.index()==0)
      //	cerr<<*itera<<" "<<*iterb<<" "<<itera.index()<<" "<<iterb.index()<<endl;
      int ia=itera.index();
      int ib=iterb.index();
      double val= *itera * *iterb;
      double newVal=C(ia,ib)-1.0/detRatio*val;
      //      if (abs(val)>cutoff)
      if (abs(newVal)>=cutoff)
	C(ia,ib) =newVal;
      else
	C.erase(ia,ib);  //=0;
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
  //  cerr<<"1"<<endl;
  int numPerDim = (int) ceil (pow(0.5*(double)num, 1.0/3.0)-1.0e-6);
  cerr<<"2"<<endl;
  double delta = box[0] / numPerDim;
  cerr<<"3"<<endl;
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
    // 	fprintf (stderr, "BCC ptcl %d position = [%8.4f %8.4f %8.4f]\n",
    // 		 ptcl, r[0], r[1], r[2]);
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
  //  cerr<<"A"<<endl;
  BuildBCC(Box,numParticles,config);
  //  cerr<<"B"<<endl;
  BuildBCC(Box,numParticles,lattice);
  //  cerr<<"C"<<endl;
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
  //  cerr<<dist2<<endl;
  return exp(-k*dist2);
  

}
void BuildInitialM(Array<double,2> &M,Array<dVec,1> &lattice, Array<dVec,1> &config)
{
   //Build initial matrix
   for (int ptclNum=0;ptclNum<config.size();ptclNum++)
     for (int orbitalNum=0;orbitalNum<lattice.size();orbitalNum++)
       M(ptclNum,orbitalNum)=calcOrbital(ptclNum, orbitalNum,lattice,config);
}


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
	//	cerr<<"Dist: "<<dist<<" "<<grid.End<<endl;
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
    //    CalcMtInverse(M,MtInverse);
    PseudoIdent.resize(size,size);
    //    MakeBlitzIdentity(Ident);
  }
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
  //This should be by particle but there is an off chance
  //that I have it backwards
//   void CheckMeByCol(Array<double,1> &frobNorm2)
//   {
//     assert(frobNorm.size()==M.extent(0));
//     Transpose(MtInverse);
//     MatMult(M,MtInverse,PseudoIdent);
//     Transpose(MtInverse);
//     for (int i=0;i<PseudoIdent.extent(0);i++)
//       for (int j=0;j<PseudoIdent.extent(1);j++)
// 	if (i!=j)
// 	  frobNorm2(i)+=PseudoIdent(i,j)*PseudoIdent(i,j);
// 	else
// 	  frobNorm2(i)+=(1-PseudoIdent(i,j))*(1-PseudoIdent(i,j));
//     return;
    
//   }
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
    //    cerr<<"What is ptcl "<<*ptclUpdateIter<<endl;

      

    //    cerr<<MSmall.extent(0)<<" "<<MSmall.extent(1)<<" "<<NSmallVec.size()
    //	<<" "<<AlmostIdentityVec.size()<<endl;
    MatVecProd(MSmall,NSmallVec,AlmostIdentityVec);
    ptclCounter=-1;
    for (ptclIter=ptclList.begin();ptclIter!=ptclList.end();ptclIter++){
      ptclCounter++;
      if (*ptclIter==*ptclUpdateIter)
	AlmostIdentityVec(ptclCounter)=1.0;
      else
	AlmostIdentityVec(ptclCounter)=0.0;
    }
    //    cerr<<"NSmallVec is "<<NSmallVec<<endl;
    //    cerr<<"AlmostIdentyVec is "<<AlmostIdentityVec<<endl;
    if (MSmall.extent(0)<MSmall.extent(1)) //Not large enough to hold answer
      AlmostIdentityVec.resizeAndPreserve(MSmall.extent(1));
    LinearLeastSquares(MSmall,NSmallVec,AlmostIdentityVec);
    if (MSmall.extent(0)<MSmall.extent(1))
            AlmostIdentityVec.resizeAndPreserve(MSmall.extent(0));
    //    cerr<<"New NSmallVec is  "<<AlmostIdentityVec<<endl;
    NSmallVec.resizeAndPreserve(MSmall.extent(0));
    AlmostIdentityVec.resizeAndPreserve(MSmall.extent(1));
    MatVecProd(MSmall,AlmostIdentityVec,NSmallVec);
    //    cerr<<"New AlmostIdentityVec is "<<NSmallVec<<endl;
    int orbitalCounter=-1;
    for (orbitalIter=orbitalList.begin();orbitalIter!=orbitalList.end();
	 orbitalIter++){
      orbitalCounter++;
      //cerr<<"AAAA "<<MtInverse(*ptclUpdateIter,*orbitalIter)-AlmostIdentityVec(orbitalCounter)<<endl;
      MtInverse(*ptclUpdateIter,*orbitalIter)=AlmostIdentityVec(orbitalCounter);
      //      MtInverse(*orbitalIter,*ptclIter)=AlmostIdentityVec(orbitalCounter);
    }
     }
  }
//   void LinearSolve(int ptcl,Array<dVec,1> &lattice, Array<dVec,1> &config)
//  {
//     list<int> orbitalList;
//     list<int> ptclList;
//     list<int> ptclUpdateList;
//     NearbyParticleToParticle(lattice, config,
// 			     ptcl, 0.1, ptclUpdateList);

//     //    cerr<<"into linear solve"<<endl;
//     //    cerr<<"Box is "<<Box<<endl;
//     list<int>::iterator ptclUpdateIter;
      
//     BuildRelevant(lattice, config,ptcl,orbitalList,ptclList);
//     cerr<<"SizesB: "<<ptclList.size()<<" "<<orbitalList.size()<<endl;
//     Array<double,2> MSmall(ptclList.size(),orbitalList.size());
//     //    Array<double,2> NSmall(orbitalList.size(),ptclList.size());
//     Array<double,2> NSmall(orbitalList.size(),1);
//     list<int>::iterator ptclIter;
//     list<int>::iterator orbitalIter;



//     for (ptclUpdateIter=ptclUpdateList.begin();ptclUpdateIter!=ptclUpdateList.end();ptclUpdateIter++){
//       cerr<<"Ptcl update list is using"<<*ptclUpdateIter<<" "<<ptcl<<" "<<ptclUpdateList.size()<<endl;
//     int ptclCounter=-1;
//     for (ptclIter=ptclList.begin();ptclIter!=ptclList.end();ptclIter++){
//       ptclCounter++;
//       int orbitalCounter=-1;
//       for (orbitalIter=orbitalList.begin();orbitalIter!=orbitalList.end();orbitalIter++){
// 	orbitalCounter++;
// 	MSmall(ptclCounter,orbitalCounter)=M(*ptclIter,*orbitalIter);
// 	if (*ptclIter==*ptclUpdateIter)
// 	  NSmall(orbitalCounter,0)=MtInverse(*ptclIter,*orbitalIter);
//       }
//     }

//     Array<double,1> NSmallVec(MSmall.extent(1));
//     NSmallVec=NSmall(Range::all(),0);
//     Array<double,1> AlmostIdentityVec(MSmall.extent(0));
//     //    cerr<<"What is ptcl "<<ptcl<<endl;

      

//     //    cerr<<MSmall.extent(0)<<" "<<MSmall.extent(1)<<" "<<NSmallVec.size()
//     //	<<" "<<AlmostIdentityVec.size()<<endl;
//     MatVecProd(MSmall,NSmallVec,AlmostIdentityVec);
//     ptclCounter=-1;
//     for (ptclIter=ptclList.begin();ptclIter!=ptclList.end();ptclIter++){
//       ptclCounter++;
//       if (*ptclIter==*ptclUpdateIter)
// 	AlmostIdentityVec(ptclCounter)=1.0;
//       else
// 	AlmostIdentityVec(ptclCounter)=0.0;
//     }
//     //    cerr<<"NSmallVec is "<<NSmallVec<<endl;
//     //    cerr<<"AlmostIdentyVec is "<<AlmostIdentityVec<<endl;
//     if (MSmall.extent(0)<MSmall.extent(1)) //Not large enough to hold answer
//       AlmostIdentityVec.resizeAndPreserve(MSmall.extent(1));
//     LinearLeastSquares(MSmall,NSmallVec,AlmostIdentityVec);
//     if (MSmall.extent(0)<MSmall.extent(1))
//             AlmostIdentityVec.resizeAndPreserve(MSmall.extent(0));
//     //    cerr<<"New NSmallVec is  "<<AlmostIdentityVec<<endl;
//     NSmallVec.resizeAndPreserve(MSmall.extent(0));
//     AlmostIdentityVec.resizeAndPreserve(MSmall.extent(1));
//     MatVecProd(MSmall,AlmostIdentityVec,NSmallVec);
//     //    cerr<<"New AlmostIdentityVec is "<<NSmallVec<<endl;
//     int orbitalCounter=-1;
//     for (orbitalIter=orbitalList.begin();orbitalIter!=orbitalList.end();
// 	 orbitalIter++){
//       orbitalCounter++;
//       //      cerr<<"AAAA "<<MtInverse(ptcl,*orbitalIter)-AlmostIdentityVec(orbitalCounter)<<endl;
//       MtInverse(*ptclUpdateIter,*orbitalIter)=AlmostIdentityVec(orbitalCounter);
//     }
//     }
//   }
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
};

class SparseUpdateInverseClass
{
public:
  double detRatio;
  Array<double,2> &M;
  sparse_matrix<double> MtInverse;
  sparse_vector<double> u;
  sparse_vector<double> MtInverseu;
  sparse_matrix<double> OuterProductTerm;
  double cutoff;
  SparseUpdateInverseClass(Array<double,2> &tempM) :
    M(tempM)
  {
    cutoff=0.0;
  }
  void Initialize(int size)
  {
    MtInverse.resize(size,size);
    MtInverseu.resize(size);
    u.resize(size);
    OuterProductTerm.resize(size,size);
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
    matrix_row<sparse_matrix<double> > mr (MtInverse, movePtcl);
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
	  MtInverse.erase(i,j);
      }
  }
};


int main(int argc, char **argv)
{
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

  BuildInitialM(M,config,lattice);

  DirectSingleParticleUpdateClass Direct(M,u);
  Direct.Initialize(config.size());
  
  DenseUpdateInverseClass Dense(M,u);
  Dense.SetCutoff(cutoff);
  Dense.Initialize(config.size());

  for (int i=0;i<config.size();i++)
    Dense.LinearSolve(i,lattice,config);
  

  SparseUpdateInverseClass Sparse(M);
  Sparse.SetCutoff(cutoff);
  Sparse.Initialize(config.size());
  
  cerr<<Box<<endl;
  cerr<<"Printing orbitals"<<endl;

  for (int i=0;i<lattice.size();i++)
    cerr<<i<<" "<<lattice(i)<<endl;

  for (int mcSteps=0;mcSteps<config.size()*config.size()*config.size();mcSteps++){
    int movePtcl=Random.LocalInt(config.size());
    double ranNumber=Random.Local();
    dVec changeGaussian;
    Random.LocalGaussianVec(sigma,changeGaussian);
    config(movePtcl)=config(movePtcl)+changeGaussian;
    PutInBox(config(movePtcl));
    for (int orbital=0;orbital<config.size();orbital++){
      u(orbital)=calcOrbital(movePtcl,orbital,lattice,config)-M(movePtcl,orbital);
    }
    Sparse.calcU(movePtcl,lattice,config);
    Direct.CalcDeterminantRatio(movePtcl);
    if (mcSteps % 10==0)
      Direct.MapM(config,lattice);
    Dense.CalcDeterminantRatio(movePtcl);


    Sparse.CalcDeterminantRatio(movePtcl);
    //    cerr<<"ratios are "<<Direct.detRatio<<" "<<Dense.detRatio<<"
    //    "<<Sparse.detRatio<<endl;

    if (Direct.detRatio>ranNumber){//accept
      //      cerr<<"ratios are "<<Dense.detRatio<<" "<<Sparse.detRatio<<endl;
      for (int orbital=0;orbital<config.size();orbital++){
	M(movePtcl,orbital)=calcOrbital(movePtcl, orbital,lattice, config);
      }
      //      Dense.UpdateInverse(movePtcl);
      //      if (Dense.detRatio>0.4)
      ///////      Dense.LinearSolve(movePtcl,lattice,config);
      //      Dense.Sparsify();

      //      cerr<<"Percent 0 is "<<(double)Dense.sparseElements/(double)(config.size()*config.size())<<endl;
      if (mcSteps % 1000==0){
	for (int i=0;i<Direct.realHist.size();i++)
	  cerr<<i<<" "<<Direct.realHist(i)/(double)Direct.timesAccumulated<<" "<<Direct.inverseHist(i)/(double)Direct.timesAccumulated<<endl;
      }
      if (1==2){//mcSteps % 1000==0){
	double theCheck=Dense.CheckMe();
	cerr<<"My check is "<<theCheck<<endl;
	if (theCheck>1.0){
	  for (int i=0;i<config.size();i++)
	    Dense.LinearSolve(i,lattice,config);
	  Direct.CalcDeterminantRatio(0);
	  double totalInverseError=0.0;
	  for (int i=0;i<config.size();i++){
	    double inverseError=0.0;
	    for (int j=0;j<config.size();j++){
	      inverseError+=abs(Direct.MtInverse(i,j)-Dense.MtInverse(i,j));
	      totalInverseError+=abs(Direct.MtInverse(i,j)-Dense.MtInverse(i,j));
	    }
	    cerr<<i<<" "<<inverseError<<endl;
	    if (inverseError>1.0){
	      cerr<<"AAAAA"<<endl;
	      for (int j=0;j<config.size();j++){
		double jinverseError=0.0;
		jinverseError+=abs(Direct.MtInverse(i,j)-Dense.MtInverse(i,j));
		cerr<<"J: "<<j<<" "<<jinverseError<<endl;
	      }
	      list<int> orbitalList;
	      list<int> ptclList;
	      BuildRelevant(lattice, config,
			    i,
			    orbitalList,ptclList);
	      list<int>::iterator grrIter;
	      for (grrIter=orbitalList.begin();grrIter!=orbitalList.end();grrIter++)
		cerr<<*grrIter<<" "<<endl;
	      cerr<<"degrr"<<endl;
	      cerr<<Direct.MtInverse(i,Range::all())<<endl;
	    }

	  }
	  if (totalInverseError>1.0){
	    cerr<<"total inverse: "<<totalInverseError<<endl;
	    cerr<<"ptcl moved is "<<movePtcl<<endl;
	    for (int i=0;i<config.size();i++)
	      cerr<<"Config is "<<config(i)<<endl;
	    cerr<<"My check is "<<theCheck<<endl;
	}
	}
	//	SaveConfiguration(outFile,config);
      }
      /////SPARSE      Sparse.UpdateInverse(movePtcl);
    }
    else {
      config(movePtcl)=config(movePtcl)-changeGaussian;
      PutInBox(config(movePtcl));
      if (1==2){ //mcSteps % 1000==0){

	cerr<<"My check is "<<Dense.CheckMe()<<endl;
	//	SaveConfiguration(outFile,config);
      }

    }  

  }
  

  
}
