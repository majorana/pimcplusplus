#ifndef ACTION_CLASS
#define ACTION_CLASS

#include "CubicSpline.h"


class SavedPairActionClass
{



};

/*! $[\frac{u_0(r;\tau)+u_0(r';\tau)}{2}+\sum_{k=1}^n \sum_{j=1}^k
u_{kj}(q;\tau)z^{2j}s^{2(k-j)} */

class PairActionClass
{
	Array<CubicSpline,2> ukj; //(level,f(k,j) )
	double calcUrrptau(double s,double q,double z,int level);
	
	
	  


};

double PairActionClass::calcUrrptau(double s,double q,double z,int level)
{
	double sum=0;
	double r=0.5*(q+z);
	double rprime=0.5*(q-z);
	sum=sum+(ukj(level,0))(r)+(ukj(level,0))(rprime);//this is the endpoint action
	zsquared=z*z;
	ssquared=s*s;
	ssquaredinverse=1/ssquared;
	Sto2k=2;
		for (int k=1;k<=n;k++){
			Zto2j=1;
			currS=Sto2k;
			for (int j=0;j<=k;j++){
				
				cof=(ukj(level,k*(k+1)/2+j))(q); //indexing into the 2darray
				sum=sum+cof*zto2j*currS;
				
				
				zto2j*=zsquared;
				currS=currS*ssquaredinverse;				
			}				
			Sto2k=Sto2k*ssquared;
		}
	
	


}

/*! This is the class that controls all of the actions and is in
  charge of calculating them. When this is initialized a pointer needs
  to be sent that has the memoizedData and IdenticalParticleClass */ 

class ActionClass
{
public:
  Array<PairActionClass,1> PairActionVector;
  Array<int,2> PairMatrix;
  Array<SavedPairActionClass,2> SavedPairActionArray;
  IdenticalParticleClass *myIdenticalParticleClass;
  MemoizedDataClass *myMemoizedDataClass;
  calcTotalAction();
private:


};

