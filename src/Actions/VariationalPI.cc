#include <gsl/gsl_linalg.h>
#include <gsl/gsl_blas.h>

#include "VariationalPI.h"
#include "../PathDataClass.h"
#include <Common/MatrixOps/MatrixOps.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
extern "C"{
#include "../det_calc_uekt.h"
}
// int testme();


//We really want to know what the changed particle is here. Currently
//I'm going to try to just store it when we are doing the singleactionbelow
void 
VariationalPIClass::AcceptCopy(int slice1, int slice2)
{
  cerr<<"I'm being accepted! YEA!"<<endl;
  for (int theRow=0;theRow<DetMatrix.extent(0);theRow++){
    DetMatrix(ChangedColumn,theRow)+=u(theRow);
  }
  for (int col=0;col<DetMatrix.extent(0);col++){
    if (col!=ChangedColumn)
      DetMatrix(col,ChangedColumn)+=u(col);
  }
  //  cerr<<DetMatrix<<endl;
  //  cerr<<"My new determinant is "<<Determinant(DetMatrix)<<endl;
  BuildDeterminantMatrix();
}

void 
VariationalPIClass::RejectCopy(int slice1, int slice2)
{
  cerr<<"I'm being rejected! YEA!"<<endl;
  
}


void 
VariationalPIClass::BuildDeterminantMatrix()
{

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
}




VariationalPIClass::VariationalPIClass (PathDataClass &pathData) :
					    
  NodalActionClass (pathData), 
  Path (pathData.Path)

{
}

void 
VariationalPIClass::Read (IOSectionClass &in)
{
  // Do nothing for now
  //  SetupFreeActions();
  cerr<<"Initializing"<<endl;
  int N=Path.NumParticles();
  cerr<<"the number of particles is"<<N<<endl;
  DetMatrix.resize(N,N);
  u.resize(N);
  BuildDeterminantMatrix();

}


void 
VariationalPIClass::calc_u(const Array<int,1> &changePtcls)
{
  double lambda=PathData.Path.Species(0).lambda;
  double T=1.0/(Path.TotalNumSlices*Path.tau);
  assert(changePtcls.size()==1);
  assert(changePtcls(0)<PathData.Path.NumParticles());
  int ptcl1=changePtcls(0);
  //  cerr<<"ptcl1 is "<<ptcl1<<" "<<u.size()<<endl;
  for (int ptcl2=0;ptcl2<Path.NumParticles();ptcl2++){
    dVec disp=Path(0,ptcl2)-Path(0,ptcl1);
    PathData.Path.PutInBox(disp);
    double dist2=dot(disp,disp);
    u(ptcl2)=exp(-T*dist2/(4.0*lambda))-DetMatrix(ptcl1,ptcl2);
    
  }
}

int matvec(double *x,double *b,int N,void *A){
  
  //  cerr<<"Inside matvec"<<endl;
  //  cerr<<A<<endl;
  //  cerr<<((double**)A)[0]<<endl;
  //  cerr<<b<<endl;
  //  cerr<<x<<endl;
  for (int counter=0;counter<N;counter++){
    b[counter]=0.0;
    for (int counter2=0;counter2<N;counter2++){
      b[counter]=b[counter]+((double*)A)[counter2*N+counter]*x[counter2];
    }
  }
  //  cerr<<"out of matvec"<<endl;
  return 0;
}

void set_gsl_vector_from_double(gsl_vector *V, double* x, int N){
  V->size=N;V->stride=1;V->data=x;V->block=NULL;V->owner=0;  
}


int matvec_gsl_matrix(double *x,double *b,int N,void *params){
  gsl_vector Xgsl,Bgsl;
  set_gsl_vector_from_double(&Xgsl,x,N);
  set_gsl_vector_from_double(&Bgsl,b,N);
  //typedef struct
  //{
  //  size_t size;
  //  size_t stride;
  //  double *data;
  //  gsl_block *block;
  //  int owner;
  //}
  // gsl_vector;

  //int gsl_blas_dgemv (CBLAS_TRANSPOSE_t TransA, double alpha, const gsl_matrix * A, const gsl_vector * x, double beta, gsl_vector * y)  
  int toReturn= gsl_blas_dgemv(CblasNoTrans,1,(gsl_matrix*)params,&Xgsl,0,&Bgsl);  
  //  cerr<<"I am returning "<<toReturn;
  return toReturn;
}/*end matvec_gsl_matrix*/


double 
VariationalPIClass::SingleAction (int startSlice, int endSlice,
				  const Array<int,1> &changePtcls,
				  int level)
{
  //  cerr<<"Calling single action"<<endl;
  //  ModeType currMode=PathData.Path.Path.GetMode();
//   SetMode(OLDMODE);
//   BuildDeterminantMatrix();
//   SetMode(NEWMODE);
  //  cerr<<"Built determinant matrix"<<endl;
  //  cerr<<"My determinanat is "<<Determinant(DetMatrix)<<endl;
  //  DetOrderN DetOrderNInstance;
  drc_uekt_vanilla_parms parms={1e-3,2000};
  calc_u(changePtcls);
  //  cerr<<"calculated change "<<changePtcls(0)<<endl;
  double det_ratio;
  // testme();
//   cerr<<"The address you shoudl look at is "<<DetMatrix.data()<<" "
//       <<u.data()<<endl;
  
//   gsl_matrix *myMatrix;
//   myMatrix=gsl_matrix_calloc(DetMatrix.extent(0),DetMatrix.extent(1));

//   for (int ptcl1=0;ptcl1<PathData.Path.NumParticles();ptcl1++){
//     for (int ptcl2=0;ptcl2<PathData.Path.NumParticles();ptcl2++){
//       gsl_matrix_set(myMatrix,ptcl1,ptcl2,DetMatrix(ptcl1,ptcl2));
//     }
//   }
    
  ChangedColumn=changePtcls(0);

  det_ratio_calculator_uekt_symmetric_vanilla_value(//matvec_gsl_matrix,
						    matvec,
						    DetMatrix.data(), 
						    //(void*)myMatrix,
						    u.size(), 
						    u.data(), changePtcls(0),
    						    0, &det_ratio, 
						    (void*)&parms);
 //  BuildDeterminantMatrix();
  //  cerr<<"My old determinanat is "<<Determinant(DetMatrix)<<endl;
  //  cerr<<"Called ratio calculator"<<endl;

  //Somewhat of a hack
  if (isnan(det_ratio)){
    cerr<<"My det ratio is "<<0.0<<endl;
    return 0.0;
  }
  //  SetMode(currMode);
  
    cerr<<"My det ratio is "<<1.0/abs(det_ratio)<<endl;
  return 1.0/abs(det_ratio);
  //  return 0.0;

}


//just here sot hatt it will compile
bool
VariationalPIClass::IsPositive(int x)
{
  return true;
}

double 
VariationalPIClass::d_dBeta (int slice1, int slice2, int level)
{ 
  return 0.0;
}

NodeType VariationalPIClass::Type()
{
  return FREE_PARTICLE;
}


bool
VariationalPIClass::IsGroundState()
{
  return (false);
}

void
VariationalPIClass::WriteInfo (IOSectionClass &out)
{
  out.WriteVar ("Type", "FREE_PARTICLE");
}

