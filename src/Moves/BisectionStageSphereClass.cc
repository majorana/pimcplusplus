#include "BisectionStageSphereClass.h"



void BisectionStageSphereClass::CartesianToSpherical(dVec &r,double &theta,double &phi)
{
  double a=sqrt(r[0]*r[0]+r[1]*r[1]+r[2]*r[2]);
  theta=atan(sqrt(r[1]*r[1]+r[0]*r[0])/r[2]);
  if (r[0]==0 && r[1]==0)
    phi=0;
  else
    phi=atan(r[1]/r[0]);
  if (r[0]<0)
    phi=M_PI+phi;
  if (theta<0)
    theta=theta+M_PI;
  return;
}

TinyVector<double,3> BisectionStageSphereClass::SphericalToCartesian(double &a,double &theta, double &phi)
{
  TinyVector<double,3> r;
  double sinTheta;
  double cosTheta;
  double sinPhi;
  double cosPhi;
  sincos(theta,&sinTheta,&cosTheta);
  sincos(phi,&sinPhi,&cosPhi);
  r[0]=a * sinTheta * cosPhi;
  r[1]=a * sinTheta * sinPhi;
  r[2]=a * cosTheta;
  return r;
}


///Projects the vector r onto the sphere of radius a
void BisectionStageSphereClass::ProjectOntoSphere(dVec &r,double a)
{
  double magR=sqrt(dot(r,r));
  dVec unitr;
  for (int dim=0;dim<NDIM;dim++){
    unitr(dim)=r(dim)/magR;
  }
  r=unitr*a;
  assert(dot(r,r)-a*a<1e-5);

}

void BisectionStageSphereClass::WriteRatio()
{ 
  AcceptRatioVar.Write((double)NumAccepted/(double)NumAttempted);
  OutSection.FlushFile();
}

void BisectionStageSphereClass::Accept()
{
  //do nothing for now
  
}

void BisectionStageSphereClass::Reject()
{
  //do nothing for now

}
double BisectionStageSphereClass::angleInBetween(dVec &r1,dVec &r2)
{
  double magr1=sqrt(dot(r1,r1));
  double magr2=sqrt(dot(r2,r2));
  double angle= acos(dot(r1,r2)/(magr1*magr2));
  //  if (angle<0)
  //    angle=angle+2*M_PI;
  return abs(angle);

}
void BisectionStageSphereClass::RotateAroundZ(dVec &vec, double phi)
{
  phi=-phi;
  if (phi<0)
    phi=phi+2*M_PI;
  double sinPhi;
  double cosPhi;
  sincos(phi,&sinPhi,&cosPhi);
  double x=vec[0];
  double y=vec[1];
  double z=vec[2];
  vec[0]=cosPhi * x + sinPhi * y;
  vec[1]=-sinPhi*x+cosPhi*y;
  vec[2]=z;
}

void BisectionStageSphereClass::RotateAroundX(dVec &vec, double theta)
{
  //  if (theta<0)
  //    theta=theta+2*M_PI;
  double sinTheta;
  double cosTheta;
  sincos(theta,&sinTheta,&cosTheta);
  double x=vec[0];
  double y=vec[1];
  double z=vec[2];
  vec[0]=x;
  vec[1]=cosTheta * y+ sinTheta * z;
  vec[2]=-(sinTheta) * y + cosTheta * z;
}

void BisectionStageSphereClass::RotateAroundY(dVec &vec, double theta)
{
  double sinTheta;
  double cosTheta;
  sincos(theta,&sinTheta,&cosTheta);
  double x=vec[0];
  double y=vec[1];
  double z=vec[2];

  vec[0]=cosTheta * x+ sinTheta * z;
  vec[1]=y;
  vec[2]=-sinTheta * x + cosTheta * z;
}

double BisectionStageSphereClass::Sample(int &slice1,int &slice2,
				   Array<int,1> &activeParticles)
{
  //  cerr<<"I'm being asked to sample"<<endl;
  PathClass &Path = PathData.Path;
  int skip = 1<<(BisectionLevel+1);
  double levelTau = 0.5*PathData.Path.tau*skip;

  int numImages = PathData.Actions.NumImages;

  double logSampleProb=0.0;
  double logOldSampleProb=0.0;

  for (int ptclIndex=0;ptclIndex<activeParticles.size();ptclIndex++){
    int ptcl=activeParticles(ptclIndex);
    double lambda=PathData.Path.ParticleSpecies(ptcl).lambda;
    double sigma2=(1.0*lambda*levelTau);
    double sigma=sqrt(sigma2);
    double prefactorOfSampleProb=0.0; //-NDIM/2.0*log(2*M_PI*sigma2);
    double a=SphereRadius;
    for (int slice=slice1;slice<slice2;slice+=skip){
      SetMode(OLDMODE);
      dVec rOld=Path(slice,ptcl);
      dVec rdiffOld=Path.Velocity(slice,slice+skip,ptcl);
      dVec rbarOld=rOld+ 0.5*rdiffOld;
      ProjectOntoSphere(rbarOld,a);
      dVec rppOld=Path(slice+(skip>>1),ptcl);
      double oldDist;
      oldDist=angleInBetween(rppOld,rbarOld)*a;
      double randomDist;
      double randomTheta;
      double randomPhi;
      double thetaBar,phiBar;
      dVec rpp,rbar;
      SetMode(NEWMODE);
      dVec r=Path(slice,ptcl);
      dVec rdiff=Path.Velocity(slice,slice+skip,ptcl);
      rbar=r+ 0.5*rdiff;

      
      ProjectOntoSphere(rbar,a);
      CartesianToSpherical(rbar,thetaBar,phiBar);
      randomDist=Path.Random.LocalGaussian(sigma);
      randomTheta=abs(randomDist/a);
      randomPhi=Path.Random.Local()*2*M_PI;
      
      
      rpp=SphericalToCartesian(a,randomTheta,randomPhi);
      RotateAroundY(rpp, thetaBar);
      RotateAroundZ(rpp, phiBar);


      /////////////////////////////
      //////////////////////////
      //Test code

      double modRandomTheta;
      modRandomTheta=randomTheta;
      if (abs(angleInBetween(rpp,rbar)-modRandomTheta)>1e-5){
	cerr<<"Angle is between: "<<angleInBetween(rpp,rbar)<<" "<<modRandomTheta<<endl;	
      dVec rbarDup;
      rbarDup[0]=0;
      rbarDup[1]=0;
      rbarDup[2]=a;
      RotateAroundY(rbarDup,thetaBar);
      cerr<<rbar<<" "<<rbarDup<<" "<<thetaBar<<" "<<phiBar<<endl;
      RotateAroundZ(rbarDup,phiBar);
      //      cerr<<rbar<<" "<<rbarDup<<" "<<thetaBar<<" "<<phiBar<<endl;
      cerr<<rbar<<" "<<rbarDup<<" "<<thetaBar<<" "<<phiBar<<endl;
      cerr<<"Random dist:  "<<randomDist<<endl;
      assert(abs(angleInBetween(rpp,rbar)-modRandomTheta)<1e-5);
      }

      ///End Test Code
      /////////////////////
      ////////////
     ///Here we've stored the new position in the path
      Path.SetPos(slice+(skip>>1),ptcl,rpp);
 
 
     
      logSampleProb += prefactorOfSampleProb + 
	-0.5*(randomDist*randomDist)/sigma2;
      logOldSampleProb += prefactorOfSampleProb + 
	-0.5*(oldDist*oldDist)/sigma2;

      
    }
  }

  return exp(-logSampleProb+logOldSampleProb);

}

