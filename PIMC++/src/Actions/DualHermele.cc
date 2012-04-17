
#include "../PathDataClass.h"
#include "Josephson.h"
#include <fftw3.h>

///This has to be called after pathdata knows how many
///particles it has
void DualHermeleClass::Read(IOSectionClass& in)
{
  cerr<<"DualHermele Class read"<<endl;

  assert(in.ReadVar("J",J));
  assert(in.ReadVar("alpha",alpha));
  TotalTime=0;

//   int N=PathData.Path.NumTimeSlices();
//   inArray = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N-1));
//   out=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N-1));
//   p=fftw_plan_dft_1d(PathData.Path.NumTimeSlices()-1,inArray,out,FFTW_BACKWARD,FFTW_MEASURE);


//   inArrayOmega = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N-1));
//   outPhi=(fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N-1));
//   pOmega2Phi=fftw_plan_dft_1d(PathData.Path.NumTimeSlices()-1,inArrayOmega,outPhi,FFTW_FORWARD,FFTW_MEASURE);

  c=1.0;
  K_s=1.0;
  //  omega_c2=(10.0/PathData.Path.tau)*(10.0/PathData.Path.tau);
  omega_c2=100.0;
  //  cerr<<"My omega_c2 is "<<omega_c2<<endl;
  
}

// void DualHermeleClass::Phi2Omega()
// {
//   int N=PathData.Path.NumTimeSlices()-1;
//   for (int slice=0;slice<N;slice++){
//     inArray[slice][0]=Path(slice,0)[0];
//     inArray[slice][1]=0.0;
//   }
//   fftw_execute(p);
//   for (int counter=0;counter<N;counter++){
//     out[counter][0]=out[counter][0]*PathData.Path.tau;
//     out[counter][1]=out[counter][1]*PathData.Path.tau;
//   }
// }

// DualHermeleClass::Omega2Phi()
// {
//   int N=PathData.Path.NumTimeSlices()-1;
//   for (int slice=0;slice<N;slice++){
//     inArrayOmega[slice][0]=Path(slice,0)[0];
//     inArrayOmega[slice][1]=0.0;
//   }
//   fftw_execute(pOmega2Phi);
//   for (int counter=0;counter<N;counter++){
//     out[counter][0]=out[counter][0]*PathData.Path.tau;
//     out[counter][1]=out[counter][1]*PathData.Path.tau;
//   }
// }


// double DualHermeleClass::ComputeFourierAction()
// {
//   double T=1.0/(PathData.Path.tau*(PathData.Path.NumTimeSlices()-1));
//   PathData.Path.Phi2Omega();
//   for (int slice=1;slice<PathData.Path.NumTimeSlices()-1;slice++){
//     double omega_n=2.0*M_PI*slice*T;
//     total+=T*1.0/J*omega_n*omega_n*(out[slice][0]*out[slice][0]+out[slice][1]*out[slice][1]);
//   }
//  return total;
// }




DualHermeleClass::DualHermeleClass(PathDataClass &pathData) : 
  ActionBaseClass (pathData)
{
  

}

double 
DualHermeleClass::SingleAction (int slice1, int slice2,
			    const Array<int,1> &changedParticles,
			    int level)
{
  double total=0.0;
#ifdef ORDER_N_FERMIONS
#if 1==2
  double T=1.0/(PathData.Path.tau*(PathData.Path.NumTimeSlices()-1));
  PathData.Path.Phi2Omega();
  //  for (int slice=0;slice<PathData.Path.NumTimeSlices()-1;slice++)
  //    cerr<<" "<<PathData.Path.outOmega[slice][0]
	
//   PathData.Path.Omega2Phi();
//   for (int slice=0;slice<PathData.Path.NumTimeSlices()-1;slice++)
//     if (abs(PathData.Path(slice,0)[0]-PathData.Path.outPhi[slice][0])>1e-10)
//       cerr<<PathData.Path(slice,0)[0]<<" "<<PathData.Path.outPhi[slice][0]
// 	  <<" "<<PathData.Path.outPhi[slice][1]
// 	  <<" "<<PathData.Path(slice,0)[0]-PathData.Path.outPhi[slice][0]<<endl;
  
  double c=1.0;
  double K_s=1.0;
  double omega_c2=50.0*50.0;

  for (int slice=1;slice<=(PathData.Path.NumTimeSlices()-1)/2;slice++){
    //  for (int slice=1;slice<=(PathData.Path.NumTimeSlices()-1);slice++){
    double omega_n=2.0*M_PI*(slice)*T;
    total+=T*2.0*(1.0/J+1.0*M_PI*M_PI/(4.0*c*K_s)*log((omega_n*omega_n+omega_c2)/(omega_n*omega_n)))
      *omega_n*omega_n*
      (PathData.Path.outOmega[slice][0]*PathData.Path.outOmega[slice][0]+
       PathData.Path.outOmega[slice][1]*PathData.Path.outOmega[slice][1]);
  }
   int slice=(PathData.Path.NumTimeSlices()-1)/2;
     double omega_n=2.0*M_PI*(slice)*T;
     total-=T*1.0*(1.0/J+1.0*M_PI*M_PI/(4.0*c*K_s)*log((omega_n*omega_n+omega_c2)/(omega_n*omega_n)))
       *omega_n*omega_n*
       (PathData.Path.outOmega[slice][0]*PathData.Path.outOmega[slice][0]+
        PathData.Path.outOmega[slice][1]*PathData.Path.outOmega[slice][1]);


//   total-=(T/J+M_PI*M_PI/(4.0*c*K_s)*log((omega_n*omega_n+omega_c2)/(omega_n*omega_n)))*omega_n*omega_n*
//     (PathData.Path.outOmega[slice][0]*PathData.Path.outOmega[slice][0]+
//      PathData.Path.outOmega[slice][1]*PathData.Path.outOmega[slice][1]);

  double totalTau=0.0;
  for (int slice=1;slice<PathData.Path.NumTimeSlices();slice++){
    //    totalTau+=1.0/J*dot(Path(slice,0)-Path(slice-1,0),Path(slice,0)-Path(slice-1,0))/PathData.Path.tau;
    dVec vel=Path(slice,0)-Path(slice-1,0);
    
    totalTau+=1.0/J*(vel[0]*vel[0])/PathData.Path.tau;
  }
  //  cerr<<"T: "<<total<<" "<<totalTau<<endl;
#endif
#endif
  return total;

}


// double 
// DualHermeleClass::SingleAction (int slice1, int slice2,
// 			       const Array<int,1> &changedParticles,
// 			       int level)
// { 
// //   ifstream infile;
// //   infile.open("myPath");
// //   for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++)
// //     infile>>Path(slice,0)[0];
// //   infile.close();
// //   slice1=0;
// //   slice2=PathData.Path.NumTimeSlices()-2;
// //   level=0;
//   double levelTau=Path.tau* (1<<level);
//   double totalU=0.0;
//   int skip = 1<<level;
//   PathClass &Path=PathData.Path;
//   //  cerr<<"My total U starts at "<<totalU<<endl;
//   //  double J=0.1;
//   double J=5.0;
//   for (int slice=slice1;slice<slice2;slice+=skip){
//     //    totalU -= (J*cos(Path(slice,0)[0])*levelTau/2.0+J*cos(Path(slice+skip,0)[0])*levelTau/2.0);
//     totalU -= (0.1*cos(2*M_PI*Path(slice,0)[0])*levelTau/2.0+0.1*cos(2*M_PI*Path(slice+skip,0)[0])*levelTau/2.0);
//     //    totalU -= (2.0*(Path(slice,0)[0]*Path(slice,0)[0]/2.0)
//     //	       *levelTau/2.0+
//     //	       2.0*(Path(slice+skip,0)[0]*Path(slice+skip,0)[0]/2.0)
//     //	       *levelTau/2.0);
//     //    cerr<<"ACTION: "<<slice<<" "<<slice+skip<<" "<<totalU<<" "<<Path(slice,0)[0]<<" "<<Path(slice+skip,0)[0]<<endl;
    
//   }
  
//   double totalPotential=ComputeFourierAction();
//   //  cerr<<"Start paths"<<endl;
//   //  for (int slice=0;slice<PathData.Path.NumTimeSlices();slice++)
//   //    cerr<<Path(slice,0)[0]<<endl;
//   //  cerr<<"End paths "<<totalU<<" "<<totalPotential<<endl;
//   //  cerr<<"My totalU is "<<totalU<<endl;
//   //  cerr<<totalU<<" "<<totalPotential<<endl;
//   //  return totalU+totalPotential;
//   return totalPotential;
// }



double 
DualHermeleClass::d_dBeta (int slice1, int slice2, int level)
{
  return 0.0;
}


string
DualHermeleClass::GetName()
{
  return "DualHermele";
}
