#ifdef BUILD_DEV
  #include "PathClassDev.h"
#else
  #include "PathClass.h"
#endif

// void 
// PathClass::InitializeJosephsonCode()
// {
// //   int N=NumTimeSlices();
// //   inPhi    = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N-1));
// //   outOmega = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N-1));
// //   phi2omegaPlan=fftw_plan_dft_1d(NumTimeSlices()-1,inPhi,outOmega,FFTW_BACKWARD,FFTW_MEASURE);

// //   inOmega    = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N-1));
// //   outPhi = (fftw_complex*)fftw_malloc(sizeof(fftw_complex)*(N-1));
// //   omega2phiPlan=fftw_plan_dft_1d(NumTimeSlices()-1,inOmega,outPhi,FFTW_FORWARD,FFTW_MEASURE);
// //   fftw_destroy_plan(phi2omegaPlan);
// //   fftw_destroy_plan(omega2phiPlan);
// }


// void 
// PathClass::Phi2Omega()
// {

// //   int N=NumTimeSlices()-1;
// //   for (int slice=0;slice<N;slice++){
// //     inPhi[slice][0]=Path(slice,0)[0];
// //     inPhi[slice][1]=0.0;
// //   }
// //   fftw_execute(phi2omegaPlan);

// //   for (int counter=0;counter<N;counter++){
// //     outOmega[counter][0]=outOmega[counter][0]*tau;
// //     outOmega[counter][1]=outOmega[counter][1]*tau;
// //   }

// }


// void
// PathClass::Omega2Phi()
// {
// //   double T=1.0/(TotalNumSlices*tau);
// //   int N=NumTimeSlices()-1;
// //   for (int slice=0;slice<N;slice++){
// //     inOmega[slice][0]=outOmega[slice][0];   //Path(slice,0)[0];
// //     inOmega[slice][1]=outOmega[slice][1];   //0.0;
// //   }
// //   fftw_execute(omega2phiPlan);
// //   for (int counter=0;counter<N;counter++){
// //     outPhi[counter][0]=outPhi[counter][0]*T;
// //     outPhi[counter][1]=outPhi[counter][1]*T;
// //   }
// }


