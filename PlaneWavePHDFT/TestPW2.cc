#include "PlaneWaves.h"

main()
{
  //  Vec3 box(5.32117439923148, 5.32117439923148, 5.32117439923148);
  Vec3 box(10.0, 11.0, 12.0);
  //  Vec3 box(6.5, 6.500, 6.501);
  double kcut = 3.0;
  IOSectionClass in;
  in.OpenFile ("NaUnscreenedPH_Feb18_05.h5");
  //in.OpenFile ("NaPH_US_March1_05b.h5");
  Potential *ph = ReadPotential(in);
  in.CloseFile();

  SystemClass system(10);
  Vec3 k(0.0, 0.0, 0.0);
  //system.Setup (box, k, kcut, *ph, true); 
  system.Setup (box, k, kcut, 1.0, true);
  Array<Vec3,1> rions(16);
  rions(0)  = Vec3 ( -1.017425796249996E-005, -2.030943049556812E-005,  4.920717706347283E-005);
  rions(1)  = Vec3 ( -0.500075239272832     , 3.968267796505207E-005 , -6.913079201222320E-005);
  rions(2)  = Vec3 (  1.669465221995777E-004,-0.500018171330705      , -9.856628047657987E-005);
  rions(3)  = Vec3 ( -0.499851202317107     ,-0.500096683844841      , -1.700448925197800E-004);
  rions(4)  = Vec3 (  1.361041380271397E-005,-1.369924089405148E-004 , -0.500063325084747 );
  rions(5)  = Vec3 ( -0.500119521962410     ,-1.829550258971563E-004 , -0.500123087616722);
  rions(6)  = Vec3 ( -3.677051036856182E-005,-0.499819600220644      , -0.499876068418090);
  rions(7)  = Vec3 ( -0.499902317553432     ,-0.499841168043744      , -0.500185624696847);
  rions(8)  = Vec3 (  0.250010721736608     , 0.249897889290116      ,  0.249837400918737);
  rions(9)  = Vec3 ( -0.250092738705326     , 0.249852899399392      ,  0.250141073709984);
  rions(10) = Vec3 (  0.249906141622914     ,-0.249978185426867      ,  0.249827459571552);
  rions(11) = Vec3 ( -0.250088430489722     ,-0.249918688645085      ,  0.249820301147556);
  rions(12) = Vec3 (  0.250078917520793     , 0.250124992394034      , -0.250149452551939);
  rions(13) = Vec3 ( -0.250030715721095     , 0.249866475422417      , -0.250163701202200);
  rions(14) = Vec3 (  0.249903446595725     ,-0.250033815724168      , -0.250138869986818);
  rions(15) = Vec3 ( -0.249894729888688     ,-0.249903398844062      , -0.249977204406117);
  for (int i=0; i<16; i++)
    for (int j=0; j<3; j++) {
      rions(i)[j] = box[j] * rions(i)[j];
      //rions(i)[j] = box[j]*(drand48()-0.5);
    }
  
  system.SetIons (rions);
  system.DiagonalizeH();
  

//   FILE *fout = fopen ("HRho_k10.0.dat", "w");
//   Vec3 r(0.0, 0.0, 0.0);
//   for (double x=-0.5*box[0]; x<=0.5*box[0]; x+=0.01) {
//     r[0] = x;
//     fprintf (fout, "%1.12e ", x);
//     for (int band=0; band<numBands; band++) {
//       complex<double> psi(0.0, 0.0);
//       for (int i=0; i<H.GVecs.size(); i++) {
// 	double phase = dot (H.GVecs(i), r);
// 	double s,c;
// 	sincos (phase,&s, &c);
// 	complex<double> z(c, s);
// 	psi += z * CG.Bands(band,i);
//       }
    
//       psi /= sqrt(H.GVecs.GetBoxVol());
//       fprintf (fout, "%1.12e ", real(conj(psi)*psi));
//     }
//     fprintf (fout, "\n");
//   }
//   fclose (fout);
  

//   fprintf (stderr, "Time = %1.3f\n", 
// 	   (double)(end-start)/(double)CLOCKS_PER_SEC);

}
