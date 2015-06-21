
```
///Langevin Dynamics example
double tau=0.125;

Section (Parallel)
{
  int ProcsPerClone = 1;
}

Section (System)
{
  int NumTimeSlices=200;
  Array<bool,1> IsPeriodic(3) = [true,true,true];
  Array<double,1> Box(3)=[20.0,20.0,20.0];
  // double kCutoff = 0.5;
  Section (Particles)
  {
    Section (Species)
    {
      string Name="e";
      string Type="e";
      double lambda=0.5;
      string Statistics="BOLTZMANNON";
      int NumParticles=2;
      int NumDim=3;
      string InitPaths="LEVIFLIGHT";
      Array<double,2> Positions(2,3)=[-0.5, 0.0, 0.0,
				       0.5, 0.0, 0.0 ];
    }
    Section (Species)
    {
      string Name="p";
      string Type="p";
      double lambda=0.0; //classical particle
      string Statistics="BOLTZMANNON";
      int NumParticles=2;
      int NumDim=3;
      string InitPaths="FIXED";
      Array<double,2> Positions(2,3)=[-0.7, 0.0, 0.0, 
				      0.7, 0.0, 0.0 ];
    }
  }
}

Section (Action)
{
  int NumImages = 0;
  int MaxLevels = 5;
  Array<string,1> PairActionFiles(3) = ["../e-p_beta8.0.coulombBCfit.h5",
                                        "../e-e_beta8.0.coulombBCfit.h5",
                                        "../p-p.PairAction"];
  int NumBreakupKnots = 20;
  bool UseLongRange = false;
  bool UseBackground = false;
   Section (Action)
   {
       string Name="ShortRange";
       string Type="ShortRange";
   }

}

Section (Observables)
{
  string OutFileBase = "H2";
  Section (Observable)
  {
    string Type = "Energy";
    string Name = "Energy";
    string Description="Total Energy";
    int Frequency = 1;
  }
  Section (Observable) {
    string Type = "Forces";
    string Name = "Forces";
    string Species = "p";
    int Frequency = 1;
  }
}

Section (Moves) {
  Section (Move) {
    string Type="ShiftMove";
    string Name="Shift";
  }
  Section (Move) {
    string Type="BisectionBlock";
    Array<string,1> HigherLevelActions(1)=["ShortRange"];
    Array<string,1> SamplingActions(1)=["ShortRange"];

    string Name="electronBlock";
    string PermuteType="NONE";
    int NumLevels=5;
    string Species = "e";
    int StepsPerBlock=3;    
  }
}



Section (Algorithm)
{
   Section (Loop){ // Before Accumulating Observables
     int Steps=1000; 
     Section (Loop){ 
       int Steps=3;
       Section (Move)  { string Name="electronBlock";   }
     }
     Section (Move)    { string Name = "Shift";         }
   }
  Section (Loop){ // Start Accumulating observables
    int Steps=3000000;
    Section (Loop){ //Write Loop
      int Steps=500;
      Section (Loop){
	int Steps=7;
	Section (Move)  { string Name = "electronBlock"; } 
	Section (Move)  { string Name = "Shift";         }
      }

      Section (Observe) { string Name = "Energy";        }
      Section (Observe) { string Name = "Forces";        }
      //Section (Observe) {string Name = "PathDump";}  //As a good exercise, add this back! 
    }
    Section (WriteData){}
  }
}
```