#ifndef DISTANCE_TABLE_CLASS_H
#define DISTANCE_TABLE_CLASS_H

/// This class stores the distances and displacements between particles 
/// of the same time slice
class DistanceTableClass
{
 private:
#ifdef PBC
  dVec Box;
#endif
  ///Pointer to the SpeciesArray
  SpeciesArrayClass &SpeciesArray;
  ///Table of distances.  This is stored in a two-d array where the first 
  ///dimension of the array is indexed by an integer that you get by 
  ///mapping two particleID's to it through our ArrayIndex Function.
  ///(The second dimensions is timeslices).
  MirroredArrayClass<double> DistTable; ///<(Global ptclXptcl number,timeslice)
  ///Table of displacements. Set up in same way as distTable.
  MirroredArrayClass<dVec> DispTable; ///< (Global ptclXptcl number, timeslice)

  ///Table of Image numbers (stored 0-26)
  MirroredArrayClass<int> ImageNumTable; ///<(Global ptclXptclnumber,timeslice)
  Array<dVec,1> ImageVectors;
  ///Maps particleID x particleID to a location in our distance or 
  ///displacement table.
  inline int ArrayIndex(ParticleID particle1,ParticleID particle2);
  ///Maps particleInt x particleInt to a location in our distance 
  ///or displacement table.
  ///particleInt gotten from SpeciesArray.ParticleID2Int
  inline int ArrayIndex(int particle1, int  particle2);
 public:
  ///Constructor. Initializes the speciesArray reference
  DistanceTableClass(SpeciesArrayClass p_speciesArray) : SpeciesArray(p_speciesArray){}
  ///Returns the distance and displacement.
  inline void GetDistDisp(ParticleID particle1, ParticleID particle2,
			  int timeSlice,double &dist, dVec &disp);
  ///Returns the distance and displacement. This is for the case where
  ///you are calculating the disp/dist and want to make sure you have the 
  ///same image for a pair of particles.
  inline void GetDistDisp(ParticleID particle1, ParticleID particle2, 
			  int timeSliceA, int timeSliceB,
			  double &distA, dVec &dispA,
			  double &distB, dVec &dispB);  
  ///Updates all of the distances between "particle1 and all the other particles
  inline void Update(ParticleID particle1, int timeSlice);

  void SetBox(dVec myBox);
  DistanceTableClass(SpeciesArrayClass &mySpeciesArray, dVec myBox) :
    SpeciesArray(MySpeciesArray), Box(myBox)
  {
    int NumVectors = 1;
    for (int i=0; i<NDIM; i++)
      NumVectors *= 3;
    ImageVectors.resize (NumVectors);
    
    if (NDIM == 2)
      {
	for (int i=-1; i<=1; i++)
	  for (int j=-1; j<=1; j++)
	    {
	      int num = (i+1)*3+(j+1);
	      ImageVector(num)[0] = (double)i*Box[0];
	      ImageVector(num)[1] = (double)j*Box[1];
	    }
      }
    if (NDIM == 3)
      {
	for (int i=-1; i<=1; i++)
	  for (int j=-1; j<=1; j++)
	    for (int k=-1; k<=1; k++)
	    {
	      int num = (i+1)*9+(j+1)*3+(k+1);
	      ImageVector(num)[0] = (double)i*Box[0];
	      ImageVector(num)[1] = (double)j*Box[1];
	      ImageVector(num)[2] = (double)k*Box[2];
	    }
      }
  }

};

///Updates all of the distances between "particle1 and all the other particles
inline void DistanceTable::Update(ParticleID particle1, int timeSlice)
{
  particleID particle2;
  int particle1Int=SpeciesArray.particleID2Int(particle1);
  particle2Int=0;
  for(int speciesCounter=0;speciesCounter<SpeciesArray.Size();speciesCounter++)
    {
      particle2[0]=speciesCounter;
      for (int particleCounter=0;
	   particleCounter<SpeciesArray(speciesCounter).Size();
	   particleCounter++){
	particle2[1]=particleCounter;
	int distTableIndex=ArrayIndex(particle1Int,particle2Int);
	/// Particle 1 is always in the box, 
	/// ImageNum is box where particle 2 is in 
	SpeciesArray.DistDisp(particle1,particle2,timeSlice,timeSlice,
			      DistTable(distTableIndex,timeSlice),
			      DispTable(distTableIndex,timeSlice),
			      ImageNumTable(distTableIndex,timeSlice)); 
	
	/// This assumes that our mapping for particleID -> particleInt
	/// is monotonic in particle 
	/// number and species, so we just iterate particle2Int
	/// instead of asking for the mapping again.
	particle2Int++;
      }
    }
}

    



/// Maps particleInt x particleInt to a location in our distance or
/// displacement table. particleInt gotten from
/// SpeciesArray.ParticleID2Int 
inline int DistanceTable::ArrayIndex(int global1, int global2)
{
  //global1*(global1+1)/2+global2;  
  int totalNum= ((global1*(global1+1))>>1)+global2;   
  return totalNum;
}

/// Maps particleID x particleID to a location in our distance or
/// displacement table. 
inline int DistanceTable::ArrayIndex(ParticleID particle1,ParticleID particle2)
{
  int global1= SpeciesArray.ParticleID2Int(particle1);
  int global2 = SpeciesArray.ParticleID2Int(particle2);
  int totalNum=ArrayIndex(global1,global2);

  return totalNum;
}


/// Returns the distance and displacement.
inline void GetDistDisp(ParticleID particle1, ParticleID particle2, 
			int timeSlice, double &dist, dVec &disp)
{
  int distanceTableIndex=ArrayIndex(particle1,particle2);
  dist=DistTable(distanceTableIndex,timeSlice);
  disp=DispTable(distanceTableIndex,timeSlice);
}

/// Returns the distance and displacement. This is for the case where
/// you are calculating the disp/dist and want to make sure you have the 
/// same image for a pair of particles.
inline void GetDistDisp(ParticleID particle1, ParticleID particle2, 
			int timeSliceA, int timeSliceB,
			double &distA, dVec &dispA,
			double &distB, dVec &dispB){
  int distanceTableIndex=ArrayIndex(particle1,particle2);
 /// We're making an implicit assumption around here that the number of
 /// time slices per particle is the same  
  dispA=DispTable(distanceTableIndex,timeSliceA);
  distA=DistTable(distanceTableIndex,timeSliceA);
  if (ImageNumTable(distanceTableIndex,timeSliceB) !=
      ImageNumTable(distanceTableIndex,timeSliceA))
    {
      dispB = SpeciesArray(particle1)


  dispB=DispTable(distanceTableIndex,timeSliceB)+
    (ImageNumTable(distanceTableIndex,timeSliceA)-
     ImageNumTable(distanceTableIndex,timeSliceB))*
    
    SpeciesArray.DistDisp(particle1,particle2,timeSlice,timeSlice,
			  DistTable(distTableIndex,timeSlice),DispTable(distTableIndex,timeSlice),
			  ImageNumTable(distTableIndex,timeSlice)); ///< Particle 1 is always in the box, 
			  ///<imageNum is box where particle 2 is in 



#endif 
