#ifndef DISTANCE_TABLE_CLASS_H
#define DISTANCE_TABLE_CLASS_H


class DistanceTableClass
{
 private:
  SpeciesArrayClass &SpeciesArray;
  MirroredArrayClass<double> DistTable; /// (Global ptclXptcl number, timeslice)
  MirroredArrayClass<dVec> DispTable; /// (Global ptclXptcl number, timeslice)
  inline int ArrayIndex(ParticleID particle1,ParticleID particle2);
  inline int ArrayIndex(int particle1, int  particle2);
 public:
  DistanceTableClass(SpeciesArrayClass p_speciesArray) : SpeciesArray(p_speciesArray){}
  
  inline void GetDistDisp(ParticleID particle1, ParticleID particle2, 
			  int timeSlice, double &dist, dVec &disp);  

  inline void Update(ParticleID particle, int timeSlice);


}

inline void DistanceTable::Update(ParticleID particle, int timeSlice)
{
  ParticleID particle2;
  particle2[1]=0;
  for (int speciesCounter=0;speciesCounter<SpeciesArray.Size();speciesCounter++){
    particle2[0]=speciesCounter;
    int globalNum2=ArrayIndex(particle2);
    for (int particleCounter=0;particleCounter<SpeciesArray(speciesCounter).size();particleCounter++){
      
      globalNum2++;  
    }

  }

}

inline int DistanceTable::ArrayIndex(int global1, int global2)
{
  int totalNum= ((global1*(global1+1))>>1)+global2;   //global1*(global1+1)/2+global2;  
  return totalNum;
}

inline int DistanceTable::ArrayIndex(ParticleID particle1,ParticleID particle2)
{
  int global1= SpeciesArray.ParticleID2IntMap(particle1);
  int global2 = SpeciesArray.ParticleID2IntMap(particle2);
  int totalNum=ArrayIndex(global1,global2);

  return totalNum;

}

inline void GetDistDisp(ParticleID particle1, ParticleID particle2, 
			int timeSlice, double &dist, dVec &disp)
{
  int globalNum=getGlobalNumberMap(particle1,particle2);
  dist=DistTable(globalNum,timeSlice);
  disp=DispTable(globalNum,timeSlice);
}


#endif 
