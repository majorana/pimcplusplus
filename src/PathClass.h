#ifndef PATH_CLASS_H
#define PATH_CLASS_H

#include "Common/IO/InputOutput.h"
#include "MirroredClass.h"
#include "SpeciesClass.h"
#include "Common/Random/Random.h"
#include "CommunicatorClass.h"

///The number of time slices is the number of slices on this processor.
///In all cases this processor shares a time slice with the processor 
///ahead of it and behind it. The convention for the shared slices
///is that the processor owns its first but not its last slice.
class PathClass
{
private:  
  /// Path stores the position of all the particles at all time
  /// slices.  The order for access is timeslice, particle
  Mirrored2DClass<dVec> Path;
   /// Stores what species a particle belongs to
  Array<int,1> SpeciesNumber;
  Array<SpeciesClass *,1> SpeciesArray;
  int MyNumSlices;
  PIMCCommunicatorClass &Communicator;

  ////////////////////////////////
  /// Boundary conditions stuff //
  ////////////////////////////////
  dVec IsPeriodic;

  dVec Box, BoxInv;
  dVec kBox; //kBox(i)=2*Pi/Box(i)

  ///////////////////////////////////////////////
  /// k-space stuff for long-range potentials ///
  ///////////////////////////////////////////////
private:
  /// This is the maximum number of k vectors in each direction
  TinyVector<int,NDIM> MaxkIndex;
  /// True if we need k-space sums for long range potentials.
  bool LongRange;
  /// Stores the radius of the k-space sphere we sum over
  double kCutoff;
  /// Allocates and sets up the k-vectors.
  void SetupkVecs();
  /// Stores indices into C array for fast computation of rho_k
  Array<TinyVector<int,NDIM>,1> kIndices;
  /// This stores e^{i\vb_i \cdot r_i^\alpha}
  TinyVector<Array<complex<double>,1>,NDIM> C;
public:
  /// Stores the actual k vectors
  Array<dVec,1> kVecs;
  /// This holds the density of particles in k-space.  Indexed by
  /// (slice, species, k-vector).  Defined as
  /// \rho^\alpha_k = \sum_i e^{i\mathbf{k}\cdot\mathbf{r}_i^\alpha}
  /// Stores the kvectors needed for the reciporical space sum.
  /// Stores only half the vectors because of k/-k symmetry.
  Mirrored3DClass< complex<double> > Rho_k;
  void CalcRho_ks_Slow(int slice, int species);  
  void CalcRho_ks_Fast(int slice, int species);  

private:
  void ShiftPathData(int sliceToShift);
  void ShiftRho_kData(int sliceToShift);
public:
  Mirrored1DClass<int> Permutation;
  RandomClass Random;
  int TotalNumSlices;
  /// A scratch array to hold a boolean indicating whether we've
  /// looped over this particle yet
  Array<bool,1> DoPtcl;

  inline void  SetBox (dVec box);
  inline const dVec& GetBox();
  inline const double GetVol();
  inline void  SetPeriodic(TinyVector<bool,NDIM> period);

  /////////////////////////////////
  /// Displacements / Distances ///
  /////////////////////////////////
  inline void DistDisp (int slice, int ptcl1, int ptcl2,
			double &dist, dVec &disp);
  inline void DistDisp (int sliceA, int sliceB, int ptcl1, int ptcl2,
			double &distA, double &distB,dVec &dispA, dVec &dispB);
  inline double Distance (int slice, int ptcl1, int ptcl2);
  inline dVec Velocity (int sliceA, int sliceB, int ptcl);
  inline void PutInBox (dVec &v);


  //////////////////////////
  /// Data manipulations ///
  //////////////////////////
  inline const dVec& operator() (int slice, int ptcl) const;
  inline dVec& operator() (int slice, int ptcl);
  inline void SetPos (int slice, int ptcl, const dVec& r);
  inline int NumParticles();
  inline int NumTimeSlices();

  void MoveJoin(int oldJoin, int newJoin);      
  void ShiftData(int sliceToShift);
  void AcceptCopy(int startTimeSlice,int endTimeSlice, 
		  const Array <int,1> &activeParticle);
  void RejectCopy(int startTimeSlice,int endTimeSlice, 
		  const Array <int,1> &activeParticle );


  /////////////////////////////
  /// Species Manipulations ///
  /////////////////////////////
  inline int ParticleSpeciesNum(int ptcl);
  inline SpeciesClass& ParticleSpecies(int ptcl);
  inline SpeciesClass& Species(int speciesNum);
  inline int SpeciesNum (string name);
  inline void AddSpecies (SpeciesClass *newSpecies);
  inline int NumSpecies();


  //////////////////////////
  /// IO and allocations ///
  //////////////////////////
  void Read(IOSectionClass &inSection);
  void Allocate();

  inline PathClass(PIMCCommunicatorClass &communicator);
  friend void SetupPathNaCl(PathClass &path);
  friend void SetupPathZincBlend(PathClass &path);
};




inline int PathClass::SpeciesNum (string name)
{
  int i=0;
  while ((i<SpeciesArray.size()) && (SpeciesArray(i)->Name != name))
    i++;
  if (i == SpeciesArray.size())
    return -1;
  else
    return i;
}


/// Return what species type a particle belongs to;
inline int PathClass::ParticleSpeciesNum(int ptcl)
{
  return (SpeciesNumber(ptcl));
}
/// Return a reference to the species that particle ptcl belongs to
inline SpeciesClass& PathClass::ParticleSpecies(int ptcl) 
{
  return *(SpeciesArray(SpeciesNumber(ptcl)));
}
/// Return a species references
inline SpeciesClass& PathClass::Species(int speciesNum)
{
  return (*(SpeciesArray(speciesNum)));
}


/// Returns the number of particle Species
inline int PathClass::NumSpecies() 
{ 
  return SpeciesArray.size();
}

inline int PathClass::NumParticles() 
{ 
  return Path.cols();
}

///The number of time slices is the number of slices on this processor.
///In all cases this processor shares a time slice with the processor 
///ahead of it and behind it. The convention for the shared slices
///is that the processor owns its first but not its last slice.
inline int PathClass::NumTimeSlices() 
{ 
  return Path.rows();
}


/// Returns the position of particle ptcl at time slice timeSlice
inline const dVec& PathClass::operator() (int slice, int ptcl) const
{ 
  return Path(slice, ptcl); 
}

/// Returns the position of particle ptcl at time slice timeSlice
inline dVec& PathClass::operator() (int slice, int ptcl)
{ 
  return Path(slice, ptcl); 
}

inline void PathClass::SetPos(int slice, int ptcl, const dVec &r)
{
  (*this)(slice, ptcl) = r;
}

inline void PathClass::AddSpecies (SpeciesClass *newSpecies)
{
  int numSpecies = SpeciesArray.size();
  /// Add an element for the new species
  SpeciesArray.resizeAndPreserve(numSpecies+1);
  SpeciesArray(numSpecies) = newSpecies;
}

inline PathClass::PathClass (PIMCCommunicatorClass &communicator) : 
  Communicator(communicator), Random(Communicator)
{
  //      NumSpecies = 0;
  TotalNumSlices=0;
  Random.Init(314159);
}


////////////////////////////////
/// Boundary conditions stuff //
////////////////////////////////

inline void PathClass::SetBox (dVec box)
{
  Box = box;
  for (int i=0; i<NDIM; i++){
    BoxInv(i) = 1.0/box(i);
    kBox(i)=2*M_PI*BoxInv(i);
  }
  
}

inline const dVec& PathClass::GetBox()
{
  return Box;
}

inline const double PathClass::GetVol()
{
  double  vol=1.0;
  for (int i=0;i<NDIM;i++){
    vol*=Box(i);
  }
  return vol;
}


inline void PathClass::SetPeriodic(TinyVector<bool,NDIM> period)
{
  for (int i=0; i<NDIM; i++)
    IsPeriodic(i) = period(i) ? 1.0 : 0.0;
}

inline void PathClass::DistDisp (int slice, int ptcl1, int ptcl2,
			       double &dist, dVec &disp)
{
  disp = Path(slice, ptcl2) -Path(slice, ptcl1);
  
  for (int i=0; i<NDIM; i++) {
    double n = -floor(disp(i)*BoxInv(i)+0.5);
    disp(i) += n*IsPeriodic(i)*Box(i);
  }
  dist = sqrt(dot(disp,disp));

#ifdef DEBUG
  dVec DBdisp = Path(slice, ptcl2) -Path(slice, ptcl1);
  for (int i=0; i<NDIM; i++) {
    while (DBdisp(i) > 0.5*Box(i))
      DBdisp(i) -= Box(i);
    while (DBdisp(i) < -0.5*Box(i))
      DBdisp(i) += Box(i);
    assert (fabs(DBdisp(i)-disp(i)) < 1.0e-12);
  }
#endif
}


inline void PathClass::DistDisp (int sliceA, int sliceB, int ptcl1, int ptcl2,
			       double &distA, double &distB, 
			       dVec &dispA, dVec &dispB)
{  
  dispA = Path(sliceA, ptcl2) - Path(sliceA,ptcl1);
  dispB = Path(sliceB, ptcl2) - Path(sliceB,ptcl1);
  //  dVec tempDispB;
  //  dVec tempDispBN;
  //  dVec dispBNew;
//   cerr << "A1 = " << Path(sliceA,ptcl1) << endl;
//   cerr << "A2 = " << Path(sliceA,ptcl2) << endl;
//   cerr << "B1 = " << Path(sliceB,ptcl1) << endl;
//   cerr << "B2 = " << Path(sliceB,ptcl2) << endl;
  int m;
  for (int i=0; i<NDIM; i++) {
    double n = -floor(dispA(i)*BoxInv(i)+0.5);
    dispA(i) += n*IsPeriodic(i)*Box(i);
//     double m = -floor(dispB(i)*BoxInv(i)+0.5);
//     dispB(i) += m*IsPeriodic(i)*Box(i);
    double mNew=-floor((dispA(i)-dispB(i))*BoxInv(i)+0.5);
    dispB(i)-= mNew*IsPeriodic(i)*Box(i);
//     // HACK HACK HACK
//     m=0;
//     tempDispB(i) = dispB(i)+m*IsPeriodic(i)*Box(i);
//     tempDispBN(i)=dispB(i)-m*IsPeriodic(i)*Box(i);
//     while (fabs(dispA(i)-tempDispB(i))>Box(i)/2.0 &&
// 	   fabs(dispA(i)-tempDispBN(i))>Box(i)/2.0){
//       m=m+1;
//       tempDispB(i) = dispB(i)+m*IsPeriodic(i)*Box(i);
//       tempDispBN(i)=dispB(i)-m*IsPeriodic(i)*Box(i);
//     }
//     if (fabs(dispA(i)-tempDispB(i))<=Box(i)/2.0)
//       dispB(i)=tempDispB(i);
//     else if (fabs(dispA(i)-tempDispBN(i))<=Box(i)/2.0)
//       dispB(i)=tempDispBN(i);
//     else cerr<<"ERROR! ERROR! ERROR!"<<endl;
//     if (fabs(dispBNew(i)-dispB(i))>1e-12){
//       cerr<<"dispBNew and dispB are not the same!\n";
//     }
//     //    double m = -floor(dispB(i)*BoxInv(i)+0.5);
//     //    dispB(i) += m*IsPeriodic(i)*Box(i);
//     //    cerr << "n = " << n << endl;
  }
//   cerr << "dispA = " << dispA << endl;
//   cerr << "dispB = " << dispB << endl;
  distA = sqrt(dot(dispA,dispA));
  distB = sqrt(dot(dispB,dispB));

#ifdef GARBAGEDEBUG
  dVec DBdispA = Path(sliceA, ptcl2) -Path(sliceA, ptcl1);
  dVec DBdispB = Path(sliceB, ptcl2) -Path(sliceB, ptcl1);
  for (int i=0; i<NDIM; i++) {
    while (DBdispA(i) > 0.5*Box(i)) 
      DBdispA(i) -= Box(i);
    while (DBdispA(i) < -0.5*Box(i)) 
      DBdispA(i) += Box(i);
    while ((DBdispB(i)-DBdispA(i)) > 0.5*Box(i))
      DBdispB -= Box(i);
    while ((DBdispB(i)-DBdispA(i)) < -0.5*Box(i))
      DBdispB += Box(i);
  }
//   cerr << "DBdispA = " << DBdispA << endl;
//   cerr << "DBdispB = " << DBdispB << endl;
  for (int i=0; i<NDIM; i++) {
    assert (fabs(DBdispA(i)-dispA(i)) < 1.0e-12);
    assert (fabs(DBdispB(i)-dispB(i)) < 1.0e-12);
  }
#endif
}

inline dVec PathClass::Velocity (int sliceA, int sliceB, int ptcl)
{
  dVec vel = Path(sliceB, ptcl) - Path(sliceA,ptcl);
  PutInBox(vel);


#ifdef DEBUG
 dVec DBvel = Path(sliceB,ptcl) - Path(sliceA,ptcl);
   for (int dim=0; dim<NDIM; dim++)
     {
       while (DBvel[dim] > (0.5*Box[dim]))
  	DBvel[dim] -= Box[dim];
       while (DBvel[dim] < -(0.5*Box[dim]))
	 DBvel[dim] += Box[dim];
       assert(fabs(DBvel[dim]-vel[dim])<1e-12);
     }

#endif

  return vel;
}


inline void PathClass::PutInBox (dVec &v)
{
#ifdef DEBUG
  dVec Dv=v;
#endif
  for (int i=0; i<NDIM; i++) {
    double n = -floor(v(i)*BoxInv(i)+0.5);
    v(i) += n*IsPeriodic(i)*Box(i);
  }
#ifdef DEBUG
  for (int dim=0; dim<NDIM; dim++){
    while (Dv[dim] > (0.5*Box[dim]))
      Dv[dim] -= Box[dim];
    while (Dv[dim] < (-(0.5*Box[dim])))
      Dv[dim] += Box[dim];
    assert(fabs(Dv[dim]-v[dim])<1e-12);
  }
#endif 
  
}

#endif
