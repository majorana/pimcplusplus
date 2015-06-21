In this tutorial, we will first calculate the energy of a single particle in a box with periodic boundary conditions. This calculation will introduce you to the PIMC++ software suite, including its input format, output analysis and processing tools.

In this first exercise, you will be exposed to a great deal of terminology and syntax related to path integrals and used in PIMC++. We will alert you to important points; occasionally, links embedded in the text will take you to more detailed information about technical or expert-level topics. However, many aspects of the simulation input and software design will be explored in more detail later, so do not worry about understanding every aspect right away.

Specifically, we will calculate the energy of a single particle in a box of size 10 angstrom. at a temperature of 1.666K.  For consistency with the Helium simulations (in the next tutorial) we'll do, we're working in a system of units in which energy is in Kelvin, length is in Angstroms and $\hbar=1$ In these units, the free particle mass will be $0.5 K<sup>{-1} \AA</sup>{-2}, which corresponds to $4.034 \times 10^{-26}$ kg.
In this simple case, it will be possible to compare the simulation with the [analytical result](ParticleInABox.md). If you do this, you should get an energy of 2.51 K.

## The Basic Decisions ##

Let's outline the basic decisions you need to make when setting up a path integral simulation:
  1. The first thing you need to decide on is your system:  How many particles? What's the box size? Are you running with distinguishable particles (called Boltzmannon's) or Bosons or Fermions?
  1. What's your Hamiltonian?  PIMC is essentially sampling the integral $<R | \exp(-\beta H) |R>$.  So you need to know what H is.  In the code, this means you need to pick the pieces of H which means you need to select your actions.  Typical actions might be the kinetic action and a potential action (say a particle-particle interaction)
  1. What moves are you using to sample your system? For example, will you use a Bisection Move or a Displace Move (only works with distinguishable particles).  It is critical that you always use a ShiftMove.  When you pick out the moves, you have to tell them which actions to use.
  1. What observables are you using? In other words, what is it you want to measure?
  1. Put it all together:  You need to build up an algorithm. When do you want to use a certain move. How often do you want to write out the observables you are making/

So for the single particle in a box you might answer these questions like:

Single distinguishable particle in a box of size ... with a hamiltonian that only include the kinetic action, simple bisection moves, and measuring the energy.

## The first input ##

Once you have made these decisions, you need to put them together in an input file.
We will first describe how the input file is organized and how important aspects of the simulation are controlled. The first input file looks like this:
You can view a fully annotated version of this file here in order to see more detail about input syntax and options.

```
double tau=.025;
Section (Parallel)
{
  int ProcsPerClone = 1;
}
Section (System)
{
    int NumTimeSlices=12;
    Array<double,1> Box(3)=[10.0,10.0,10.0];
    Array<bool, 1 > IsPeriodic(3)=[true,true,true];
    Section (Particles)
    {
        Section (Species)
        {
          string Name="free";
          string Type="free";
          double lambda=1.0;
          string Statistics="BOLTZMANNON";
          int NumParticles=1;
          int NumDim=3;
          string InitPaths="BCC";
        }
   }
}
Section (Action)
{
   int NumImages=0;
   int MaxLevels = 2;
   Array<string,1> PairActionFiles(1) = ["zero.PairAction"];
   Section (Action)
   {
       string Name="ShortRange";
       string Type="ShortRange";
   }

}

Section (Observables)
{
    string OutFileBase = "SingleParticle";
    Section (Observable)
    {
        string Type = "Energy";
        string Name = "Energy";
        string Description="Total Energy";
        int Frequency=1;
    }
}
Section (Moves){
    Section (Move) {
        string Type="BisectionBlock";
        string Name="BisectionBlock";
        Array<string,1> HigherLevelActions(1)=["ShortRange"];
        Array<string,1> SamplingActions(1)=["ShortRange"];
        string PermuteType="NONE";
        string Species="free";
        int NumLevels=2;
        int StepsPerBlock=2;
    }
    Section (Move)
    {
        string Type="ShiftMove";
        string Name="Shift";
    }
}

Section (Algorithm)
{
    Section (Loop){
        int Steps=1000;
     Section (Loop){
         int Steps=50;
         Section (Move) {string Name="BisectionBlock";}
         Section (Observe) {string Name = "Energy"; }
         Section (Move) {string Name = "Shift"; }
     }
     Section (WriteData){}
   }
}
```

**Helpful editing tip for input files:** The input files are inspired by C++ syntax. If you turn on C++ syntax highlighting in your editor (Meta-x c++-mode followed by Meta-x font-lock-mode in emacs), they will be colored and the braces will be matched in a friendly way.

**A comment about additional files you need:** Although this is the main file you need for input, you will also need files to tell you what the interaction is between the particles (in this case, the interactions is zero but we have to tell the code this as well!. This includes two files:
  1. PairAction File:  This file tells you some basic information and gives you the path to the h5 file (below).  It is important to get this path right!  In this case, take http://www.physics.princeton.edu/~bclark/Tutorials/pimc++_ff/zero.PairAction
  1. h5 File:  This has all the data in it.  We will learn later how to create these but for the moment it's just important that we have it.  In this case, we don't need this file because there really isn't any PairAction.

Let's briefly examine the input file. Notice the content is grouped by Section (that largely corresponds to the decisions you need to make to set up your simulation), each of which contains other hierarchal sections or variables that define the input to the simulation.

In summary, the sections are :
  * Parallel: Specifies how many processors over which each path is distributed.
  * System: Information about the physical system to be simulated. Contains number of time slices, tau, box-size, and the particle Section.
    * Particle: Contains a section for each Species.
    * Species Section: Contains lambda (), statistics, number of particles, etc. A Species can be an element, ion, molecule, or any object with uniform physical properties.
  * Moves: Defines the moves of your Monte Carlo
  * Observables: Defines what properties of the system you are measuring
  * Action: The action (dimensionless energy) is used to evaluate whether to accept or reject a Monte Carlo move. By specifying which actions to use and their properties, the description of the system's physical interactions is fixed (i.e. free particles, Coulomb, Lennard-Jones, etc.).
  * Algorithm: Specifies the simulation sequence of events (when each Move and Observable is called)

**A note about units:** PIMC++ has no intrinsic units and any system of units can be used by specifying consistent values in the input file.  , where  is the number of time slices.   has units of inverse energy and establishes the temperature scale.

For the tutorial, input files will use units of Kelvin and Angstrom with ![http://chart.googleapis.com/chart?cht=tx&chl=\hbar=1&nonsense=something_that_ends_with.png](http://chart.googleapis.com/chart?cht=tx&chl=\hbar=1&nonsense=something_that_ends_with.png)

Be sure to understand which lines of the input file establish that you are simulating a particle of mass, temperature  and box.  The mass is defined through lambda (![http://chart.googleapis.com/chart?cht=tx&chl=\lambda=\hbar^2/2m=1&nonsense=something_that_ends_with.png](http://chart.googleapis.com/chart?cht=tx&chl=\lambda=\hbar^2/2m=1&nonsense=something_that_ends_with.png)) (Section: System/Species/Particle). The temperature is defined implicitly as
$1/(NumTimeSlices \times tau)$ with NumTimeSlices (Section:System) and tau (Section: System) specified in the input.

Let's run PIMC++ with the input file above. PIMC++ can be run by typing into the terminal,
`pimc++ SingleParticle.in`

## Analysis ##

You will see that the code has produced a file called SingleParticle.0.h5.  Where does this name come from?  The "SingleParticle" piece comes because "OutFileBase" (Section: Observables) was set to be "SingleParticle."  Each parallel clone that PIMC++ then gets a number appended to it (starting with 0).  Since you ran the code in serial mode, you only get one output file and it has a 0 appended. Finally the ".h5" indicates the output files are written in a portable, hierarchical file format known as HDF5.

Let us take a look at the output we've generated. In order to do this, we will run the analysis script Report.py on the output. On the lab machines, the script will be in your path and you can just enter

`Report.py SingleParticle`

If you have a different setup where the top-level PIMC++ installation directory is PIMC++, the analysis script can be found at PIMC++/src/analysis/Report.py.
This script will produce an HTML file that will give us some information about our run. Examine SingleParticle/index.html with a web browser. On the lab machines, we recommend using

`firefox SingleParticle/index.html &`

Your screen should be similar to the following:
![http://cms.mcc.uiuc.edu/pimcpp/images/0/07/PimcppAnalysis.jpg](http://cms.mcc.uiuc.edu/pimcpp/images/0/07/PimcppAnalysis.jpg)

and the html that is created is seen  as [index.html](http://www.physics.princeton.edu/~bclark/Tutorials/SingleParticle/index.html).


Let us take note of the following information in the html summary.
  * Run information: (Upper left table) Specifies where and when it was run, who ran it, etc.
  * System Data: (Upper right table) System Parameters including temperature, species information, etc.
  * Move information: Table specifying acceptance ratio of moves and center of mass drift (add this to your input). You will notice that in this simulation, the acceptance ratio is 1.00. Why is this?. There is additional information on the bisection move here. You can ignore this for now. Note: There is another "move" in the simulation, "ShiftMove", that isn't in the table. It is not a real move, but simply rearranges the storage of the path in memory. (It is the same as "RelabelBeads" from yesterday's Python tutorial.) Thus, it has no acceptance ratio, but it is nonethless required to ensure ergodicity. Go here to understand the technical reason behind this.
  * Observable Information: We have measured the Energy observable in this simulation. Notice the energy table with the different components of the energy and their respective variance, error and autocorrelation time. We see that the energy is close to  as anticipated above. Of course, as we lower the temperature, we would see quantum effects becoming more important as the energy approaches the ground state for a single particle in a box. (Note: there is a subtlety to getting the correct energy if your box is small compared to  in PIMC++ caused by an approximation to the kinetic action. This can be skipped on the first pass through the tutorial, but experts can explore this NumImage approximation here ) Although a summary of the energy is convenient, it is often useful to be able to see a trace plot of the energy (for example to choose an equilibration time). We can get such detailed information by clicking on the  "Detailed HTML Page" at the top of index.html. Do that now. This page also allows us to see that the data was output 1000 times. By examining the Algorithm section of our input file, we see that the data gets written every time a WriteData is encountered in the algorithm (nb: this is NOT true for the PathDump observable which is special!). Looking at the number of iterations given for each of the Loop sections, we see that WriteData will be called 1000 times.

At some point, you are going to have to get to the data yourself and do some analysis of.  Let's learn how to do that.  We will assume that you have tables, numpy, pylab, and the statistics package installed in your python installation.

```
f=tables.openFile("SingleParticle.0.h5")
data= f.root.Observables.Energy.Kinetic
data= numpy.array(data)
myStats=stats.Stats(data)
print "My kinetic energy is ",str(myStats[0])+" "+str(myStats[2])
pylab.plot(data)
pylab.xlabel("MC Time")
pylab.ylabel("Kinetic Energy")
pylab.show()
```

Check to see that you get the right answer for the energy of a particle in a box!

Other variables that the energy observable creates includes:
  * Kinetic
  * dUShort
  * dULong

Interesting things to note:
  * NumImages is important for small boxes
  * There is no time step error here!


Once you feel familiar with this introductory simulation, we can move towards a more interesting system:
[helium](Liquid.md)