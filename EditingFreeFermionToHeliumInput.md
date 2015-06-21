We are converting the Single Particle input file to one that is useful for Helium.

There are four key steps to converting the input file.
  1. Changing System Parameters (particle number, type, name, lambda and simulation box size)
  1. Switching particle statistics from boltzmannon to boson
  1. Getting the correct pair interaction for the Helium atoms (and choosing a time step )
  1. Setting up observables.

## Changing System Parameters ##

To accomplish this, first let's change the species
<font color='blue'> Type </font> (Section: System/Particles/Species) from "free" to "He". Let's also change <font color='blue'> lambda </font> (Section: System/Particle/Species) to the correct value of lambda (the inverse mass) for He-4 which is 6.059615.

We will do our simulation with 10 Helium particles (in practice, this is too small to get accurate results for the bulk, but we want our simulation to run in a reasonable time frame). Therefore, set <font color='blue'> NumParticles </font> (Section: System/Particles/Species) to 10. Standard volume and pressure (SVP) is a number density of 0.02182 $\AA^{-3}$, which corresponds to a <font color='blue'> Box </font> (Section: System) with sides 7.709894 $\AA$.

Now, we want the output file to be named differently. Therefore, change the <font color='blue'> OutFileBase </font> (Section: Observables) to Helium.

Finally, there is the <font color='blue'>Name </font> of the Species used throughout the input file. It is specified by you, the user, and can be anything (myParticle, He, free, etc.) as long as it is consistent throughout the input file. (Name of a <font color='blue'> Species</font> is distinct from <font color='blue'>Type</font> for reasons that will become apparent later). The <font color='blue'> Name</font> is currently set to "free", which is somewhat confusing for He-4 particles. Instead, change all (three) instances of "free" to "Helium4" in the input (Name (Section: System/Particles/Species), Species (Section: Moves/Move) and Species (Section: Observables/Observable))
[edit](edit.md)Particle Statistics

## Doing Bosonic Statistics ##
Helium is a boson. In path integral Monte Carlo, it is typical to sample the bosonic density matrix by explicitly integrating over permutations , or "sampling permutations." This means that the Metropolis formalism is employed to propose and probabilistically accept exchanges of beads between different particles. When an exchange (permutation) is accepted, two or more polymers become linked to form a single polymer.

We must now tell PIMC++ to use bosons. For the single free particle, the <font color='blue'> Statistics </font> variable was set to "BOLTZMANNON", meaning all particles (if there were more than one) are distinguishable. Inside the Helium species section, change Statistics (Section: System/Particles/Species) from "BOLTZMANNON" to "BOSON".  Now, there are different moves that can sample over permutations. We have to specify how we want this sampling done. For this tutorial, we will use the BisectionBlock move (like we did for free particles) with the addition that it will attempt to sample over all four particle permutations. To modify the algorithm appropriately, we must go to the BisectionBlock move and change PermuteType (Section: Moves/Move) from "NONE" to "TABLE".

When the <font color='blue'> PermuteType </font> is "TABLE", we have to specify some additional parameters. Here we will give a brief overview of these parameters, but for an in-depth explanation of them, please look here in the documentation. The table samples all possible 1,2,3, and 4 particle permutations with a probability (approximately) proportional to their kinetic action. Because 2, 3, and 4 particle permutations are important but typically have a much smaller kinetic action, we may want to increase the likelihood that they are sampled. We do this by multiplying their probability by some factor. In the input, we want to add the variable <font color='blue'> Array<double,1> Gamma(4) </font> (Section: Moves/Move) parameter and set it to `[1.0,20.0,100.0,400.0]`. This means that four particle permutations will be attempted with a likelihood 400 times greater than if only the kinetic action were used to compute the transition probability. (nb: all these numbers must be greater then 1)

The other important parameter is <font color='blue'>double epsilon </font>(Section: Moves/Move). This value specifies the smallest probability that is kept in the table. Because the table is naively $N^4$, this is a useful numerical approximation that dramatically increases the efficiency of the algorithm. Add this variable and set it to 1e-5. (Note that doing so introduces a small numerical approximation into the code)
Remember that if you are editing the input file as you go, you will need to add a semicolon ";" at the end of each new line of input to conform to standard C++ syntax!

## Generating density matrix (or how to build a PairAction) ##

Path Integral Monte Carlo calculates integrals of the following form:
$$\frac{1}{Z} \int \rho(R,R;\beta) O(R) dR$$ where $\rho(R)$ is the many body density matrix and O(R) are some observables. However, we can use a good approximation to $\rho(R,R;\tau)$, the high temperature many body density matrix (remember that $\beta$ and $\tau$ are inverse temperatures so small $\tau$ corresponds to high temperature).

This is done because we don't know how to accurately approximate , much less write an exact expression. However, we can use a good approximation for , the high temperature many body density matrix (remember that  and are inverse temperatures, so small  corresponds to a high temperature).
We will now discuss how (and why) we approximate the high temperature density matrix. The high temperature density matrix is exactly .
Because we don't have exact eigenvalues for the operator  with an arbitrary potential, we have to introduce an approximation. The standard approach is to do a Trotter breakup, disregarding the commutator between the kinetic and potential action operators, and write it approximately as
. This is exactly  where  is the (Helium) potential. Generally, this approximation introduces an error of order  but is exact in the limit  and . Therefore, we must use a very small tau to maintain the validity of the approximation.
We can produce a better approximation whose error is higher order in , allowing us to use a much larger value and hence require fewer time slices  for a given  (remember ). We do this by using the approximation that the many body density matrix is a product of pair density matrices  . We can (effectively exactly) calculate these pair density matrices between any two types of particles in our system. Because we are doing a system of bulk Helium-4, it has only one possible type of pair interaction, and we are only required to specify/create the pair density matrix for He4-He4. Constructing this interaction is an important but technically involved process. To learn how to produce the relevant files, please see Constructing a density matrix . If you want to skip this step for now, you may use the already prepared files [He4.95.OffDiag.PairAction](http://www.physics.princeton.edu/~bclark/Tutorials/liquidHelium/He4.95.OffDiag.PairAction) and [He4.95.dm](http://www.physics.princeton.edu/~bclark/Tutorials/liquidHelium/He4.95.dm)) (Note: Make sure all these files are in the in the directory from which you are running the code).

After constructing the files, we need to tell the code where to find them. We do this by changing the variable <font color='blue'> Array<string,1> PairActionFiles(1) </font> (Section: Action) and setting it to `["He4.95.OffDiag.PairAction"]`. This will ensure that the code can find the PairAction between He4 particles and other He4 particles. If we had other species Types in our system, we would have to have a PairAction for each species type.

Finally, let's set the time step . Because we have produced a particularly accurate many body density matrix, we can use a much higher value for tau then we would have otherwise. For, this tutorial let us set <font color='blue'> tau </font> (Section: System) to be 0.05. (This is actually still somewhat large, but we want our tutorial to run quickly and will check for convergence later). Remember when you choose tau that the value you select must be in the density matrix file as well as all values of tau required by the BisectionBlock Move (i.e. all powers up to `2**NumLevels*tau`).


## Adding observables ##

At this point you should have a basic helium input file and should now be able to run the Helium simulation. Although we haven't set up anything to measure yet, picked good moves, or tuned the Algorithm, run the Helium simulation now to verify we've successfully changed our input.

Now we would like to measure some aspects of our simulations. Many of these measurements (Observables) will have to do with Helium, although some of them deal with properties of our simulation, such as the efficiency of the algorithm. We will start by setting up the observables, and then we will add them to the algorithm. We will now set up a number of observables.

Let us start with the Structure Factor. We begin with adding a
```
Section (Observable)
```
inside of our Observables section. Its <font color='blue'> Type </font>(Section: Observables/Observable) will be "StructureFactor". As before, we can choose any <font color='blue'> Name </font> (Section: Observables/Observable) we desire but then must consistently use it throughout the rest of the input. For this tutorial, let's set the Name to "HeliumStructureFactor". In order for the structure factor to work appropriately, we have to specify the k points at which  will be evaluated. There are two ways to do this. To begin with, you can set a k-cutoff. If this is done, all k points whose magnitude is less than this k-cutoff will be calculated (if a long range action is being used, the k-cutoff from the long range action is used automatically). This <font color='blue'> kCutoff </font> (Section: Observables/Observable) is set by specifying
```
double kCutoff=5.0
```

The other method is to specify a list of k-vectors explicitly. See the documentation on the Structure Factor for the latter. Here we will just specify a k-cutoff of 5.0. We now must set the Frequency which specifies how many times the observable is called in the algorithm before it is measured. We will start by specifying the Frequency (Section: Observables/Observable) to 25. Later, if we find that this is taking too much time (or statistics are not good) in our simulation, we will change this number. Finally, specify the strings Species1 and Species2 (Section: Observables/Observable) equal to "Helium4". At this point your Structure Factor observable should look like this:
```
  Section (Observable)
     {
       string Type="StructureFactor";
       string Name="HeliumStructureFactor";
       double kCutoff=5.0;
       int Frequency=25; 
       string Species1="Helium4";
       string Species2="Helium4";

     }

```

By following the [documentation](http://cms.mcc.uiuc.edu/pimcpp/index.php/Observables) and this example, add observables for the Superfluid Fraction , PathDump (set the Frequency to 100), and Pair Correlation Function.

Finally, let's specify the Time Analysis observable. The code automatically calculates how much time is spent in each observable and move in the system. By adding the Time Analysis observable of our system, we are able to access this information and can then do a more effective job of tuning the simulation.

Now that we've defined our observables, we have to place them into the algorithm. To do this, let's go to the Algorithm section of our simulation. Here we can add in the new observables we have just declared underneath the current observable. Make sure you use the "Name" that you specified when you defined the observables.

Now, go back and continue with the main tutorial.