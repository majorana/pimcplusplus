Here we will do a short exploration of the properties of a hydrogen molecule using PIMC++.

Earlier in the week we computed properties of the hydrogen molecule using Variational Monte Carlo. As another example of the sort of problems PIMC++ can do, we can compute the Energy of a hydrogen molecule. Look in the file H2.in
This file has been set up to calculate the energy and pair correlation function of an H2 molecule. Skim through this file. There are a couple of key differences between this file and the helium file. Notice the way that two species (protons and electrons) are setup. In this case, we have set the protons to be classical particles (by setting lambda to 0.0). Also, inside the action section, a series of PairActions exist for the coulumb action.
To run this code, we should call

```
pimc++ H2.in
```

(do this now)

Once the code has run, you can again analyze the data by calling

Report.py H2.0.h5

For this simulation we are using a bondlength of approximately 1.4 although it should be reasonably clear how to modify this.

You should notice the following things:

  1. the energy you get is within error bars of the correct answer of -1.1744.
  1. that the kinetic and potential energy obey the virial theorem:  . Note that this should only be true for coloumb systems

Note, the difference in the input between Variational Monte Carlo and Path Integral Monte Carlo. In the VMC calculation from Monday, we needed to guess a trial wave function and the quality of our result was dependent on this trial function. In PIMC, the only input is the coloumb potential (1/r), the bond length of the particles, and the temperature (which we currently have set sufficiently small to ensure we are in the ground state.) From this information, we are able to get the exact result. (Of course, if you were doing a simulation that involved many fermions, this is no longer true as you would need a trial wave function (density matrix) to specify the nodal restriction. Nonetheless, the dependence will still be smaller then in VMC.)

Now, let us visualize our simulation.

Call

```
pathvis++ H2.0.h5
```

You should see a graphical representation of the paths as well as the location of the ions. One interesting thing to note is that at this bond length, one can't tell from looking at the electron paths that there are two point charges instead of one point charge. In fact, everything looks largely spherical. This should give us good indication that using a trial function with a bond centered gaussian is physically very reasonable.

At this point, you are free to explore different things. Possible things to look at include

  1. energy as a function of bond length
  1. what temperature the molecule breaks down. This requires:
    * making the protons quantum (i.e. putting the correct value for lambda in the proton species section)
    * adding moves for the protons (Displace and BisectionBlock)
    * calling these moves in the Algorithm section.
If you are interested, we would be happy to help get you started in implementing these steps.

# Introduction #

Add your content here.


# Details #

Add your content here.  Format your content with:
  * Text in **bold** or _italic_
  * Headings, paragraphs, and lists
  * Automatic links to other wiki pages