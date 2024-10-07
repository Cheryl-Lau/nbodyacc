
NbodyAcc 
--------

A simple N-body code for modelling stellar dynamics in cores/clusters and gas accretion.  \
Point masses can be treated as binaries. 

Installation 
-------------

To get a copy, fork the repo. Then on your terminal, navigate to home directory and enter: \
`git clone https://github.com/<your git username>/nbodyacc.git`

It is recommended to add the following lines into your `.bashrc` (linux) or `.bash_profile` (mac): \
`export SYSTEM=gfortran` (or `export SYSTEM=ifort`) \
`export OMP_SCHEDULE="dynamic"` \
`export OMP_STACKSIZE=512M` \
`export KMP_STACKSIZE=128m` \
`ulimit -s unlimited` \
`alias getmake='~/nbodyacc/scripts/writemake.sh > Makefile'`

Don't forget to type `source .bashrc` after editing the file, or restart the terminal.  

Compiling the code  
------------------

Go to `nbodyacc/build/Makefile` and enter the name of your setup file in: \
`SETUPFILE= <name of setup>.f90` \
This setup file should be stored in `nbodyacc/src/setup/`.

Choose whether you want self-gravity to be computed: `GRAVITY=yes/no`, \
and whether the point masses are binary pairs: `BINARY=yes/no`.

Now create a work-directory anywhere _outside_ this code's directory. The work directory is where the simulation outputs are stored. Create a new work-directory for every new run. 

Navigate to a work-directory of your choice and enter: \
`getmake` \
to create the local Makefile.

To compile, enter: \
`make; make setup` 

To re-compile after modifying the source code, enter: \
`rm nbody*; make clean; make; make setup` \
which clears the previously generated executables. 

Running a simulation 
--------------------

First, in your work-directory, run the command: \
`./nbodyaccsetup` \
This creates the initial dumpfile `ptmass_00000.tmp` and an input file `input_params.in` which contains the defaults of the runtime parameters. Modify if necessary. 

To begin a simulation with this input file, run: \
`./nbodyacc` 


**HAPPY DEBUGGING**

