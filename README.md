
NbodyAcc 
--------

A simple N-body code for modelling stellar dynamics in cores/clusters and gas accretion.  \
Point masses can be treated as binaries. 

Installation 
-------------

To get a copy, fork to create a clone in your github, then enter the following into your terminal: \
`git clone https://github.com/USERNAME/nbodyacc.git`

It is recommended to add the following lines into your .bashrc or .bash-profile: \
`export SYSTEM=gfortran` (or `export SYSTEM=ifort`) \
`export OMP_SCHEDULE="dynamic"` \
`export OMP_STACKSIZE=512M` \
`export KMP_STACKSIZE=128m` \
`ulimit -s unlimited` \
`alias getmake='~/nbodyacc/scripts/writemake.sh > Makefile'`

Don't forget to type `source .bashrc` after editing the file, or restart the terminal.  

Running a simulation 
--------------------

Go to `nbodyacc/build/Makefile` and enter the name of your setup file in: \
`SETUPFILE= <name of setup>.f90`

Choose whether you want self-gravity to be computed: `GRAVITY=yes/no`, \
and whether the point masses are binary pairs `BINARY=yes/no`.

Now create a work-directory anywhere _outside_ this code's directory. The work directory is where the simulation outputs are stored. Create a new work-directory for every new run. 

Navigate to a work-directory of your choice and enter: \
`getmake` \
to create the local Makefile.

To compile, enter: \
`make; make setup` 

To clear the executables and re-compile after modifying the code, enter: \
`rm nbody*; make clean`

Modify/create your setup file in `src/setup/setup_*.f90`, then use it to generate the initial dumpfile by running the command: \
`./nbodyaccsetup` \
in your work-directory.
This creates a dumpfile `ptmass_00000.tmp` and an input file `input_params.in` which contains the defaults of the runtime parameters. Modify if necessary. 

