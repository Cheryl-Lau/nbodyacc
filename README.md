
NbodyAcc 
--------

A simple N-body code for modelling stellar dynamics in cores/clusters and gas accretion.  \
Point masses can be treated as binaries. 

Installation 
-------------

To get a copy, fork to create a clone then enter the following into your terminal: \
`git clone https://github.com/USERNAME/nbodyacc.git`

It is recommended to add the following lines into your .bashrc or .bash-profile: \
`export SYSTEM=gfortran` (or `export SYSTEM=ifort`) \
`export OMP_SCHEDULE="dynamic"` \
`export OMP_STACKSIZE=512M` \
`export KMP_STACKSIZE=128m` \
`ulimit -s unlimited` \
`alias getmake='~/nbodyacc/scripts/writemake.sh > Makefile'`

Don't forget to `source .bashrc`

Running a simulation 
--------------------

Navigate to the work-directory, and enter: \
`getmake` \
to create the local Makefile.

To compile, enter: \
`make; make setup` 

To clear the executables after modifying the code, enter: \
`make clean`
