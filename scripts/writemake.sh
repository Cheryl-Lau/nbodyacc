#!/bin/bash

mynbodyaccdir=${0%/*/*}

echo 'NBODYDIR='$mynbodyaccdir;

makeflags='RUNDIR=${PWD}';

echo 'again:'
echo '	cd ${NBODYDIR}; make '$makeflags'; cd -; cp ${NBODYDIR}/bin/nbodyacc .'

echo '.PHONY: nbodyacc nbodyaccsetup'
echo 'nbodyacc         : again'
echo 'nbodyaccsetup    : setup'

echo 'setup:'
echo '	cd ${NBODYDIR}; make '$makeflags' setup; cd -; cp ${NBODYDIR}/bin/nbodyaccsetup .'
echo 'clean:'
echo '	cd ${NBODYDIR}; make clean'
