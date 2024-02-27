#!/bin/bash
k=$1
kb=$2
outdir=$3
ncores=$4
pfile=$5

mpirun -np ${ncores} cheartsolver.out ${pfile} -\#k=${k} -\#kb=${kb} -\#outdir=${outdir}
