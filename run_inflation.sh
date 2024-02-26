#!/bin/bash
k=$1
kb=$2
outdir=$3
ncores=$4
ncores=$5

mpirun -np ${ncores} cheartsolver.out ${5} -\#k=${k} -\#kb=${kb} -\#outdir=${outdir}
