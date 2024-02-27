#!/bin/bash
k=$1
kb=$2
outdir=$3
ncores=$4
folder=$5
pfile=$6

cd ${folder}
mpirun -np ${ncores} cheartsolver.out ${pfile} -\#k=${k} -\#kb=${kb} -\#outdir=${outdir}
