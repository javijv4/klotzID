#!/bin/bash
k=$1
kb=$2
outdir=$3
ncores=$4

mpirun -np ${ncores} cheartsolver.out inflation_mh.P -\#k=${k} -\#kb=${kb} -\#outdir=${outdir}
