#!/bin/bash
k=$1
kb=$2
outdir=$3
ncores=$4
folder=$5
pfile1=$6
pfile2=$7

cd ${folder}
mkdir -p ./tmp3/
mkdir -p ./tmp3/init/

mpirun -np ${ncores} cheartsolver.out ${pfile1} -\#k=${k} -\#kb=${kb}

mpirun -np ${ncores} cheartsolver.out ${pfile2} -\#k=${k} -\#kb=${kb} -\#outdir=${outdir}
