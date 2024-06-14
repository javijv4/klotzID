#!/bin/bash
k=$1
kb=$2
outdir=$3
ncores=$4
folder=$5
pfile1=$6
pfile2=$7
edv_lv=$8
edp_lv=$9
edv_rv=$10
edp_rv=$11


cd ${folder}

mpirun -np ${ncores} cheartsolver.out ${pfile1} -\#k=${k} -\#kb=${kb} -\#outdir=${outdir}  --pedantic-printing

mpirun -np ${ncores} cheartsolver.out ${pfile2} -\#k=${k} -\#kb=${kb} -\#outdir=${outdir}  --pedantic-printing