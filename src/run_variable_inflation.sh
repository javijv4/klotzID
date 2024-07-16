#!/bin/bash
k=$1
kb=$2
parlv=$3
parrv=$4
outdir=$5
ncores=$6
folder=$7
pfile1=$8
pfile2=$9
edv_lv=$10
edp_lv=$11
edv_rv=$12
edp_rv=$13


cd ${folder}

mpirun -np ${ncores} cheartsolver.out ${pfile1} -\#k=${k} -\#kb=${kb} -\#par_LV=${parlv} -\#par_RV=${parrv} -\#outdir=${outdir}  --pedantic-printing

mpirun -np ${ncores} cheartsolver.out ${pfile2} -\#k=${k} -\#kb=${kb} -\#outdir=${outdir}  --pedantic-printing