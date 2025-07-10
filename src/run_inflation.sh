#!/bin/bash
k=$1
kb=$2
par_LV=$3
par_RV=$4
outdir=$5
ncores=$6
folder=$7
pfile=$8
lv_edv=$9
rv_edv=${10}
lv_v0=${11}
rv_v0=${12}
meshdir=${13}
datadir=${14}

cd ${folder}
mpirun -np ${ncores} cheartsolver.out ${pfile} -\#k=${k} -\#kb=${kb} -\#par_LV=${par_LV} -\#par_RV=${par_RV} -\#lv_edv=${lv_edv} -\#rv_edv=${rv_edv} -\#lv_v0=${rv_v0} -\#lv_v0=${rv_v0} -\#outdir=${outdir} --pedantic-printing
