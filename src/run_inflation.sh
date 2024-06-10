#!/bin/bash
k=$1
kb=$2
par_LV=$3
par_RV=$4
outdir=$5
ncores=$6
folder=$7
pfile=$8
ed_volume=$9

cd ${folder}
mpirun -np ${ncores} cheartsolver.out ${pfile} -\#k=${k} -\#kb=${kb} -\#outdir=${outdir} -\#par_LV=${par_LV} -\#par_RV=${par_RV} -\#EDV=${ed_volume} 
