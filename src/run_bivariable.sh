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
mkdir -p ./tmp3/
mkdir -p ./tmp3/init/

mpirun -np ${ncores} cheartsolver.out ${pfile1} -\#k=${k} -\#kb=${kb} -\#EDV=${edv_lv} -\#EDV_rv=${edv_rv} 

mpirun -np ${ncores} cheartsolver.out ${pfile2} -\#k=${k} -\#kb=${kb} -\#outdir=${outdir} -\#EDV=${edv_lv} -\#EDV_rv=${edv_rv} -\#pres_lv=${edp_lv} -\#pres_rv=${edp_rv}
