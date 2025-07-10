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
lv_edv=${10}
lv_edp=${11}
rv_edv=${12}
rv_edp=${13}
lv_v0=${14}
rv_v0=${15}
meshdir=${16}
datadir=${17}


cd ${folder}

mpirun -np ${ncores} cheartsolver.out ${pfile1} -\#k=${k} -\#kb=${kb} -\#lv_edv=${lv_edv} -\#rv_edv=${rv_edv} -\#lv_v0=${lv_v0} -\#rv_v0=${rv_v0} -\#outdir=${outdir} -\#meshdir=${meshdir} -\#datadir=${datadir} --pedantic-printing

mpirun -np ${ncores} cheartsolver.out ${pfile2} -\#k=${k} -\#kb=${kb} -\#lv_edp=${lv_edp} -\#lv_edv=${lv_edv} -\#rv_edp=${rv_edp} -\#rv_edv=${rv_edv} -\#lv_v0=${lv_v0} -\#rv_v0=${rv_v0}  -\#outdir=${outdir} -\#meshdir=${meshdir} -\#datadir=${datadir} --pedantic-printing