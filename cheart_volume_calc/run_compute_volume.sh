#!/bin/bash
path='../../../Desmoplakin/Models/DSPPatients/'

for patient in AB-17 AS-10 AV-19 BI-18 CA-15 CW-21 DM-23 JL-3 JN-8 KL-4 KL-5 KR-13 MB-16 SL-16 TS-9 VB-1 ZS-11
do
    meshdir=${path}/${patient}/mesh/
    cheartsolver.out cheart_compute_bv_volume.P -\#meshdir=$meshdir
done
