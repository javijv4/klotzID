#!/bin/bash

cheartsolver.out prep.P -#meshdir=$1 --prep
cheartsolver.out cheart_compute_bv_volume.P -#meshdir=$1