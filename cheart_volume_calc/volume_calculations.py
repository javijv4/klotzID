#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 19 14:40:18 2024

@author: Javiera Jilberto Vallejos
"""

import numpy as np
import matplotlib.pyplot as plt
import os
import shutil
from subprocess import Popen

def run_cheart_volume_compute(mesh_fldr):
    # Copy boundary file because CH is dumb
    shutil.copy(mesh_fldr + 'boundaries.P', '.')

    # Run cheart
    with open('cheart.log', 'w') as ofile:
        p = Popen(['bash', 'run_compute_volume.sh', mesh_fldr], stdout=ofile, stderr=ofile)
        p.wait()

    # Reading results
    lv_volume = np.loadtxt('volume_lv.norm', usecols=1)
    rv_volume = np.loadtxt('volume_rv.norm', usecols=1)

    # Cleaning up
    os.remove("boundaries.P")
    os.remove("volume_lv.norm")
    os.remove("volume_rv.norm")

    return lv_volume, rv_volume

def klotz_v0(lv_volume, lv_pressure):  # in mm3, KPa
    Vm = lv_volume/1000
    Pm = lv_pressure*7.50062

    return Vm*(0.6 - 0.006*Pm)*1000



mesh_fldr = '../work/test_fine/mesh/'
lv_pressure = 1.33
rv_pressure = 0.4
lv_volume, rv_volume = run_cheart_volume_compute(mesh_fldr)
lv_v0 = klotz_v0(lv_volume, lv_pressure)
rv_v0 = rv_volume-(lv_volume - lv_v0)*rv_pressure/lv_pressure          # We assume the RV volume will decrease by the same factor as the LV

# Saving results
f = open(mesh_fldr + 'volumes.P', "w")
f.write('#LV_Ved={}\n'.format(lv_volume))
f.write('#RV_Ved={}\n'.format(rv_volume))
f.write('#LV_V0={}\n'.format(lv_v0))
f.write('#RV_V0={}\n'.format(rv_v0))
f.close()