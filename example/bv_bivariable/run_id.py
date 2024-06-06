#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 11:29:50 2024

@author: Javiera Jilberto Vallejos
"""

from KlotzID import KlotzID

ncores = 8                              # WARNING!!! The program uses (ncores x 2) cores.
out_fldr = 'mesh/out'                   # Output folder
out_fname = 'mesh/optimize_params.P'    # Final parameters will be saved in this file
inflation_type = 'volume'               # Method: "volume" or "inverse_volume"
pfile = 'biv_volume_inflation.P'      # P-file
pfile_bv = 'biv_bivariable_inflation.P'      # P-file
pfile_bv_init = 'biv_bivariable_inflation_runinit.P'      # P-file

pressure_var = 'LV_LM'                  # Name of CH pressure variable
volume_var = 'LV_Vol'                   # Name of CH volume variable

times = (1,100,1)                       # Time at which pressure_var and volume_var are output

# ED values, left ventricle
ed_pressure = 2.480742718597106 # kPa
ed_volume = 222334.751  # mm3

# Initial guess
k0 = 1
kb0 = 1

par_lv=0
par_rv=0

params0 = (k0, kb0,par_lv,par_rv)

# Optimization
klotzopt = KlotzID(pfile, pressure_var, volume_var, out_fldr, times, 
                   ed_pressure, ed_volume, inflation_type, ncores, 
                   pfile_bv_init,pfile_bv,
                   plot_intermediate=True,      # Save a png of the pressure volume curve every iteration
                   save_intermediate=True,      # Save intermediate pressure, volume and parameters in an .npz file
                   lvrv=True)     
klotzopt.max_iterations = 20                    # Setting max number of iterations
params = klotzopt.optimize(params0)
klotzopt.run_last_simulation(params)
klotzopt.write_params(out_fname, params)