#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 11:29:50 2024

@author: Javiera Jilberto Vallejos
"""

from KlotzID import KlotzID

ncores = 5   # WARNING!!! The program uses ncores x 2 cores.
out_fldr = 'out_mh'
out_fname = 'optimize_params_mh.P'
inflation_type = 'volume'
pfile = 'volume_inflation.P'

pressure_var = 'LV_LM'
volume_var = 'LV_Vol'

times = (1,100,1)

ed_pressure = 2  # kPa
ed_volume = 131153.376791549788322  # mm3

k0 = 0.688901
kb0 = 0.419691
params0 = (k0, kb0)

klotzopt = KlotzID(pfile, pressure_var, volume_var, out_fldr, times, 
                   ed_pressure, ed_volume, inflation_type, ncores, plot_intermediate=True)
klotzopt.max_iterations = 20
params = klotzopt.optimize(params0)
klotzopt.run_last_simulation(params)
klotzopt.write_params(out_fname, params)