#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 11:29:50 2024

@author: Javiera Jilberto Vallejos
"""

from KlotzID import KlotzID

ncores = 3   # WARNING!!! The program uses ncores x 2 cores.
out_fldr = 'out_mh'
out_fname = 'optimize_params.P'
inflation_type = 'inverse_volume'
pfile = 'inverse_volume_lv.P'

pressure_var = 'LV_LM'
volume_var = 'LV_Vol'

times = (1,100,1)

ed_pressure = 1.333  # kPa
ed_volume = 129659.88361282746  # mm3

k0 = 0.010179
kb0 = 0.696370
params0 = (k0, kb0)

klotzopt = KlotzID(pfile, pressure_var, volume_var, out_fldr, times,
                   ed_pressure, ed_volume, inflation_type, ncores, plot_intermediate=True)
klotzopt.max_iterations = 10
params = klotzopt.optimize(params0)
klotzopt.run_last_simulation(params)
klotzopt.write_params(out_fname, params)

# import cheartio as chio
# import matplotlib.pyplot as plt
# self = klotzopt
# pres = chio.read_scalar_dfiles('{}/{}/{}'.format(self.cheart_folder, 'tmp1', self.pressure_var), self.times)[::-1]
# vol = chio.read_scalar_dfiles('{}/{}/{}'.format(self.cheart_folder, 'tmp1', self.volume_var), self.times)
# pres_eps = chio.read_scalar_dfiles('{}/{}/{}'.format(self.cheart_folder, 'tmp2', self.pressure_var), self.times)[::-1]

# plt.figure(1, clear=True)
# plt.plot(vol/1000, self.klotz_function(vol)*7.50062, 'k', label='klotz')
# plt.plot(vol/1000, pres*7.50062, 'r', label='it={:d}'.format(self.it))
# plt.plot(vol/1000, pres_eps*7.50062, 'r', label='it={:d}'.format(self.it))
# plt.show()
