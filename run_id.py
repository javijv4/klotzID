#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 11:29:50 2024

@author: Javiera Jilberto Vallejos
"""


from KlotzID import KlotzID

ncores = 3
out_fldr = 'out_mh'
inflation_type = 'volume'

pressure_var = 'LV_LM'
volume_var = 'LV_Vol'

times = (1,10,1)

ed_pressure = 1.33  # kPa
ed_volume = 131153.376791549788322  # mm3

k0 = 1.0
kb0 = 1.0
params0 = (k0, kb0)

klotzopt = KlotzID(out_fldr, pressure_var, volume_var, times, ed_pressure, ed_volume, inflation_type, ncores)
klotzopt.max_iterations = 2
params = klotzopt.optimize(params0)

