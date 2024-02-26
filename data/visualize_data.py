#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Aug  7 13:04:52 2023

@author: Javiera Jilberto Vallejos
"""

import cheartio as chio
import meshio as io

# mesh_path =
mesh = chio.read_mesh('../mesh/lv_model', meshio = True)
# mesh.point_data['fiber60'], mesh.point_data['sheet60'], mesh.point_data['sheetnormal60'] = chio.read_fibers('../data/fiber60.field')
# mesh.point_data['fiber75'], _, _ = chio.read_fibers('fiber75.field')
mesh.point_data['at60'] = chio.read_dfile('../data/at60.FE')
# mesh.point_data['at75'] = chio.read_dfile('at75.FE')
mesh.point_data['trans'] = chio.read_dfile('../data/trans.FE')
mesh.cell_data['endo'] = [chio.read_dfile('../data/endo.FE')]
mesh.cell_data['mid'] = [chio.read_dfile('../data/mid.FE')]
mesh.cell_data['epi'] = [chio.read_dfile('../data/epi.FE')]
mesh.cell_data['aha'] = [chio.read_dfile('../data/aha.FE')]
# mesh.cell_data['layers'] = [chio.read_dfile('layers.FE')]

# # Apex distance
mesh.point_data['apex'] = chio.read_dfile('../data/apex.FE')

io.write('data.vtu', mesh)

