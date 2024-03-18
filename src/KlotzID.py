#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Mon Feb 26 11:29:50 2024

@author: Javiera Jilberto Vallejos
"""

import os
import numpy as np
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
import cheartio as chio
from subprocess import Popen
import shutil
import time


class KlotzID:
    def __init__(self, pfile, pressure_var, volume_var, out_fldr, sim_times, ed_pressure, ed_volume, inflation_type, ncores, 
                 plot_intermediate=False, save_intermediate=False):
        self.self_path = os.path.dirname(os.path.abspath(__file__))
        self.cheart_folder = os.path.dirname(pfile)
        if self.cheart_folder == '': self.cheart_folder = '.'
        self.ncores = ncores
        self.inflation_type = inflation_type
        self.pfile = os.path.basename(pfile)

        # Output path
        self.out_fldr = out_fldr
        if not os.path.exists('{}/{}'.format(self.cheart_folder, out_fldr)):
            os.mkdir('{}/{}'.format(self.cheart_folder, out_fldr))

        self.pressure_var = pressure_var
        self.volume_var = volume_var

        self.times = sim_times

        # Compute target curve
        self.ed_pressure = ed_pressure  # kPa
        self.ed_volume = ed_volume # mm3
        self.klotz_volume, self.klotz_pressure = klotz_curve(ed_volume, ed_pressure)
        self.klotz_function = interp1d(self.klotz_volume, self.klotz_pressure, fill_value='extrapolate')

        # Optimization parameters
        self.eps = 1e-2
        self.lam = 1.0
        self.lam_rand = [1.0, 5.0]
        self.max_iterations = 10
        self.gnorm = 1e3              # Pressure error norm
        self.it = 0

        # Others
        self.plot_intermediate = plot_intermediate
        self.save_intermediate = save_intermediate
        self.intermediate = {'vol': [], 'pres': [], 'params': []}

        self.logfile = 'klotzid.log'


    def pre_clean(self):
        # Creating folders
        if not os.path.exists('{}/tmp1'.format(self.cheart_folder)):
            os.mkdir('{}/tmp1'.format(self.cheart_folder))
        if not os.path.exists('{}/tmp2'.format(self.cheart_folder)):
            os.mkdir('{}/tmp2'.format(self.cheart_folder))

        # Deleting log
        try:
            os.remove(self.logfile)
        except OSError:
            pass


    def post_clean(self):
        # Deleting temporaty folders.
        shutil.rmtree('{}/tmp1'.format(self.cheart_folder))
        shutil.rmtree('{}/tmp2'.format(self.cheart_folder))
        os.remove('tmp1.log')
        os.remove('tmp2.log')


    def optimize(self, params, post_clean=False):
        self.pre_clean()

        # Initializing variables
        error = 1e3
        self.it = 0
        start_time = time.time()
        while (error > 1e-3) and (self.it < self.max_iterations):
            new_params, error = self.optimize_iteration_volume(params)

            self.write_log(params, error)

            # update
            params = new_params
            self.it += 1

        if self.it < self.max_iterations:
            end_time = time.time()
            print('Optimization succesful in {:d} iterations and {:2.3f} s'.format(self.it, end_time-start_time))
            print('Parameters found are k={:f} and kb={:f}'.format(params[0], params[1]))
        else:
            print('Optimization failed.')

        if self.save_intermediate:
            int_vol = np.vstack(self.intermediate['vol'])
            int_pres = np.vstack(self.intermediate['pres'])
            int_params = np.vstack(self.intermediate['params'])
            np.savez('{}/{}/intermediate_results.npz'.format(self.cheart_folder, self.out_fldr), volume=int_vol, pressure=int_pres, params=int_params)

        if post_clean:
            self.post_clean()
        return params


    def optimize_iteration_volume(self, params):
        k, kb = params

        p1 = self.run_cheart_inflation((k, kb), 'tmp1')
        p2 = self.run_cheart_inflation((k, kb+self.eps), 'tmp2')
        exit_codes = [p.wait() for p in (p1, p2)]

        # Load results
        pres = chio.read_scalar_dfiles('{}/{}/{}'.format(self.cheart_folder, 'tmp1', self.pressure_var), self.times)
        vol = chio.read_scalar_dfiles('{}/{}/{}'.format(self.cheart_folder, 'tmp1', self.volume_var), self.times)
        pres_eps = chio.read_scalar_dfiles('{}/{}/{}'.format(self.cheart_folder, 'tmp2', self.pressure_var), self.times)

        if self.inflation_type == 'inverse_volume':
            vol = vol[::-1]
            vol = np.append(vol, self.ed_volume)
            pres = np.append(0., pres)
            pres_eps = np.append(0., pres_eps)

        if self.plot_intermediate:
            self.plot_inflation_curve('{}/{}/klotz_it{:d}.png'.format(self.cheart_folder, self.out_fldr, self.it), vol, pres, params)

        if self.save_intermediate:
            self.intermediate['vol'].append(vol)
            self.intermediate['pres'].append(pres)
            self.intermediate['params'].append(params)
            
        # Parameter update
        k, kb = params

        # Linear correction
        k_delta = self.ed_pressure/pres[-1]
        k = k*k_delta
        pres = pres*self.ed_pressure/pres[-1]
        pres_eps = pres_eps*self.ed_pressure/pres_eps[-1]

        # Levenber-marquadt iteration
        g = self.compute_curve_difference(vol, pres)
        kb_delta = self.LM_update(g, pres, pres_eps)
        kb += kb_delta

        # Compute error
        error = np.abs(kb_delta) + np.abs(k_delta-1)

        # Return new parameters
        params = (k, kb)
        return params, error


    def LM_update(self, g, pres, pres_eps):
        # Compute jacobian (forward differences)
        J = (pres_eps - pres)/self.eps
        J = J[self.data_mask,None]       # Making it a matrix

        # Solve LM system
        mat = J.T@J + self.lam*np.eye(J.shape[1])
        rhs = J.T@g

        delta = np.linalg.solve(mat, rhs)

        # Update lam
        gnorm = np.linalg.norm(g)
        if self.it > 0:
            if gnorm < self.gnorm:
                self.lam *= gnorm/self.gnorm
            else:
                fac = np.random.uniform(self.lam_rand[0], self.lam_rand[1])
                print("Warning: No descent direction! Applying random lambda scaling. Factor: %.4f" % (fac))
                self.lam *= fac

        self.gnorm = gnorm

        return delta[0]


    def run_last_simulation(self, params, plot=True):
        print('Running simulation with optimized parameters')
        p1 = self.run_cheart_inflation(params, self.out_fldr)
        p1.wait()

        # Load results
        pres = chio.read_scalar_dfiles('{}/{}/{}'.format(self.cheart_folder, self.out_fldr, self.pressure_var), self.times)
        vol = chio.read_scalar_dfiles('{}/{}/{}'.format(self.cheart_folder, self.out_fldr, self.volume_var), self.times)

        if self.inflation_type == 'inverse_volume':
            vol = vol[::-1]
            vol = np.append(vol, self.ed_volume)
            pres = np.append(0., pres)

        ed_error = np.abs(self.ed_pressure/pres[-1]-1)

        g = self.compute_curve_difference(vol, pres)
        curve_error = np.linalg.norm(g)

        if plot:
            self.plot_inflation_curve('{}/{}/klotz_fit.png'.format(self.cheart_folder, self.out_fldr), vol, pres, params)
        print('Final simulation results: ED error = {:f}, curve error {:f}'.format(ed_error, curve_error))


    def compute_curve_difference(self, vol, pres):
        vol0 = vol[0]
        voled = vol[-1]
        vol_lim = (voled-vol0)/3 + vol0
        
        self.data_mask = vol > vol_lim
        vol_cut = vol[self.data_mask]
        pres_cut = pres[self.data_mask]

        pres_klotz = self.klotz_function(vol_cut)

        return pres_klotz - pres_cut


    def write_log(self, params, error):
        with open(self.logfile, "a") as file:
            # Writing data to a file
            file.write("Iteration {:d}, k={:f}, kb={:f}, error = {:e}\n".format(self.it, params[0], params[1], error))


    @staticmethod
    def write_params(fname, params):
        with open(fname, "w") as file:
            # Writing data to a file
            file.write("#k={:12.12f}\n".format(params[0]))
            file.write("#kb={:12.12f}".format(params[1]))


    def run_cheart_inflation(self, params, outdir):
        k, kb = params

        # Run cheart
        with open('{}.log'.format(outdir), 'w') as ofile:
            p = Popen(['bash', '{}/run_inflation.sh'.format(self.self_path), '{:f}'.format(k), '{:f}'.format(kb),
                       outdir, '{:d}'.format(self.ncores), self.cheart_folder, self.pfile],
                       stdout=ofile, stderr=ofile)

        return p


    def plot_inflation_curve(self, fname, volume, pressure, params):
        plt.figure(1, clear=True)
        plt.plot(volume/1000, self.klotz_function(volume)*7.50062, 'k', label='klotz')
        plt.plot(volume/1000, pressure*7.50062, 'r', label='it={:d}'.format(self.it))
        plt.annotate('k={:2.3e}\nkb={:2.3e}'.format(params[0], params[1]), (0.1,0.9), xycoords='axes fraction')
        plt.xlabel('Volume [ml]')
        plt.ylabel('Pressure [kPa]')
        plt.savefig(fname, bbox_inches='tight')



def klotz_V0(Vm, Pm):
    return Vm*(0.6 - 0.006*Pm)

def klotz_V30(V0, Vm, Pm, An = 27.78, Bn = 2.76):
    return V0 + (Vm-V0)/(Pm/An)**(1/Bn)

def klotz_ab(V30, Vm, Pm):
    beta = np.log(Pm/30)/np.log(Vm/V30)
    alpha = 30/V30**beta
    return alpha, beta

def klotz_curve(Vm, Pm, nv=100):  # inputs are in kPa mm3
    # Tranform inputs
    Vm = Vm/1000
    Pm = Pm*7.50062

    # Calculating curve
    V0 = klotz_V0(Vm, Pm)
    V30 = klotz_V30(V0, Vm, Pm)
    alpha, beta = klotz_ab(V30, Vm, Pm)

    V = np.linspace(0, Vm, nv)
    P = alpha*V**beta

    # Return in mm3, kpa
    return V*1000, P*0.133322
