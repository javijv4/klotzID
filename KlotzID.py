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
    def __init__(self, pfile, pressure_var, volume_var, out_fldr, sim_times, ed_pressure, ed_volume, inflation_type, ncores, plot_intermediate=False):
        self.ncores = ncores
        self.inflation_type = inflation_type
        self.pfile = pfile

        # Output path
        self.out_fldr = out_fldr
        if not os.path.exists(out_fldr): os.mkdir(out_fldr)

        self.pressure_var = pressure_var
        self.volume_var = volume_var

        self.times = sim_times

        # Compute target curve
        self.ed_pressure = ed_pressure  # kPa
        self.ed_volume = ed_volume # mm3
        self.klotz_volume, self.klotz_pressure = klotz_curve(ed_volume, ed_pressure)
        self.klotz_function = interp1d(self.klotz_volume, self.klotz_pressure)

        # Optimization parameters
        self.eps = 1e-2
        self.lam = 1.0
        self.lam_rand = [1.0, 5.0]
        self.max_iterations = 10
        self.gnorm = 1e3              # Pressure error norm
        self.it = 0

        # Others
        self.plot_intermediate = plot_intermediate
        self.logfile = 'klotzid.log'


    def pre_clean(self):
        # Creating folders
        if not os.path.exists('tmp1'): os.mkdir('tmp1')
        if not os.path.exists('tmp2'): os.mkdir('tmp2')

        # Deleting log
        try:
            os.remove(self.logfile)
        except OSError:
            pass


    def post_clean(self):
        # Deleting temporaty folders.
        shutil.rmtree('tmp1/')
        shutil.rmtree('tmp2/')
        os.remove('tmp1.log')
        os.remove('tmp2.log')


    def optimize(self, params):
        self.pre_clean()

        # Initializing variables
        error = 1e3
        self.it = 0
        start_time = time.time()
        while (error > 1e-3) and (self.it < self.max_iterations):
            new_params, error = self.optimize_iteration(params)

            self.write_log(params, error)

            # update
            params = new_params
            self.it += 1

        if self.it < self.max_iterations:
            end_time = time.time()
            print('Optimization succesful in {:d} iterations and {:2.3f} s'.format(self.it, end_time-start_time))
            print('Parameters found are k={:f} and kb={:f}'.format(params[0], params[1]))
            self.write_params(params)
            self.run_last_simulation(params)
        else:
            print('Optimization failed.')

        self.post_clean()
        return params


    def optimize_iteration(self, params):
        k, kb = params

        p1 = self.run_cheart_inflation((k, kb), 'tmp1')
        p2 = self.run_cheart_inflation((k, kb+self.eps), 'tmp2')
        exit_codes = [p.wait() for p in (p1, p2)]

        # Load results
        pres = chio.read_scalar_dfiles('{}/{}'.format('tmp1', self.pressure_var), self.times)
        vol = chio.read_scalar_dfiles('{}/{}'.format('tmp1', self.volume_var), self.times)
        pres_eps = chio.read_scalar_dfiles('{}/{}'.format('tmp2', self.pressure_var), self.times)
        
        # Parameter update
        k, kb = params
        
        # Linear correction
        k_delta = self.ed_pressure/pres[-1]
        k = k*k_delta
        pres = pres*self.ed_pressure/pres[-1]
        pres_eps = pres_eps*self.ed_pressure/pres_eps[-1]

        if self.plot_intermediate:
            self.plot_inflation_curve(vol, pres)

        # Levenber-marquadt iteration
        pres_klotz = self.klotz_function(vol)
        g = pres_klotz-pres
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
        J = J[:,None]       # Making it a matrix

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
    

    def run_last_simulation(self, params):
        print('Running simulation with optimized parameters')
        p1 = self.run_cheart_inflation(params, self.out_fldr)
        p1.wait()

        # Load results
        pres = chio.read_scalar_dfiles('{}/{}'.format(self.out_fldr, self.pressure_var), self.times)
        vol = chio.read_scalar_dfiles('{}/{}'.format(self.out_fldr, self.volume_var), self.times)
        ed_error = np.abs(self.ed_pressure/pres[-1]-1)

        pres_klotz = self.klotz_function(vol)
        g = pres_klotz-pres
        curve_error = np.linalg.norm(g)

        print('Final simulation results: ED error = {:f}, curve error {:f}'.format(ed_error, curve_error))


    def write_log(self, params, error):
        with open(self.logfile, "a") as file:
            # Writing data to a file
            file.write("Iteration {:d}, k={:f}, kb={:f}, error = {:e}\n".format(self.it, params[0], params[1], error))


    def write_params(self, params):
        with open('optimized_params.P', "w") as file:
            # Writing data to a file
            file.write("#k={:f}".format(params[0]))
            file.write("#kb={:f}".format(params[1]))


    def run_cheart_inflation(self, params, outdir):
        k, kb = params

        # Run cheart
        with open('{}.log'.format(outdir), 'w') as ofile:
            p = Popen(['bash', 'run_inflation.sh', '{:f}'.format(k), '{:f}'.format(kb), outdir, '{:d}'.format(self.ncores), self.pfile], stdout=ofile, stderr=ofile)

        return p
    

    def plot_inflation_curve(self, volume, pressure):
        plt.figure(1, clear=True)
        plt.plot(volume/1000, self.klotz_function(volume)*7.50062, 'k', label='klotz')
        plt.plot(volume/1000, pressure*7.50062, 'r', label='it={:d}'.format(self.it))
        plt.xlabel('Volume [ml]')
        plt.ylabel('Pressure [kPa]')
        plt.savefig('{}/klotz_it{:d}.png'.format(self.out_fldr, self.it), bbox_inches='tight')



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
