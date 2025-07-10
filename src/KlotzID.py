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
    def __init__(self, pfile, pressure_var, volume_var, sim_folder, sim_times, 
                 lv_ed_pressure, lv_ed_volume, lv_v0, 
                 inflation_type, ncores, meshdir, datadir,
                 rv_ed_pressure=0.0, rv_ed_volume=0.0, rv_v0=0.0, 
                 constraint_vars=None, pfile_bv_init=None,pfile_bv=None,alternate_export=False,
                 plot_intermediate=False, save_intermediate=False,use_inflate_pfile=False):
        self.self_path = os.path.dirname(os.path.abspath(__file__))
        self.cheart_folder = os.path.dirname(pfile)
        if self.cheart_folder == '': self.cheart_folder = '.'
        self.ncores = ncores
        self.inflation_type = inflation_type
        assert (inflation_type == 'volume') or (inflation_type == 'inverse_volume') or \
            (inflation_type == 'volume_bivariable') or (inflation_type == 'volume_variable'), 'Non-recognized inflation type'

        self.pfile = os.path.basename(pfile)

        if inflation_type=='volume_bivariable':
            # Correct times to start from the second time point
            sim_times = (sim_times[0]+sim_times[2], sim_times[1], sim_times[2])

            

            if pfile_bv_init is None or pfile_bv is None:
                raise ValueError('For bivariable inflation, both bivariable P files must be provided')
            self.pfile_bv_init = pfile_bv_init
            self.pfile_bv = pfile_bv
            self.rv_ed_pressure = rv_ed_pressure
            self.rv_ed_volume = rv_ed_volume
            self.rv_v0 = rv_v0
            if constraint_vars is None:
                raise ValueError("For bivariable inflation, constraint variables must be provided")
            if len(constraint_vars) != 2:
                raise ValueError("For bivariable inflation, two constraint variables must be provided")
            self.pars=constraint_vars
            if rv_ed_pressure == 0.0 or rv_ed_volume == 0.0:
                raise ValueError("For bivariable inflation, ED pressure and volume for the RV must be provided")

        elif inflation_type=='volume_variable':
            # Correct times to start from the second time point
            sim_times = (sim_times[0]+sim_times[2], sim_times[1], sim_times[2])

            if pfile_bv_init is None or pfile_bv is None:
                raise ValueError('For variable inflation, both bivariable P files must be provided')
            self.pfile_bv_init = pfile_bv_init
            self.pfile_bv = pfile_bv
            self.rv_ed_pressure = rv_ed_pressure
            self.rv_ed_volume = rv_ed_volume
            self.rv_v0 = rv_v0
            if len(constraint_vars) != 1:
                raise ValueError("For variable inflation, one constraint variable must be provided")
            self.pars=constraint_vars
        else:
            pass


        # Output path
        self.sim_folder = sim_folder
        if not os.path.exists('{}'.format(sim_folder)):
            os.mkdir('{}'.format(sim_folder))
        self.meshdir = meshdir
        self.datadir = datadir

        self.pressure_var = pressure_var
        self.volume_var = volume_var
        self.alternate_export = alternate_export

        self.times = sim_times

        # Compute target curve
        self.lv_ed_pressure = lv_ed_pressure  # kPa
        self.lv_ed_volume = lv_ed_volume # mm3
        self.lv_v0 = lv_v0 # mm3
        self.klotz_volume, self.klotz_pressure = klotz_curve(lv_ed_volume, lv_ed_pressure)
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

        self.logfile = f'{self.sim_folder}/klotzid.log'
        self.inflate_file='inflate_values.P'

        self.use_inflate_pfile=use_inflate_pfile


    def prepare(self):
        # Creating folders
        if not os.path.exists('{}/tmp1'.format(self.sim_folder)):
            os.mkdir('{}/tmp1'.format(self.sim_folder))
        if not os.path.exists('{}/tmp2'.format(self.sim_folder)):
            os.mkdir('{}/tmp2'.format(self.sim_folder))
        if not os.path.exists('{}/tmp3'.format(self.sim_folder)):
            os.mkdir('{}/tmp3'.format(self.sim_folder))

        # Deleting log
        try:
            os.remove(self.logfile)
        except OSError:
            pass

        # Removing old inflate value file

        if self.use_inflate_pfile:
            try:
                os.remove(self.logfile)
            except OSError:
                pass

            # write values to new inflate value file
            inf_file=open(self.inflate_file,"w")
            inf_file.write("#LV_EDP="+str(self.lv_ed_pressure)+'\n')
            inf_file.write("#LV_EDV="+str(self.lv_ed_volume)+'\n')
            inf_file.write("#RV_EDP="+str(self.rv_ed_pressure)+'\n')
            inf_file.write("#RV_EDV="+str(self.rv_ed_volume)+'\n')
            inf_file.close()



    def post_clean(self):
        # Deleting temporaty folders.
        shutil.rmtree('{}/tmp1'.format(self.sim_folder))
        shutil.rmtree('{}/tmp2'.format(self.sim_folder))
        shutil.rmtree('{}/tmp3'.format(self.sim_folder))
        shutil.rmtree('{}/out'.format(self.sim_folder))
        os.remove(f'{self.sim_folder}/tmp1.log')
        os.remove(f'{self.sim_folder}/tmp2.log')
        os.remove(f'{self.sim_folder}/tmp3.log')


    def optimize(self, params, post_clean=False):
        self.prepare()

        # Checking number of parameters
        if self.inflation_type == 'volume_bivariable':
            assert len(params) == 4, 'For bivariable inflation, 4 parameters must be provided'
        elif self.inflation_type == 'volume_variable':
            assert len(params) == 3, 'For variable inflation, 3 parameters must be provided'
        else:
            assert len(params) == 2, 'For volume inflation, 2 parameters must be provided'

        # Initializing variables
        error = 1e3
        self.it = 0
        start_time = time.time()
        while (error > 1e-3) and (self.it < self.max_iterations):
            
            self.print_iteration(self.it)
            if self.inflation_type == 'volume_bivariable' or self.inflation_type == 'volume_variable':
                #print('Entering bivariable inflation')
                new_params_1= self.optimize_linear_variable(params)
                params=new_params_1

            new_params, error = self.optimize_non_linear(params)
            print('Finished iteration {:d} with update norm {:e}'.format(self.it, error))
            print('Parameters found are k={:f} and kb={:f}'.format(new_params[0], new_params[1]))

            self.write_log(params, error)

            # update
            params = new_params
            self.it += 1

            print()


        # For bvariable inflation I need to run one more inflation to account for LV-RV interaction
        if self.inflation_type=='volume_bivariable':
            k, kb, par_lv, par_rv = params
            print('Running variable inflation with optimized parameters k={:f} and kb={:f}'.format(k, kb))
            p_bv = self.run_cheart_variable_inflation(params,f'{self.sim_folder}/tmp3') #uses new k and kb
            p_bv.wait()

            times = (self.times[1], self.times[1], self.times[2])
            par_lv = chio.read_scalar_dfiles('{}/{}/{}'.format(self.sim_folder, 'tmp3', self.pars[0]), times)[-1]
            par_rv = chio.read_scalar_dfiles('{}/{}/{}'.format(self.sim_folder, 'tmp3', self.pars[1]), times)[-1]
            print('Variable inflation found par_lv={:f} and par_rv={:f}'.format(par_lv, par_rv))

            params = (k, kb, par_lv, par_rv)

        if self.it < self.max_iterations:
            end_time = time.time()
            print('Optimization succesful in {:d} iterations and {:2.3f} s'.format(self.it, end_time-start_time))
            print('Parameters found are k={:f} and kb={:f}'.format(params[0], params[1]))

            if self.inflation_type=='volume_bivariable':
                print('LV/RV params are par_lv={:f} and par_rv={:f}'.format(params[2], params[3]))
            elif self.inflation_type=='volume_variable':
                print('LV params are par_lv={:f}'.format(params[2]))

        else:
            print('Optimization failed.')

        if self.save_intermediate:
            int_vol = np.vstack(self.intermediate['vol'])
            int_pres = np.vstack(self.intermediate['pres'])

            try:
                int_params = np.vstack(np.array(self.intermediate['params']))
            except TypeError:
                print('Type error in saving intermediates!')
                int_params=np.array([0,0,0,0])

            np.savez('{}/intermediate_results.npz'.format(self.sim_folder), volume=int_vol, pressure=int_pres, params=int_params)

        print()

        if post_clean:
            self.post_clean()
        return params


    def optimize_non_linear(self, params):

        # Grabbing parameters
        if self.inflation_type=='volume_bivariable':
            k, kb, par_lv, par_rv = params
        elif self.inflation_type=='volume_variable':
            k, kb, par_lv = params
            par_rv=0
        else:
            k,kb = params
            par_lv=0
            par_rv=0

        # Update k if variable used
        if self.inflation_type=='volume_bivariable' or self.inflation_type=='volume_variable':
            print('Updating k using variable inflation par_LV')
            kold = k
            k = k*(1+par_lv)
            par_lv = 0
            if self.inflation_type=='volume_bivariable':
                par_rv_eff=kold*(1+par_rv)
                par_rv = par_rv_eff/k - 1

        print('Running simulation with parameters k='+str(k)+' and kb='+str(kb))
        if self.inflation_type=='volume_bivariable':
            print('Bivariable parameters par_lv='+str(par_lv)+' and par_rv='+str(par_rv))
        elif self.inflation_type=='volume_variable':
            print('Variable parameters par_lv='+str(par_lv))

        p1 = self.run_cheart_inflation((k, kb, par_lv, par_rv), f'{self.sim_folder}/tmp1')
        p2 = self.run_cheart_inflation((k, kb+self.eps, par_lv, par_rv), f'{self.sim_folder}/tmp2')
        exit_codes = [p.wait() for p in (p1, p2)]

        # Check if the simulations were successful
        check1 = check_simulation_log(f'{self.sim_folder}/tmp1.log')
        check2 = check_simulation_log(f'{self.sim_folder}/tmp2.log')
        if not check1 or not check2:
            raise RuntimeError('Simulation failed. Check the log files for more details.')

        # Load results
        pres = chio.read_scalar_dfiles('{}/{}/{}'.format(self.sim_folder, 'tmp1', self.pressure_var), self.times)
        vol = chio.read_scalar_dfiles('{}/{}/{}'.format(self.sim_folder, 'tmp1', self.volume_var), self.times)
        pres_eps = chio.read_scalar_dfiles('{}/{}/{}'.format(self.sim_folder, 'tmp2', self.pressure_var), self.times)

        if self.inflation_type == 'inverse_volume':
            vol = vol[::-1]
            vol = np.append(vol, self.lv_ed_volume)
            pres = np.append(0., pres)
            pres_eps = np.append(0., pres_eps)

        if self.plot_intermediate:
            self.plot_inflation_curve('{}/klotz_it{:d}.png'.format(self.sim_folder, self.it), vol, pres, (k, kb, par_lv, par_rv))

        if self.save_intermediate:
            self.intermediate['vol'].append(vol)
            self.intermediate['pres'].append(pres)
            self.intermediate['params'].append(params)

        # Linear correction
        k_delta = self.lv_ed_pressure/pres[-1]
        k = k*k_delta
        pres = pres*self.lv_ed_pressure/pres[-1]
        pres_eps = pres_eps*self.lv_ed_pressure/pres_eps[-1]

        # Levenber-marquadt iteration
        g = self.compute_curve_difference(vol, pres)
        kb_delta = self.LM_update(g, pres, pres_eps, kb)
        kb += kb_delta

        # Compute error
        error = np.abs(kb_delta) + np.abs(k_delta-1)

        # Return new parameters
        if self.inflation_type=='volume_bivariable':
            params = (k, kb, par_lv, par_rv)
        elif self.inflation_type=='volume_variable':
            params = (k, kb, par_lv)
        else:
            params = (k, kb)

        return params, error


    def optimize_linear_variable(self, params):
        if self.inflation_type == 'volume_variable':
            print('Running variable inflation')
            k, kb, par_lv = params
        elif self.inflation_type == 'volume_bivariable':
            print('Running bivariable inflation')
            k, kb, par_lv, par_rv = params
        else:
            raise ValueError('This function is only for variable inflation')

        p_bv = self.run_cheart_variable_inflation(params,f'{self.sim_folder}/tmp3') #uses new k and kb
        exit_code_bv = [p.wait() for p in ([p_bv])]

        times = (self.times[1], self.times[1], self.times[2])

        check = check_simulation_log(f'{self.sim_folder}/tmp3.log')
        if check:
            par_lv = chio.read_scalar_dfiles('{}/{}/{}'.format(self.sim_folder, 'tmp3', self.pars[0]), times)[-1]

            if self.inflation_type == 'volume_bivariable':
                par_rv = chio.read_scalar_dfiles('{}/{}/{}'.format(self.sim_folder, 'tmp3', self.pars[1]), times)[-1]
                params_full=(k, kb, par_lv, par_rv)
                print('Variable inflation found par_lv={:f} and par_rv={:f}'.format(par_lv, par_rv))
            else:
                params_full=(k, kb, par_lv)
                print('Variable inflation found par_lv={:f}'.format(par_lv))
        else:
            print('Variable inflation failed. Keeping previous parameters')
            if self.inflation_type == 'volume_bivariable':
                params_full=(k, kb, par_lv, par_rv)
            else:
                params_full=(k, kb, par_lv)



        return params_full



    def LM_update(self, g, pres, pres_eps, kb):
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

        # If np.abs(delta) > 1, divide by half
        while kb + delta < 0:
            print('Warning: Negative kb found. Dividing delta by 2')
            delta = delta/2


        return delta[0]


    def run_last_simulation(self, params, plot=True):
        # Grabbing parameters
        if self.inflation_type=='volume_bivariable':
            k, kb, par_lv, par_rv = params
        elif self.inflation_type=='volume_variable':
            k, kb, par_lv = params
            par_rv=0
        else:
            k,kb = params
            par_lv=0
            par_rv=0

        print('Running simulation with optimized parameters k={:f}, kb={:f}, par_lv={:f}, par_rv={:f}.'.format(k, kb, par_lv, par_rv))
        p1 = self.run_cheart_inflation((k, kb, par_lv, par_rv), f'{self.sim_folder}/out')
        p1.wait()

        # Load results
        pres = chio.read_scalar_dfiles('{}/{}'.format(self.sim_folder, self.pressure_var), self.times)
        vol = chio.read_scalar_dfiles('{}/{}'.format(self.sim_folder, self.volume_var), self.times)

        if self.inflation_type == 'inverse_volume':
            vol = vol[::-1]
            vol = np.append(vol, self.lv_ed_volume)
            pres = np.append(0., pres)

        ed_error = np.abs(self.lv_ed_pressure/pres[-1]-1)

        g = self.compute_curve_difference(vol, pres)
        curve_error = np.linalg.norm(g)

        if plot:
            self.plot_inflation_curve('{}/klotz_fit.png'.format(self.sim_folder), vol, pres, (k, kb, par_lv, par_rv))
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
            if self.inflation_type == 'volume_bivariable':

                klv=params[0]*(1+params[2])
                krv=params[0]*(1+params[3])
                kb=params[1]


                if self.alternate_export:
                    file.write("Iteration {:d}, k={:f}, kb={:f}, par_LV={:f}, par_RV={:f}\n".format(self.it,params[0],params[1],params[2],params[3]))

                else:
                    file.write("Iteration {:d}, k_lv={:f}, k_rv={:f}, kb={:f}, error = {:e}\n".format(self.it, klv,krv,kb, error))



            else:
                file.write("Iteration {:d}, k={:f}, kb={:f}, error = {:e}\n".format(self.it, params[0], params[1], error))

    def write_params(self,fname, params):
        with open(fname, "w") as file:

            if self.inflation_type == 'volume_bivariable':

                klv=params[0]*(1+params[2])
                krv=params[0]*(1+params[3])
                kb=params[1]

                if self.alternate_export:
                # Writing data to a file
                    file.write("#k={:12.12f}\n".format(params[0]))
                    file.write("#kb={:12.12f}\n".format(params[1]))
                    file.write("#par_LV={:12.12f}\n".format(params[2]))
                    file.write("#par_RV={:12.12f}".format(params[3]))

                else:
                    file.write("#k_lv={:12.12f}\n".format(klv))
                    file.write("#k_rv={:12.12f}\n".format(krv))
                    file.write("#kb={:12.12f}".format(kb))

            elif self.inflation_type == 'volume_variable':
                klv=params[0]*(1+params[2])
                kb=params[1]

                if self.alternate_export:
                # Writing data to a file
                    file.write("#k={:12.12f}\n".format(params[0]))
                    file.write("#kb={:12.12f}\n".format(params[1]))
                    file.write("#par_LV={:12.12f}".format(params[2]))

                else:
                    file.write("#k={:12.12f}\n".format(klv))
                    file.write("#kb={:12.12f}".format(kb))

            else:

                file.write("#k={:12.12f}\n".format(params[0]))
                file.write("#kb={:12.12f}".format(params[1]))


    def run_cheart_inflation(self, params, outdir):
        k, kb, par_lv, par_rv = params

        # Run cheart
        with open('{}.log'.format(outdir), 'w') as ofile:
            p = Popen(['bash', '{}/run_inflation.sh'.format(self.self_path), 
                       '{:f}'.format(k),                    # 1    
                       '{:f}'.format(kb),                   # 2
                       '{:f}'.format(par_lv),               # 3
                       '{:f}'.format(par_rv),               # 4
                       outdir,                              # 5   
                       '{:d}'.format(self.ncores),          # 6
                       self.cheart_folder,                  # 7
                       self.pfile,                          # 8  
                       '{:f}'.format(self.lv_ed_volume),    # 9
                       '{:f}'.format(self.rv_ed_volume),    # 10
                       '{:f}'.format(self.lv_v0),           # 11
                       '{:f}'.format(self.rv_v0),           # 12
                       self.meshdir, self.datadir,],        # 13, 14
                       stdout=ofile, stderr=ofile)

        return p


    def run_cheart_variable_inflation(self, params, outdir):
        k, kb, par_lv, par_rv = params
        
        # Run cheart
        with open('{}.log'.format(outdir), 'w') as ofile:
            p = Popen(['bash', '{}/run_variable_inflation.sh'.format(self.self_path), 
                      '{:f}'.format(k),                     # 1
                      '{:f}'.format(kb),                    # 2 
                        '{:f}'.format(par_lv),              # 3
                        '{:f}'.format(par_rv),              # 4 
                       outdir,                              # 5
                       '{:d}'.format(2*self.ncores),        # 6
                       self.cheart_folder,                  # 7
                       self.pfile_bv_init,                  # 8
                       self.pfile_bv,                       # 9
                       '{:f}'.format(self.lv_ed_volume),    #10 
                       '{:f}'.format(self.lv_ed_pressure),  #11
                       '{:f}'.format(self.rv_ed_volume),    #12
                       '{:f}'.format(self.rv_ed_pressure),  #13
                       '{:f}'.format(self.lv_v0),           #14
                       '{:f}'.format(self.rv_v0),       #15
                       self.meshdir, self.datadir,],    # 16, 17
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


    def print_iteration(self, it):
        print()
        print('-------------------------------------------')
        print('ITERATION {:d}'.format(self.it))



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


def check_simulation_log(log_file):
    """
    Check the simulation log file for errors.
    Returns True if no errors are found, False otherwise.
    """
    if not os.path.exists(log_file):
        raise FileNotFoundError(f"Log file {log_file} does not exist.")

    with open(log_file, 'r') as f:
        lines = f.readlines()
    
    for line in lines[-3:]:
        if 'Program complete ...' in line:
            print("Simulation completed successfully.")
            return True
    
    return False