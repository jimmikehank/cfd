### NEED TO ADD DOUBLE ITERATION FOR FLOW ON VS FLOW OFF

from processing import *
from bezier_foil import *
import numpy as np
import os
import time
from multiprocessing import Process


mdots = np.array([0,0.003])
meanflow = 30
c = 0.3
starting_iteration = 0
max_iterations = 2000
run_name = "opt_fd5"
cp_size = 5
targs = np.arange(1,cp_size-1,1)
targs = np.append(targs, targs+cp_size)
direcs = [1,2]
delta = 1e-5
gamma = 2e-7
block = "blockMesh"
run = "rhoSimpleFoam"
force = "rhoSimpleFoam -postProcess -func forces"
foils_dir = '/home/james/Documents/research/cfd/airfoils/'
airfoil = 'naca0015'
vol_init = 0.00725
#
# window = np.ones([cp_size,2])
# window[0,0] = 0
# window[0,1] = 0
# window[1,0] = 0
# window[-1,0] = 0
# window[-1,1] = 0
# window = np.vstack([window,window])


for ii in range(starting_iteration, starting_iteration + max_iterations):
    dFdX_base = gradient_sweep(mdots[0],meanflow,run_name,cp_size,targs,ii,delta)
    dFdX_act  = gradient_sweep(mdots[1],meanflow,run_name,cp_size,targs,ii,delta)
    dFdX = dFdX_act - dFdX_base
    # dFdX = dFdX * window
    if ii == 0:
        print(dFdX)
        input("\n\nContinue?")
    np.savetxt('./output/gradients/{}_dFdX.txt'.format(ii),dFdX)
    # if ii == 0:
    #     input("Shape Gradient: {}\n\nContinue?".format(dFdX))
    cpu_init = np.loadtxt('{}/optimization/{}/{}_cpu.txt'.format(foils_dir,run_name,ii))
    cpl_init = np.loadtxt('{}/optimization/{}/{}_cpl.txt'.format(foils_dir,run_name,ii))
    cp = np.vstack([cpu_init,cpl_init])
    cp_new = cp + dFdX * gamma
    cpu_new = cp_new[:cp_size,:]
    cpl_new = cp_new[cp_size:,:]
    # cpu_new[1,0] = (((cpl_new[1,0]+c)/cpl_new[1,1]) * cpu_new[1,1]) - c
    cpu_new, cpl_new, dFdX, dVdX = constrain_volume(cpu_init, cpl_init, cpu_new, cpl_new, vol_init, constraint=0.35)
    cpu_new[0,0] = np.around(cpu_new[0,0],1)
    cpu_new[1,0] = np.around(cpu_new[1,0],1)
    cpl_new[0,0] = np.around(cpu_new[0,0],1)
    cpl_new[1,0] = np.around(cpu_new[1,0],1)
    np.savetxt('./output/gradients/{}_dVdX.txt'.format(ii),dFdX)
    np.savetxt('{}/optimization/{}/{}_cpu.txt'.format(foils_dir,run_name,ii+1),cpu_new)
    np.savetxt('{}/optimization/{}/{}_cpl.txt'.format(foils_dir,run_name,ii+1),cpl_new)
