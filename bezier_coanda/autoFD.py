### NEED TO ADD DOUBLE ITERATION FOR FLOW ON VS FLOW OFF

from processing import *
import numpy as np
import os
import time

mdot = 0.003
meanflow = 30
max_iterations = 2
run_name = "opt_fd1"
cp_size = 8
targs = np.array([1,2,3,4,5,6])
targs = np.append(targs, targs+cp_size)
direcs = [1,2]
delta = 1e-6
gamma = 1e-10
print(targs)

block = "blockMesh"
run = "rhoSimpleFoam"
force = "rhoSimpleFoam -postProcess -func forces"
foils_dir = '/home/james/Documents/research/cfd/airfoils/'
airfoil = 'naca0015'

for ii in range(max_iterations):
    dFdX = np.zeros(2*cp_size,2)
    if ii == 0:
        mesh = "python3 autoMesh.py --clean true --airfoil naca0015 --mdot {} --meanFlow {} --runName {}".format(mdot, meanflow, run_name)
        os.system(mesh)
        os.system(block)
        os.system(run)
        os.system(force)
        forces,moments,tt = retrieve_lift('./')
        L_base = forces[-1,1]

        cpu = np.loadtxt(foils_dir + 'control_points/{}_cpu.txt')
        cpl = np.loadtxt(foils_dir + 'control_points/{}_cpl.txt')
        points = np.vstack([cpu,cpl])
        for j in targs:
            for k in direcs:
                mesh = "python3 autoMesh.py --clean true --airfoil naca0015 --mdot {} --meanFlow {} --runName {} --FD {} {}".format(mdot, meanflow, run_name, j, k)
                os.system(mesh)
                os.system(block)
                os.system(run)
                os.system(force)
                forces,moments,tt = retrieve_lift('./')
                dFdX[j,k-1] = (forces[-1,1]-L_base)/delta
        print(dFdX)
        input("Do you want to continue?")

    else:
