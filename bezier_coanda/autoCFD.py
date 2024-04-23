from processing import *
import numpy as np
import os
import time

airfoil = 'naca0015'
rho = 1.17
mu  = 1.82e-5
c   = 0.3
b   = 0.156
cmu = 0.0075
re = 3e6
# re_range = np.arange(0.25e6,1.75e6+1,250000)
# cmu_range = np.arange(0.005,0.035,0.0025)
# # cmu_range = [0]
# # re_range = [200000]
# cmu_range= [0]
#
# continue_bool = input("Reynolds number sweep: {}\nCmu sweep: {}\nProceed (y/n): ".format(re_range,cmu_range)).lower()
# if continue_bool == 'y':
#     for cmu in cmu_range:
#         for re in re_range:
#             try:
#                 single_run(cmu, cmu, re, rho, mu, c, b, .005*cmu, urf=0.7, airfoil='naca0015', parallel=False)
#             except(KeyboardInterrupt):
#                 print("Run aborted!")
#                 exit()
#             except:
#                 pass
#
# elif continue_bool == 'n':
#     print("\nSimulation exited\n")
# else:
#     print("Input must be y or n, simulation exited")
# #


##################################
### For single modifiable runs ###
##################################

aoas = np.arange(0,10.1,2.5)
aoas = [0]
mdot_range = np.arange(0.,0.0031,0.0005)
# mdot_range = [0,0.003]
re_range = [15,20,25]
for re in re_range:
    meanflow = re
    for aoa in aoas:
        for i in mdot_range:
            mesh = 'python3 autoMesh.py --clean True --airfoil naca0015 --mdot {} --meanFlow {} --aoa {}'.format(i, meanflow, aoa)
            block = 'blockMesh'
            run = 'rhoSimpleFoam'
            save = 'python3 autoMesh.py --store true --runName /u_{}/mdot_{}'.format(int(re),np.around(i,5))

            try:
                os.system(mesh)
                os.system(block)
                os.system(run)
                os.system(save)
            except(KeyboardInterrupt):
                exit()
