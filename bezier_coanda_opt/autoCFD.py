from processing import *
import numpy as np
import os
import time


# # cmu = 0.015
# rho = 1.17
# mu  = 1.82e-5
# c   = 0.5
# b   = 0.156
# cmu = 0.0125
#
# # re = 3e6
# re_range = np.arange(.25e6,4.5e6+1,250000)
# # cmu_range = np.arange(0.001,0.0031,0.0005)/2
#
# continue_bool = input("Reynolds number sweep: {}\n Proceed (y/n): ".format(re_range)).lower()
# if continue_bool == 'y':
#     for re in re_range:
#         single_run(cmu, cmu, re, rho, mu, c, b, .005*cmu, urf=0.7, parallel=True)
#
# # continue_bool = input("Cmu sweep: {}\n Proceed (y/n): ".format(cmu_range)).lower()
# # if continue_bool == 'y':
# #     for cmu in cmu_range:
# #         single_run(cmu, cmu, re, rho, mu, c, b, .005*cmu, urf=0.5, parallel=True)
#
# elif continue_bool == 'n':
#     print("\nSimulation exited\n")
# else:
#     print("Input must be y or n, simulation exited")

### For single modifiable runs

mdot = 0.007
meanflow = 60


mesh = "python3 autoMesh.py --clean true --airfoil naca0015 --mdot {} --meanFlow {}".format(mdot, meanflow)
block = "blockMesh"
run = "rhoSimpleFoam"
force = "rhoSimpleFoam -postProcess -func forces"



os.system(mesh)
os.system(block)
os.system(run)
os.system(force)

U = max_velocity()
T = U * mdot
cmu = T / (0.5 * 1.17 * meanflow**2 * .3 * .156)
c = 0.3
b = 0.156
QS = 0.5 * 1.17 * meanflow**2 * c * b
forces, moments, time = retrieve_lift('./')
L = forces[-1,1]
D = forces[-1,0]
My = moments[-1,2]
print("CD: {} \nCL: {} \nCm: {}\n".format(D/QS, L/QS, My/(QS*c)))

clcmu = L / T
print("dCL/dCmu = {}".format(clcmu))
