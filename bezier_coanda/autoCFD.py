from processing import *
import numpy as np
import os
import time

airfoil = 'opt1'
rho = 1.17
mu  = 1.82e-5
c   = 0.3
b   = 0.156
cmu = 0.0075
re = 3e6
re_range = np.arange(0.25e6,3.0e6+1,250000)
cmu_range = np.arange(0.005,0.035,0.0025)
# cmu_range = [0]
re_range = [200000]
cmu_range= [0.02]

# continue_bool = input("Reynolds number sweep: {}\nCmu sweep: {}\nProceed (y/n): ".format(re_range,cmu_range)).lower()
# if continue_bool == 'y':
#     for cmu in cmu_range:
#         for re in re_range:
#             try:
#                 single_run(cmu, cmu, re, rho, mu, c, b, .005*cmu, urf=0.7, airfoil='opt1', parallel=False)
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
#


##################################
### For single modifiable runs ###
##################################

U = 15
cmu = 0.0
mdot = cmu_openloop(cmu, U, c, b, rho)
print('Mass flow selected: {}kg/s'.format(mdot))
time.sleep(1.0)
mesh = 'python3 autoMesh.py --clean True --airfoil {} --mdot {} --meanFlow {}'.format(airfoil, mdot, U)

os.system(mesh)
os.system('blockMesh')
os.system('rhoSimpleFoam')

f,m,tt = retrieve_lift('./')
cs = f[-1,:] / (0.5 * rho * U**2 * b * c)

print('Forces = {}'.format(cs))
