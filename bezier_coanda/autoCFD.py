from processing import *
import numpy as np
import os
import time


# rho = 1.17
# mu  = 1.82e-5
# c   = 0.6
# b   = 0.156
# cmu = 0.0075
# re = 3e6
# re_range = np.arange(0.25e6,3.0e6+1,250000)
# cmu_range = np.arange(0.005,0.035,0.0025)
# # cmu_range = [0]
#
# continue_bool = input("Reynolds number sweep: {}\nCmu sweep: {}\nProceed (y/n): ".format(re_range,cmu_range)).lower()
# if continue_bool == 'y':
#     for cmu in cmu_range:
#         for re in re_range:
#             try:
#                 single_run(cmu, cmu, re, rho, mu, c, b, .005*cmu, urf=0.7, parallel=False)
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



##################################
### For single modifiable runs ###
##################################

meanflow = 30
aoas = np.arange(0,14.1,2)
mdot = 0.0075
mdotu = np.linspace(0,mdot,21)
mdotl = mdot - mdotu
# mdots = np.arange(.0005,0.0081,0.0005)
rho = 1.17
c = 0.3
b = 0.156
M = len(mdotu)
N = len(aoas)
print("AoA sweep: {}\n".format(aoas))
conty = input("Commencing {} runs, continue (y/n): ".format(M*N))

# if os.path.exists('./output/mdot_to_Ujet.txt'):
#     print("File exists")
#     pass
# else:
#     print("File does not exist, creating data file")
#     np.savetxt('./output/mdot_to_Ujet.txt',np.array([0,0]))

if conty.lower() == 'n':
    exit()

for j in range(N):
    aoa = aoas[j]
    for i in range(M):
        mesh = 'python3 autoMesh.py --clean true --airfoil naca0015 --mdot {} --mdot_lower {} --meanFlow {} --aoa {}'.format(mdotu[i],mdotl[i],meanflow,aoa)
        block = 'blockMesh'
        run = 'rhoSimpleFoam'

        try:
            os.system(mesh)
            os.system(block)
            os.system(run)
            forces,moments,time = retrieve_lift('./')
            U = max_velocity('./')
            T = U * mdot
            QS = 0.5 * 1.17 * meanflow**2 * b * c
            drag = forces[-1,0] / QS
            lift = forces[-1,1] / QS
            point = np.array([mdot,U])
            print("Max velocity: {}\nThrust AFC: {}\nDrag: {}\nLift: {}".format(U,T,drag-T,lift))
            print(0.5 * U**2 * 0.009*.0254 * 0.156 *.0254 * 1.17)
            save = 'python3 autoMesh.py --store True --runName /aoa{}/u{}_l{}'.format(aoa,mdotu[i],mdotl[i])
            os.system(save)

        except:
            exit()
