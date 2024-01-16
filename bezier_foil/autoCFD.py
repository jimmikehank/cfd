import numpy as np
import os
import sys
from processing import *

# aoa_sweep = np.linspace(-12,12,25)
# parallel = False
# flow_speed = 50
#
# blockmesh = "blockMesh"
# run_command = "rhoSimpleFoam"
# save_folder = "/home/james/Documents/research/completed_cases"
# testdir = '/home/james/Documents/research/completed_cases/pitch_airfoils/steady/NACA0015/{}mps'.format(flow_speed)
#
# if os.path.isdir(testdir):
#     print('Directory exists')
#     pass
# else:
#     print('Creating directory')
#     os.system('mkdir {} -p'.format(testdir))
#
#
# print("AoA sweep values: {}".format(aoa_sweep))
# runval = input("Run full sweep (y/n): ")
#
# if runval == 'y':
#     for j in aoa_sweep:
#         file_ext = str(int(np.around(j,4))).replace('-','n')
#         mesh_command = "python3 autoMesh.py --clean true --aoa {}".format(j)
#         save_command = "python3 autoMesh.py --store True --runName {}/{}".format(testdir,file_ext)
#         os.system(mesh_command)
#         os.system(blockmesh)
#         os.system(run_command)
#         os.system(save_command)

meanflow =120
aoa = 0
alpha = aoa * np.pi / 180

mesh = "python3 autoMesh.py --clean true --airfoil naca0006 --meanFlow {} --aoa {}".format(meanflow, aoa)
block = "blockMesh"
run = "rhoSimpleFoam"
force = "rhoSimpleFoam -postProcess -func forces"

os.system(mesh)
os.system(block)
os.system(run)
os.system(force)

c = 1.0
b = 1.0
QS = 0.5 * 1.17 * meanflow**2 * c * b
forces, moments, time = retrieve_lift('./')
L = forces[-1,1] * np.cos(alpha) - forces[-1,0] * np.sin(alpha)
D = forces[-1,0] * np.cos(alpha) + forces[-1,1] * np.sin(alpha)
My = moments[-1,2]

CD = D / QS
CL = L / QS
Cm = My / (QS*c)

print("CD: {} \nCL: {} \nCm: {}\n".format(D/QS, L/QS, My/(QS*c)))

print("CL / AoA: {}".format(CL / (alpha)))
