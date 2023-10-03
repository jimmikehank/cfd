import numpy as np
import os
import sys

aoa_sweep = np.linspace(-12,12,25)
parallel = False
flow_speed = 50

blockmesh = "blockMesh"
run_command = "rhoSimpleFoam"
save_folder = "/home/james/Documents/research/completed_cases"
testdir = '/home/james/Documents/research/completed_cases/pitch_airfoils/steady/NACA0015/{}mps'.format(flow_speed)

if os.path.isdir(testdir):
    print('Directory exists')
    pass
else:
    print('Creating directory')
    os.system('mkdir {} -p'.format(testdir))


print("AoA sweep values: {}".format(aoa_sweep))
runval = input("Run full sweep (y/n): ")

if runval == 'y':
    for j in aoa_sweep:
        file_ext = str(int(np.around(j,4))).replace('-','n')
        mesh_command = "python3 autoMesh.py --clean true --aoa {}".format(j)
        save_command = "python3 autoMesh.py --store True --runName {}/{}".format(testdir,file_ext)
        os.system(mesh_command)
        os.system(blockmesh)
        os.system(run_command)
        os.system(save_command)
