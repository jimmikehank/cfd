import numpy as np
import os

freq_sweep = np.array([5,6,7,8,9,10,15,20])
flow_speed = 25
print(freq_sweep)

for i in freq_sweep:
    mesh_command = "python3 autoMesh.py -clean -freq {}".format(i)
    blockmesh = "blockMesh"
    decompose = "decomposePar"
    reconstruct = "reconstructPar"
    run_command = "mpirun -np 6 rhoPimpleFoam -parallel"
    save_command = "python3 autoMesh.py -save /{}/{}hz_{}mps".format(flow_speed,i,flow_speed)
    os.system(mesh_command)
    os.system(blockmesh)
    os.system(decompose)
    os.system(run_command)
    os.system(reconstruct)
    os.system(save_command)
