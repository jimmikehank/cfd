import numpy as np
import os

aoa_sweep = np.linspace(-10,20,31)
print(aoa_sweep)

for i in aoa_sweep:
    mesh_command = "python3 autoMesh.py -clean -aoa {}".format(i)
    blockmesh = "blockMesh"
    decompose = "decomposePar"
    reconstruct = "reconstructPar -latestTime"
    run_command = "mpirun -np 6 rhoSimpleFoam -parallel"
    save_command = "python3 autoMesh.py -save {}".format(i)
    os.system(mesh_command)
    os.system(blockmesh)
    os.system(decompose)
    os.system(run_command)
    os.system(reconstruct)
    os.system(save_command)
