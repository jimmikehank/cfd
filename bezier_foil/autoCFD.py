import numpy as np
import os
import sys

aoa_sweep = np.linspace(11,13,3)
temp_sweep = np.linspace(300,450,16)
args = sys.argv

blockmesh = "blockMesh"
decompose = "decomposePar"
reconstruct = "reconstructPar -latestTime"
run_command = "mpirun -np 6 rhoSimpleFoam -parallel"


if len(args) > 1:
    if args[1] == '-single':
        mesh_command = 'python3 autoMesh.py -clean'
        os.system(mesh_command)
        os.system(blockmesh)
        os.system(decompose)
        os.system(run_command)
        os.system(reconstruct)

else:
    print("Temp sweep values: {}\n AoA sweep values: {}".format(temp_sweep, aoa_sweep))
    runval = input("Run full sweep (y/n): ")
    print(runval)
    if runval == 'y':
        for j in aoa_sweep:
            for i in temp_sweep:
                mesh_command = "python3 autoMesh.py -clean -temp {} -aoa {}".format(i,j)
                save_command = "python3 autoMesh.py -save {}/{}".format(j,i)
                os.system(mesh_command)
                os.system(blockmesh)
                os.system(decompose)
                os.system(run_command)
                os.system(reconstruct)
                os.system(save_command)
    else:
        pass
