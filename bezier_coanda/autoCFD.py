import numpy as np
import os

mdot_sweep = np.array([0,.0005, .001, .0015, 0.002,0.0025,0.003,0.0035,0.004,0.0045,0.005,0.0055,0.006])
flow_speed = 25
print(mdot_sweep)

for i in mdot_sweep:
    mesh_command = "python3 autoMesh.py -clean -mdot {}".format(i)
    blockmesh = "blockMesh"
    decompose = "decomposePar"
    reconstruct = "reconstructPar -latestTime"
    run_command = "mpirun -np 6 rhoSimpleFoam -parallel"
    save_command = "python3 autoMesh.py -save /steady/RAE2822/upward/{}kgps_{}mps".format(i,flow_speed)
    os.system(mesh_command)
    os.system(blockmesh)
    os.system(decompose)
    os.system(run_command)
    os.system(reconstruct)
    os.system(save_command)
