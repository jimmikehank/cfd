import numpy as np
import os


variable_selected = 'MDOT'
control_sweep = np.linspace(0, 0.008, 81)
input("{} control sweep: {}\n Press ENTER to continue, Ctrl-C to exit".format(variable_selected,control_sweep))
flow_speed = 25
parallel = False

testdir = '/home/james/Documents/research/completed_cases/coanda_airfoils/steady/NACA0015/downward/{}mps'.format(flow_speed)

if os.path.isdir(testdir):
    print('Directory exists')
    pass
else:
    print('Creating directory')
    os.system('mkdir {} -p'.format(testdir))



for i in control_sweep:
    if variable_selected.lower() == 'mdot':
        output = np.around(i,4)
    mesh_command = "python3 autoMesh.py --clean True --{} {} --airfoil NACA0015 --runName opt_base_sub_3".format(variable_selected.lower(),output)
    file_ext = str(output).replace('-','n')
    blockmesh = "blockMesh"
    if parallel:
        run_command = "mpirun -np 6 rhoSimpleFoam -parallel"
        decompose = "decomposePar"
        reconstruct = "reconstructPar -latestTime"
    else:
        run_command = "rhoSimpleFoam"
    save_command = "python3 autoMesh.py --store True --runName /steady/NACA0015/downward/{}mps/{}kgps".format(flow_speed,file_ext)
    os.system(mesh_command)
    os.system(blockmesh)
    if parallel:
        os.system(decompose)
        os.system(run_command)
        os.system(reconstruct)
    else:
        os.system(run_command)
    os.system(save_command)
