import numpy as np
import os
from processing import *

U = 25
c = 0.3
tau = c/U
mdot = 0.004
ks = np.arange(np.pi,np.pi*2+.1,np.pi)

# print("Frequency range to run: {}".format(ks))

def get_lines_freq(U,mdot,freq):
    with open('./system/inits/sine','r') as f:
        lines = f.readlines()
        f.close()

    lines[8] = lines[8].format(U)
    lines[10] = lines[10].format(U)
    lines[21] = lines[21].format(mdot/2)
    lines[23] = lines[23].format(freq)
    lines[28] = lines[28].format(mdot/2)
    lines[31] = lines[31].format(mdot)

    lines[42] = lines[42].format(-mdot/2)
    lines[44] = lines[44].format(freq)
    lines[49] = lines[49].format(mdot/2)
    lines[52] = lines[52].format(mdot)

    return lines

def get_lines_imp(U,mdot,start,dur):
    with open('./system/inits/imp','r') as f:
        lines = f.readlines()
        f.close()

    lines[8] = lines[8].format(U)
    lines[10] = lines[10].format(U)
    lines[21] = lines[21].format(mdot)
    lines[31] = lines[31].format(start)
    lines[32] = lines[32].format(dur)

    return lines






def initializer(mesh,block,folder = './',text_write=''):
    import os
    import numpy as np
    os.system(mesh)
    os.system(block)
    init_dir = '/home/james/Documents/research/cfd/bezier_coanda_init/'
    target_dir = '/home/james/Documents/research/cfd/bezier_coanda_transient/0.0001'
    files = os.listdir(init_dir)
    cleaned = [int(x) for x in files if check_float(x)]
    print(cleaned)
    ic_file = np.max(cleaned)
    os.system('cp -r /home/james/Documents/research/cfd/bezier_coanda_init/{} {}'.format(ic_file, target_dir))
    print('Attempting to open')
    with open('./0.0001/U','r') as f:
        print('Did open')
        lines = f.readlines()
        sel = [i for i in range(len(lines)) if lines[i] == "boundaryField\n"][0]
        print(sel)
        newlines = lines[:sel]
        f.close()
    with open('./0.0001/U','w') as f:
        f.writelines(newlines)
        f.writelines(text_write)
        f.close()
    os.remove('./0.0001/uniform/time')
    os.system('cp ./system/inits/time ./0.0001/uniform/time')



for k in ks:
    f = k * U / 0.3 / np.pi
    mesh = "python3 autoMesh.py --clean true --mdot {} --frequency {} --aoa 0 --airfoil naca0015 --meanFlow {}".format(mdot,f,U)
    block = "blockMesh"
    decomp = "decomposePar"
    run = "mpirun -np 4 rhoPimpleFoam -parallel"
    recon = "reconstructPar"
    forces = "rhoPimpleFoam -postProcess -func forces"
    save = "python3 autoMesh.py --store true --runName /freq3/k_{}".format(k)

    f = k / np.pi / tau
    text_lines = get_lines_freq(U,mdot,f)

    initializer(mesh,block,text_write = text_lines)
    os.system(decomp)
    os.system(run)
    os.system(recon)
    os.system(forces)
    os.system(save)








# mdot = 0.003
# pulse_start = 0.05
# pulse_duration = 0.0005
# meanflow = 20
# parallel = True
# name = 'imp'
#
#
# mesh_command = 'python3 autoMesh.py --clean true --airfoil naca0015 --mdot {} --meanFlow {}'.format(mdot, meanflow)
# blockmesh = 'blockMesh'
#
# text_lines = get_lines_imp(meanflow,mdot,pulse_start,pulse_duration)
#
# initializer(mesh_command,blockmesh,folder = './', text_write = text_lines)
#
# if parallel:
#     decompose = 'decomposePar'
#     reconstruct = 'reconstructPar'
#     run_command = 'mpirun -np 4 rhoPimpleFoam -parallel'
#     force_command = 'rhoSimpleFoam -postProcess -func forces'
#     save_command = 'python3 autoMesh.py --store true --runName {} --airfoil naca0015'.format(name)
#     os.system(decompose)
#     try:
#         os.system(run_command)
#     except:
#         input("System failed before completion, press ENTER to reconstruct.")
#         os.system(reconstruct)
#         exit()
#     os.system(reconstruct)
#     os.system(force_command)
#     # save_data = input("Save data (y/n): ").lower()
#     # save_data = 'y'
#     # if save_data == 'y':
#     os.system(save_command)
#     # else:
#     #     exit()
#     # os.system('rhoPimpleFoam -postProcess -func forces -case /media/james/Data/james/completed_cases/coanda_airfoils/era/{}/'.format(name))
# else:
#     run_command = 'rhoPimpleFoam'
#     os.system(run_command)
