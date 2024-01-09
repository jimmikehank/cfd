import numpy as np
import os
from processing import *

mdot = 0.006
pulse_start = 0.2
pulse_duration = 0.0012
meanflow = 25
parallel = True
name = 'std_init_5'

mesh_command = 'python3 autoMesh.py --clean true --airfoil naca0015 --mdot {} --meanFlow {}'.format(mdot, meanflow)
blockmesh = 'blockMesh'

text_lines = ["boundaryField\n","{\n","\tinOut\n","\t\t{\n","\t\ttype            uniformFixedValue;\n","\t\tuniformValue\n","\t\t{\n","\t\t\ttype            constant;\n","\t\t\tvalue           ({} 0 0);\n".format(meanflow),"\t\t}\n","\t\tvalue           uniform ({} 0 0);\n".format(meanflow),"\t}\n","\tcoandaUpper\n","\t{\n","\t\ttype            flowRateInletVelocity;\n","\t\tmassFlowRate    \n","\t\t{\n","\t\t\ttype            scale;\n","\t\t\tscale\n","\t\t\t{\n","\t\t\t\ttype            constant;\n","\t\t\t\tvalue           {};\n".format(mdot),"\t\t\t}\n","\t\t\txScale          \n","\t\t\t{\n","\t\t\t\ttype            constant;\n","\t\t\t\tvalue           1;\n","\t\t\t}\n","\t\t\tvalue\n","\t\t\t{\n","\t\t\t\ttype            squarePulse;\n","\t\t\t\tstart           {};\n".format(pulse_start),"\t\t\t\tduration        {};\n".format(pulse_duration),"\t\t\t}\n","\t\t}\n","\t\tvalue           uniform (0 -0 -0);\n","\t}\n","\tcoandaLower\n","\t{\n","\t\ttype            noSlip;\n","\t}\n","\tsurface\n","\t{\n","\t\ttype            noSlip;\n","\t}\n","\tdefaultFaces\n","\t{\n","\ttype            empty;\n","\t}\n","}"]

def initializer(folder = './',text_write=''):
    import os
    import numpy as np
    os.system(mesh_command)
    os.system(blockmesh)
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



initializer(text_write = text_lines)

if parallel:
    decompose = 'decomposePar'
    reconstruct = 'reconstructPar'
    run_command = 'mpirun -np 4 rhoPimpleFoam -parallel'
    save_command = 'python3 autoMesh.py --store true --runName {} --airfoil naca0015'.format(name)
    os.system(decompose)
    try:
        os.system(run_command)
    except:
        input("System failed before completion, press ENTER to reconstruct.")
        os.system(reconstruct)
        exit
    os.system(reconstruct)
    # save_data = input("Save data (y/n): ").lower()
    # save_data = 'y'
    # if save_data == 'y':
    os.system(save_command)
    # else:
    #     exit()
    os.system('rhoPimpleFoam -postProcess -func forces -case /media/james/Data/james/completed_cases/coanda_airfoils/era/{}/'.format(name))
else:
    run_command = 'rhoPimpleFoam'
    os.system(run_command)
