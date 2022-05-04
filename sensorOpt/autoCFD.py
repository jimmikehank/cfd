import os
import shutil
import numpy as np

retain = ['0', 'constant', 'system','autoCFD.py','con','sensorOptimization.ipynb','.ipynb_checkpoints','autoWing.py','autoVTK.py','VTK', 'vtkTransfer.py']

boop = []

timestep = .00010;

max_angle = 16
N = 17
angles = np.linspace(0,max_angle,N)

def check_float(textin):
    try:
        float(textin)
        return True

    except ValueError:
        return False

def rotrans(angle):
    command1 = 'surfaceTransformPoints -rollPitchYaw "(0 0 {})" constant/triSurface/wing_base.stl constant/triSurface/wing.stl'.format(-angle)
    command2 = 'surfaceTransformPoints -translate "(-10 0 0)" constant/triSurface/wing.stl constant/triSurface/wing.stl'
    os.system(command1)
    os.system(command2)

def setup():
    os.system('blockMesh')
    os.system('decomposePar')
    os.system('snappyHexMesh -overwrite')
    os.system('extrudeMesh')
    os.system('deconstructPar')

def simpleclean(retain):
    dirs = os.listdir()
    for item in dirs:
        if item not in retain:
            shutil.rmtree(item)


def cleanhouse(copy,retain,target,angle):
    dirs = os.listdir()
    for item in dirs:
        if item not in retain:
            if check_float(item):
                if float(item)%1 == 0:
                    it = int(item)
                else:
                    it = float(item)
                boop = boop + [it]
            else:
                shutil.rmtree(item)
    new_converged = str(max(boop))
    new_name = '{}_deg'.format(angle)
    os.system('mkdir {}/{}'.format(target,new_name))
    os.system('mv {} {}/{}'.format(new_converged,target,new_name))

    for item in copy:
        os.system('cp {} {}/{} -r'.format(item,target,new_name))

    for item in dirs:
        if item not in retain:
            shutil.rmtree(item)



for i in range(N):
    boop = []
    angle = angles[i]
    copy = ['0','system','constant']
    target = '/home/james/Documents/research/converged_cases/sensorOptOutput'
    rotrans(angle)
    setup()

    os.system('mpirun -np 6 simpleFoam -parallel')

    cleanhouse(copy,retain,target,angle)
