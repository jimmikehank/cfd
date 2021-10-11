import os
import shutil
import numpy as np

retain = ['0', 'constant', 'system','autoCFD.py']


boop = []

timestep = .00025;
start_folder = str(3*timestep)
print(start_folder)

N = 6
angles = np.linspace(9,15,N)

def check_float(textin):
    try:
        float(textin)
        return True

    except ValueError:
        return False

def rotrans(angle):
    command1 = 'surfaceTransformPoints -rollPitchYaw "(0 0 {})" constant/triSurface/wing_base.stl constant/triSurface/wing.stl'.format(angle)
    command2 = 'surfaceTransformPoints -translate "(-10 0 0)" constant/triSurface/wing.stl constant/triSurface/wing.stl'
    os.system(command1)
    os.system(command2)

def setup():
    os.system('blockMesh')
    os.system('snappyHexMesh')
    os.system('cp 0/* {}/'.format(start_folder))


for i in range(N):
    boop = []
    angle = angles[i]
    rotrans(angle)
    setup()

    os.system('simpleFoam')

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
    new_name = 'con_{}'.format(angle)
    os.rename(new_converged,new_name)

    retain = retain + [new_name]
    dirs = os.listdir()

    for item in dirs:
        if item not in retain:
            shutil.rmtree(item)
