import os
import shutil
import numpy as np

def rotrans(angle):
    # Note, Global origin located at leading edge of wing_base.stl
    # Therefore: Leading edge located at <-10 0 0>
    command1 = 'surfaceTransformPoints -rollPitchYaw "(0 0 {})" constant/triSurface/wing_base.stl constant/triSurface/wing.stl'.format(angle)
    command2 = 'surfaceTransformPoints -translate "(-10 0 0)" constant/triSurface/wing.stl constant/triSurface/wing.stl'
    os.system(command1)
    os.system(command2)

def ofSetup():
    os.system('blockMesh')
    os.system('snappyHexMesh -overwrite')
    os.system('extrudeMesh')

def check_float(textin):
    try:
        float(textin)
        return True

    except ValueError:
        return False

def closeout(angle):

    os.system('mkdir con/{}_deg_turb'.format(angle))
    dirs = os.listdir()
    copy = ['0','constant','system']

    boop = []
    for item in dirs:
        if check_float(item):
            if float(item)%1 == 0:
                it = int(item)
            else:
                it = float(item)
            boop = boop + [it]

    converged_soln = str(max(boop))
    if converged_soln != '0':
        os.system('mv {} con/{}_deg'.format(converged_soln,angle))
    for item in copy:
        os.system('cp {} con/{}_deg -r'.format(item,angle))

angle = 15
rotrans(-angle)
ofSetup()
os.system('simpleFoam')
closeout(angle)
