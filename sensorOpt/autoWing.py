import os
import shutil
import numpy as np


def rotrans(angle):
    command1 = 'surfaceTransformPoints -rollPitchYaw "(0 0 {})" constant/triSurface/wing_base.stl constant/triSurface/wing.stl'.format(angle)
    command2 = 'surfaceTransformPoints -translate "(-10 0 0)" constant/triSurface/wing.stl constant/triSurface/wing.stl'
    os.system(command1)
    os.system(command2)

def setup():
    os.system('blockMesh')
    os.system('snappyHexMesh -overwrite')
    os.system('extrudeMesh')

rotrans(-4.5)
setup()

os.system('simpleFoam')
