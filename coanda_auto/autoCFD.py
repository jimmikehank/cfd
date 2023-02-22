import os
import shutil
import numpy as np

mdot = 5.5

os.system('python3 autoMesh.py -clean')
os.system('blockMesh')
os.system('decomposePar')
os.system('mpirun -np 6 rhoSimpleFoam -parallel')
os.system('reconstructPar')
os.system('python3 autoMesh.py -save {}'.format(mdot))
