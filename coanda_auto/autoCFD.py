import os
import shutil
import numpy as np

mdot = 4.0

os.system('python3 autoMesh.py -clean')
os.system('blockMesh')
os.system('decomposePar')
os.system('mpirun -np 6 rhoSimpleFoam -parallel')
os.system('reconstructPar -latestTime')
os.system('python3 autoMesh.py -save {}'.format(mdot))
