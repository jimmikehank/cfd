import os
import shutil
import numpy as np

retain = ['0', 'constant', 'system','autoCFD.py','converged','test_bench','autoBC.py','autoMesh.py','images']
copy = ['0','system','constant']
boop = []

N = 21
start = 100000
step_size = 5000

pressures = np.linspace(start,(N-1)*step_size+start,N)

input('Your pressure steps are: {} \n Press ENTER to confirm'.format(pressures))

def check_float(textin):
    try:
        float(textin)
        return True

    except ValueError:
        return False

for i in range(N):
    boop = []
    pressure = pressures[i]
    os.system('python3 autoBC.py {}'.format(pressure))
    os.system('decomposePar')
    os.system('mpirun -np 4 rhoSimpleFoam -parallel')
    os.system('reconstructPar -latestTime')

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

    if pressure-100000 < 100000:
        new_con = '0{}'.format(pressure-100000)
    else:
        new_con = '{}'.format(pressure-100000)

    foldname = 'afc_{}kPag'.format(new_con[0:3])
    os.system('mkdir converged/{}'.format(foldname))
    os.system('mv {} converged/{}'.format(new_converged,foldname))

    for item in copy:
        os.system('cp {} converged/{} -r'.format(item,foldname))

    dirs = os.listdir()

    for item in dirs:
        if item not in retain:
            shutil.rmtree(item)
