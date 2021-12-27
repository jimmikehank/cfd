import os
import shutil
import numpy as np

os.system('cd con')
dirs = os.listdir('con/')
print(dirs)


for item in dirs:
    print(item)
    dir1 = os.listdir('con/{}/VTK/'.format(item))
    dir2 = os.listdir('con/{}/VTK/wing/'.format(item))
    for sub in dir1:
        if '.vtk' in sub:
            os.system('mv con/{}/VTK/{} VTK/{}.vtk'.format(item,sub,item))
        else:
            continue
    os.system('mv con/{}/VTK/wing/{} VTK/wing/{}.vtk'.format(item,dir2[0],item))


print('VTK Conversion complete!')
