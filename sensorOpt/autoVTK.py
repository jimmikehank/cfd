import os
import shutil
import numpy as np

os.system('cd con')
dirs = os.listdir('con/')
print(dirs)


for item in dirs:
    print(item)

    os.system('foamToVTK -case con/{} -latestTime'.format(item))

print('VTK Conversion complete!')
