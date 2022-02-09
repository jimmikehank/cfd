import numpy as np
import os

os.system('rhoSimpleFoam -postProcess -func forces -latestTime')

forceFold = os.listdir('./postProcessing/forces')
forceFold = forceFold[0]
print(forceFold)
data = np.genfromtxt(
    './postProcessing/forces/{}/forces.dat'.format(forceFold),
    skip_header=1,
    skip_footer=1,
    names=True,
    dtype=None,
    delimiter=' '
)

print(type(data[0]))
