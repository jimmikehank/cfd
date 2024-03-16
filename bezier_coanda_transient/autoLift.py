import numpy as np
import os
from processing import *
from bezier_foil import *

casedir = '/media/james/Data/james/completed_cases/coanda_airfoils/era/era_std_2/'

forces, moments, time = retrieve_lift(casedir)

lift = forces[:,1]
drag = forces[:,0]

np.savetxt('./output/lift0.txt',lift)
np.savetxt('./output/drag0.txt',drag)
np.savetxt('./output/time0.txt',time)
