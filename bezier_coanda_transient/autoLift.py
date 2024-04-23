import numpy as np
import os
from processing import *
from bezier_foil import *
from matplotlib import pyplot as plt


casedir = '/media/james/Data/james/completed_cases/coanda_airfoils/era/imp'
targets = ['./']
x = input("Target folders: {}".format(targets))


for target in targets:
    forces, moments, time = retrieve_lift(target)
    l = forces[1:,1]
    d = forces[1:,1]
    t = time
    plt.figure()
    plt.plot(t,l)
    plt.savefig('./output/test.png')
    np.savetxt('./output/limp.txt',l)
    np.savetxt('./output/timp.txt',t)
    
