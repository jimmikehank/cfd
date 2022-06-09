import sys
import numpy as np

folder = '/Documents/research/cfd/coanda_auto/forces/0.010582'

G = np.load('forces.dat'.format(folder),allow_pickle=True)

print('\n Test end')
