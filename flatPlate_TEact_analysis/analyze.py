import os

os.system('mapFields ~/Documents/research/cfd/flatPlate_TEact')

os.system("simpleFoam -postProcess -func 'forces'")
