import os
import shutil

retain = ['0', 'constant', 'system','autoCFD.py']
dirs = os.listdir()

boop = []

def check_float(textin):
    try:
        float(textin)
        return True

    except ValueError:
        return False

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

retain = retain + [str(max(boop))]
dirs = os.listdir()

print(retain)
for item in dirs:
    if item not in retain:
        shutil.rmtree(item)
