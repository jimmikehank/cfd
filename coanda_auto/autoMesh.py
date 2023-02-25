import numpy as np
import matplotlib.pyplot as plt
import sys

# List of files / folders to protect from deletion
retain = ['0','constant','system','test_bench','gradOpt.py','autoCFD.py','autoBCu.py','autoBC.py','autoMesh.py']


# Argument handling functions:

def cleanup(retain):
    print(retain)
    import shutil
    import os
    dirs = os.listdir()
    for item in dirs:
        if item not in retain:
            shutil.rmtree(item)
        else:
            continue

def store(retain, target):
    import os
    import shutil
    ignore = ['data_process.ipynb','autoMesh.py','test_bench','gradOpt.py','autoCFD.py','autoBCu.py','autoBC.py']
    copy = ['0','system','constant','dynamicCode']
    dirs = os.listdir()
    delete = []
    for item in dirs:
        if item[0:3] == 'pro':
            delete.append(item)
    casefile = "/home/james/Documents/research/completed_cases/coanda_plenum/{}/".format(target)
    if os.path.exists(casefile):
        existing = os.listdir(casefile)
    else:
        os.mkdir(casefile)
        existing = []
    for item in dirs:
        if item in delete:
            shutil.rmtree(item)
        elif item not in retain and item not in ignore and item not in existing:
            shutil.move(item,casefile)
        elif item in copy:
            shutil.copytree(item,casefile+item)
#---------- Control Variables handled by this block ----------#

# Command line argument handler:

args = sys.argv

# This function parses input arguments to generate mesh allowing external programs to iterate through different mesh configurations

def arg_handle(args):
    # Default Values
    Rc = 0.0036   # Coanda cylinder radius
    te = 0.000178  # R/te = 20
    tu = 0.000229  # Upper surface thickness at exit
    ru = 0.5    # Contour radius of upper surface to exit
    ai = 60.0   # Contraction angle inside plenum

    for i in range(np.size(args)):
        current = args[i].lower()
        if current == '-rc':
            Rc = float(args[i+1])
        elif current == '-te':
            te = float(args[i+1])
        elif current == '-tu':
            tu = float(args[i+1])
        elif current == '-ru':
            ru = float(args[i+1])
        elif current == '-ai':
            ai = float(args[i+1])
        elif current == '-clean':
            cleanup(retain)
        elif current == '-save':
            target = args[i+1]
            try:
                store(retain,target)
            except:
                "Invalid target!"

    return Rc, te, tu, ru, ai

Rc, te, tu, ru, ai = arg_handle(args)

#-------------------------------------------------------------#

print('\nCoanda Surface Radius = {}\nSlot Exit Width = {}\nExit Surface Thickness = {}\nUpper Surface Contour Radius = {}\nInner Surface Contraction Angle = {}\n'.format(Rc, te, tu, ru, ai))

# blockMeshDefaults
scale = 1

# Physical Non control Shape Parameters

pw =  0.0015    # Plenum wall thickness
pl =  0.3       # Plenum external length
p0 =  0.05      # Plenum internal length

# Bounding Box:

bU =  1.0 # Upper bound
bD = -1.0 # Lower bound
bL = -1.0 # Left bound
bR =  1.0 # Right bound
fr =  0.078 # Front bound
bk = -0.078 # Back bound

# Mesh expansion for coanda cylinder

rL = -0.25          # Lower point for cylinder expansion
rU =  0.25          # Upper point for cylinder expansion
sU =  rU + 8 * te   # Upper point for slot expansion
gU =  sU + 4 * tu   # Upper point for top plate expansion

# Dependent Dimensions

hi = 2*Rc + te - pw  # Plenum height
us = hi + pw         # Upper surface of plenum
eu = hi + tu         # Upper surface of slot exit
le = pl + Rc         # Plenum exit x coordinate
ly = hi - te         # Plenum exit y coordinate
pm = (hi)/2          # Middle of plenum for leading edge curve
pn = -(hi/2 + pw)    # LE of curved surface
CC = Rc - pw         # Cylinder center height
ri = (CC + Rc * np.cos(np.radians(ai)))/(1 - np.cos(np.radians(ai)))

# Fillet calculations:

fil = np.zeros([2,2])

fil[0,0] = le - (np.sin(np.radians(ai)) * (Rc+ri))
fil[0,1] = 0
fil[1,0] = le - Rc * np.sin(np.radians(ai))
fil[1,1] = CC + Rc * np.cos(np.radians(ai))

# Additional Coanda Dimensions:


pi =  le - p0/2  # Plenum additional mesh length
pst =  le - p0

hcon = pw - tu
phi = np.arccos(1-(hcon/ru))
fil2 = le - ru*np.sin(phi/2)


# print(fil,p0)
topang = np.degrees(np.arctan(abs(us-eu)/abs(le-fil[0,0]))) # Fillet for upper surface to exit

pback = np.zeros([26,2])

# Define Point Matrix:

pback[ 0,:] = np.array([pst,0])
pback[ 1,:] = fil[0,:]
pback[ 2,:] = np.array([pi,hi])
pback[ 3,:] = np.array([pst,hi])
pback[ 4,:] = fil[1,:]
pback[ 5,:] = np.array([pl,hi])
pback[ 6,:] = np.array([le,ly])
pback[ 7,:] = np.array([le,hi])
pback[ 8,:] = np.array([le,-pw])
pback[ 9,:] = np.array([0,-pw])
pback[10,:] = np.array([0,bD])
pback[11,:] = np.array([le,bD])
pback[12,:] = np.array([bR,bD])
pback[13,:] = np.array([bR,rL])
pback[14,:] = np.array([bR,rU])
pback[15,:] = np.array([bR,sU])
pback[16,:] = np.array([bR,gU])
pback[17,:] = np.array([bR,bU])
pback[18,:] = np.array([le,bU])
pback[19,:] = np.array([fil[0,0],bU])
pback[20,:] = np.array([0,bU])
pback[21,:] = np.array([0,us])
pback[22,:] = np.array([fil2,us])
pback[23,:] = np.array([le,eu])
pback[24,:] = np.array([bL,pm])
pback[25,:] = np.array([pn,pm])

# print(pback[0:10,:])
# plt.plot(pback[0:10,0],pback[0:10,1])
# plt.show()
#pback = np.around(pback,5)

# Finally: Define the blocking and grading parameters!
blocks_x_in  = 15
blocks_y_in  = 15
blocks_x_exf = 5
blocks_x_out = 80
blocks_y_out = 80

grade_x_flat = 1
grade_y_flat = 400

grade_x_curve = 1
grade_y_curve = 400

header = [
    '/*---------------------------------*- C++ -*-----------------------------------*/\n',
    '//    //////////////      ///////////////\n',
    '//    //////////////      //////////////\n',
    '//         ////           ////             /////       ////  //  /// //// ////\n',
    '//         ////     ////  ///////////   ///    ///  ///    ///   ////  ///  ///\n',
    '//  ////   ////           ///////////   ///    ///  ///    ///   ///   ///  ///\n',
    '//   //// ////            ////          ///    ///  ///    ///   ///   ///  ///\n',
    '//     /////              ////            /////       ///// ///  ///   ///  ///\n',
    '/*---------------------------------*- C++ -*-----------------------------------*/\n\n'
]

opener = [
    'FoamFile',
    '\n',
    '{',
    '\n',
    '\t\tversion\t\t\t2.0;',
    '\n',
    '\t\tformat\t\t\tascii;',
    '\n',
    '\t\tclass\t\t\t\tdictionary;',
    '\n',
    '\t\tobject\t\t\tblockMeshDict;',
    '\n',
    '}\n\n',
    '// ***************************************************************************** //\n\n'
]

scale = 'scale\t\t\t\t1;\n\n'

pointsdef = []

def makepoints(points_back,back,front):
    pointsdef = ['vertices\n','(\n']
    for i in range(np.shape(pback)[0]*2):
        if i % 2 == 0:
            pointsdef.append('\t\t({} {} {})\t//\t{}\n'.format(points_back[int(i/2),0],points_back[int(i/2),1],back,i))

        else:
            pointsdef.append('\t\t({} {} {})\t//\t{}\n'.format(points_back[int((i-1)/2),0],points_back[int((i-1)/2),1],front,i))

    pointsdef.append(');\n\n')
    return pointsdef

# Block order is inside coanda plenum 0, 1, 2 then 3 is the block directly below the surface, the rest follow counter-clockwise

blocks = [
    'blocks\n(\n',
    '\t\thex ( 0  2  4  6  1  3  5  7) ({} {} 1) simpleGrading ( 1   1  1) // 0\n'.format(blocks_x_in, blocks_y_in),
    '\t\thex ( 2  8 10  4  3  9 11  5) ({} {} 1) simpleGrading (.25  1  1) // 1\n'.format(blocks_x_in, blocks_y_in),
    '\t\thex ( 8 12 14 10  9 13 15 11) ({} {} 1) simpleGrading ( 1   1  1) // 2\n'.format(blocks_x_in, blocks_y_in),
    '\t\thex (20 22 16 18 21 23 17 19) ({} {} 1) simpleGrading ( 1  {}  1) // 3\n'.format(blocks_x_out+70, blocks_y_out, 1/grade_y_flat),
    '\t\thex (24 26 16 22 25 27 17 23) ({} {} 1) simpleGrading ({}  {}  1) // 4\n'.format(blocks_x_out, blocks_y_out, 1/grade_y_flat, 1/grade_y_curve),
    '\t\thex (26 28 12 16 27 29 13 17) ({} {} 1) simpleGrading ( 1  {}  1) // 5\n'.format(blocks_x_out, blocks_y_out, 1/grade_y_curve),
    '\t\thex (28 30 14 12 29 31 15 13) ({} {} 1) simpleGrading ( 1  {}  1) // 6\n'.format(blocks_y_in, blocks_y_out, 1/grade_y_curve),
    '\t\thex (30 32 46 14 31 33 47 15) ({} {} 1) simpleGrading ( 1  {}  1) // 7\n'.format(blocks_x_exf, blocks_y_out, 1/grade_y_curve),
    '\t\thex (32 34 36 46 33 35 37 47) ({} {} 1) simpleGrading ({}  {}  1) // 8\n'.format(blocks_x_out, blocks_y_out, grade_y_flat, 1/grade_y_curve),
    '\t\thex (36 38 44 46 37 39 45 47) ({} {} 1) simpleGrading ( 1  {}  1) // 9\n'.format(blocks_x_out+30, blocks_y_out, 1/grade_y_flat),
    '\t\thex (38 40 42 44 39 41 43 45) ({} {} 1) simpleGrading ( 1  {}  1) // 10\n'.format(blocks_x_out, blocks_y_out, 1/grade_y_flat),
    '\t\thex (40 48 50 42 41 49 51 43) ({} {} 1) simpleGrading ( 1  {}  1) // 11\n'.format(blocks_x_out, blocks_y_out, 1/grade_y_flat),
    '\t\thex (48 20 18 50 49 21 19 51) ({} {} 1) simpleGrading ( 1  {}  1) // 12\n'.format(blocks_x_out, blocks_y_out, 1/grade_y_flat),
    ');\n\n'
]

edges = [
    'edges\n(\n',
    '\t\tarc   2  8  {} (0 0  1)\n'.format(ai),
    '\t\tarc   3  9  {} (0 0  1)\n'.format(ai),
    '\t\tarc   8 12  {} (0 0 -1)\n'.format(ai),
    '\t\tarc   9 13  {} (0 0 -1)\n'.format(ai),
    '\t\tarc  12 16 180.0 (0 0 -1)\n',
    '\t\tarc  13 17 180.0 (0 0 -1)\n',
    '\t\tarc  44 46  {} (0 0 -1)\n'.format(topang),
    '\t\tarc  45 47  {} (0 0 -1)\n'.format(topang),
    '\t\tarc  40 48  90.0 (0 0  1)\n',
    '\t\tarc  41 49  90.0 (0 0  1)\n',
    '\t\tarc  48 20  90.0 (0 0  1)\n',
    '\t\tarc  49 21  90.0 (0 0  1)\n',
    '\t\tarc  42 50  90.0 (0 0  1)\n',
    '\t\tarc  43 51  90.0 (0 0  1)\n',
    '\t\tarc  50 18  90.0 (0 0  1)\n',
    '\t\tarc  51 19  90.0 (0 0  1)\n',
    ');\n\n'
]

pointsdef = makepoints(pback,bk,fr)


with open('./system/blockMeshDict','w') as f:
    f.writelines(header)
    f.writelines(opener)
    f.writelines(scale)
    f.writelines(pointsdef)
    f.writelines(blocks)
    f.writelines(edges)
    with open('./system/faces','r') as g:
        f.writelines(g.readlines())
