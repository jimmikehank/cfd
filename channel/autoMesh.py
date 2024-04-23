### Generates channel mesh in blockMeshDict

import numpy as np
import argparse


# Dimensions of Channel

x = 2.5
y = 1
z = 1

block_x = 100
block_y = 40
block_z = 40

grade_y = 8
grade_z = 8

parser = argparse.ArgumentParser()
parser.add_argument('--clean', default = False, type = bool, help = 'Run clean function. \nDefault: False')
args = parser.parse_args()
clean_bool = args.clean

retain = ['autoMesh.py','system','0','constant']


def cleanup(retain):
    import shutil
    import os
    dirs = os.listdir()
    for item in dirs:
        if item not in retain and os.path.isdir(item):
            shutil.rmtree(item)
        elif item not in retain and os.path.isfile(item):
            os.remove(item)
        else:
            continue


if clean_bool:
    cleanup(retain)
    
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


pointsdef = [
    'vertices\n',
    '(\n'
    '\t\t( {}\t{}\t{} ) // 0\n'.format(0,-y/2,-z/2),
    '\t\t( {}\t{}\t{} ) // 1\n'.format(x, -y/2,-z/2),
    '\t\t( {}\t{}\t{} ) // 2\n'.format(x,  y/2,-z/2),
    '\t\t( {}\t{}\t{} ) // 3\n'.format(0,  y/2,-z/2),
    '\t\t( {}\t{}\t{} ) // 4\n'.format(0, -y/2, z/2),
    '\t\t( {}\t{}\t{} ) // 5\n'.format(x, -y/2, z/2),
    '\t\t( {}\t{}\t{} ) // 6\n'.format(x,  y/2, z/2),
    '\t\t( {}\t{}\t{} ) // 7\n'.format(0,  y/2, z/2),
    ');\n\n'
]

blocksdef = [
    'blocks\n',
    '(\n',
    '\t\thex (0 1 2 3 4 5 6 7) ({} {} {}) simpleGrading (1 ((.5 .5 {}) (.5 .5 {})) ((.5 .5 {}) (.5 .5 {})))\n'.format(block_x, block_y, block_z, grade_y, 1/grade_y, grade_z, 1/grade_z),
    ');\n\n'
]

boundarydef = [
    'boundary\n',
    '(\n',
    '\t\tinlet\n',
    '\t\t{\n',
    '\t\t\t\ttype patch;\n',
    '\t\t\t\tfaces ((0 4 7 3));\n',
    '\t\t}\n',
    '\t\toutlet\n',
    '\t\t{\n',
    '\t\t\t\ttype patch;\n',
    '\t\t\t\tfaces ((1 5 6 2));\n',
    '\t\t}\n',
    '\t\tleft\n',
    '\t\t{\n',
    '\t\t\t\ttype wall;\n',
    '\t\t\t\tfaces ((0 1 2 3));\n',
    '\t\t}\n',
    '\t\tdown\n',
    '\t\t{\n',
    '\t\t\t\ttype wall;\n',
    '\t\t\t\tfaces ((0 4 5 1));\n',
    '\t\t}\n',
    '\t\tright\n',
    '\t\t{\n',
    '\t\t\t\ttype wall;\n',
    '\t\t\t\tfaces ((4 5 6 7));\n',
    '\t\t}\n',
    '\t\tup\n',
    '\t\t{\n',
    '\t\t\t\ttype wall;\n',
    '\t\t\t\tfaces ((7 6 2 3));\n',
    '\t\t}\n',
    ');\n\n'

]


with open('./system/blockMeshDict','w') as f:
    f.writelines(header)
    f.writelines(opener)
    f.writelines(scale)
    f.writelines(pointsdef)
    f.writelines(blocksdef)
    f.writelines(boundarydef)
