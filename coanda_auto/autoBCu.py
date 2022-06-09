import sys

pBC = sys.argv[1]

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
    'FoamFile\n',
    '{\n',
    '\t\tversion\t\t\t2.0;\n',
    '\t\tformat\t\t\tascii;\n',
    '\t\tclass\t\t\t\tvolVectorField;\n',
    '\t\tlocation\t\t"0";\n'
    '\t\tobject\t\t\tU;\n',
    '}\n',
    '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n'
]

contents = [
    'dimensions\t\t\t[0 1 -1 0 0 0 0];\n\n',
    'internalField\t\tuniform (30  0  0);\n\n',
    'boundaryField\n',
    '{\n',
    '\t\tdown\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tzeroGradient;\n',
    '\t\t}\n\n',
    '\t\tright\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tzeroGradient;\n',
    '\t\t}\n\n',
    '\t\tup\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tzeroGradient;\n',
    '\t\t}\n\n',
    '\t\tmeanFlow\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tuniformFixedValue;\n',
    '\t\t\t\tuniformValue\t\tconstant (30  0  0);\n',
    '\t\t\t\tvalue\t\t\t\t\t\tuniform (30  0  0);\n',
    '\t\t}\n\n',
    '\t\tsurface\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tnoSlip;\n',
    '\t\t}\n\n',
    '\t\tsurfaceInner\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tnoSlip;\n',
    '\t\t}\n\n',
    '\t\tcoandaInlet\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tflowRateInletVelocity;\n',
    '\t\t\t\tmassFlowRate\t\t{};\n'.format(pBC),
    '\t\t\t\tvalue\t\t\t\t\t\tuniform ({} 0 0);\n'.format(pBC),
    '\t\t}\n\n',
    '\t\tdefaultFaces\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tempty;\n',
    '\t\t}\n',
    '}\n\n',
    '// ************************************************************************* //'
]

with open('./0/U','w') as f:
    f.writelines(header)
    f.writelines(opener)
    f.writelines(contents)
