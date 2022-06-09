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
    '\t\tclass\t\t\t\tvolScalarField;\n',
    '\t\tobject\t\t\tp;\n',
    '}\n',
    '// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //\n\n'
]

contents = [
    'dimensions\t\t\t[1 -1 -2 0 0 0 0];\n\n',
    'internalField\t\tuniform 100000;\n\n',
    'boundaryField\n',
    '{\n',
    '\t\tdown\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tfixedValue;\n',
    '\t\t\t\tvalue\t\t\t\t\t\tuniform 100000;\n',
    '\t\t}\n\n',
    '\t\tright\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tfixedValue;\n',
    '\t\t\t\tvalue\t\t\t\t\t\tuniform 100000;\n',
    '\t\t}\n\n',
    '\t\tup\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tfixedValue;\n',
    '\t\t\t\tvalue\t\t\t\t\t\tuniform 100000;\n',
    '\t\t}\n\n',
    '\t\tmeanFlow\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tfixedValue;\n',
    '\t\t\t\tvalue\t\t\t\t\t\tuniform 100000;\n',
    '\t\t}\n\n',
    '\t\tsurface\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tzeroGradient;\n',
    '\t\t}\n\n',
    '\t\tsurfaceInner\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tzeroGradient;\n',
    '\t\t}\n\n',
    '\t\tcoandaInlet\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\ttotalPressure;\n',
    '\t\t\t\tvalue\t\t\t\t\t\tuniform {};\n'.format(pBC),
    '\t\t\t\tp0\t\t\t\t\t\t\tuniform {};\n'.format(pBC),
    '\t\t}\n\n',
    '\t\tdefaultFaces\n',
    '\t\t{\n',
    '\t\t\t\ttype\t\t\t\t\t\tempty;\n',
    '\t\t}\n',
    '}\n\n',
    '// ************************************************************************* //'
]

with open('./0/p','w') as f:
    f.writelines(header)
    f.writelines(opener)
    f.writelines(contents)
