/*---------------------------------------------------------------------------*  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) YEAR OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "fixedValueFvPatchFieldTemplate.H"
#include "addToRunTimeSelectionTable.H"
#include "fvPatchFieldMapper.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "unitConversion.H"
//{{{ begin codeInclude

//}}} end codeInclude


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 3e04f7a92a5b89d44c98f1ecde4cc203dfbfcc1b
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void swirl_3e04f7a92a5b89d44c98f1ecde4cc203dfbfcc1b(bool load)
    {
        if (load)
        {
            // code that can be explicitly executed after loading
        }
        else
        {
            // code that can be explicitly executed before unloading
        }
    }
}

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

makeRemovablePatchTypeField
(
    fvPatchVectorField,
    swirlFixedValueFvPatchVectorField
);


const char* const swirlFixedValueFvPatchVectorField::SHA1sum =
    "3e04f7a92a5b89d44c98f1ecde4cc203dfbfcc1b";


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

swirlFixedValueFvPatchVectorField::
swirlFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(p, iF)
{
    if (false)
    {
        Info<<"construct swirl sha1: 3e04f7a92a5b89d44c98f1ecde4cc203dfbfcc1b"
            " from patch/DimensionedField\n";
    }
}


swirlFixedValueFvPatchVectorField::
swirlFixedValueFvPatchVectorField
(
    const swirlFixedValueFvPatchVectorField& ptf,
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const fvPatchFieldMapper& mapper
)
:
    fixedValueFvPatchField<vector>(ptf, p, iF, mapper)
{
    if (false)
    {
        Info<<"construct swirl sha1: 3e04f7a92a5b89d44c98f1ecde4cc203dfbfcc1b"
            " from patch/DimensionedField/mapper\n";
    }
}


swirlFixedValueFvPatchVectorField::
swirlFixedValueFvPatchVectorField
(
    const fvPatch& p,
    const DimensionedField<vector, volMesh>& iF,
    const dictionary& dict
)
:
    fixedValueFvPatchField<vector>(p, iF, dict)
{
    if (false)
    {
        Info<<"construct swirl sha1: 3e04f7a92a5b89d44c98f1ecde4cc203dfbfcc1b"
            " from patch/dictionary\n";
    }
}


swirlFixedValueFvPatchVectorField::
swirlFixedValueFvPatchVectorField
(
    const swirlFixedValueFvPatchVectorField& ptf
)
:
    fixedValueFvPatchField<vector>(ptf)
{
    if (false)
    {
        Info<<"construct swirl sha1: 3e04f7a92a5b89d44c98f1ecde4cc203dfbfcc1b"
            " as copy\n";
    }
}


swirlFixedValueFvPatchVectorField::
swirlFixedValueFvPatchVectorField
(
    const swirlFixedValueFvPatchVectorField& ptf,
    const DimensionedField<vector, volMesh>& iF
)
:
    fixedValueFvPatchField<vector>(ptf, iF)
{
    if (false)
    {
        Info<<"construct swirl sha1: 3e04f7a92a5b89d44c98f1ecde4cc203dfbfcc1b "
            "as copy/DimensionedField\n";
    }
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

swirlFixedValueFvPatchVectorField::
~swirlFixedValueFvPatchVectorField()
{
    if (false)
    {
        Info<<"destroy swirl sha1: 3e04f7a92a5b89d44c98f1ecde4cc203dfbfcc1b\n";
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

void swirlFixedValueFvPatchVectorField::updateCoeffs()
{
    if (this->updated())
    {
        return;
    }

    if (false)
    {
        Info<<"updateCoeffs swirl sha1: 3e04f7a92a5b89d44c98f1ecde4cc203dfbfcc1b\n";
    }

//{{{ begin code
    #line 30 "/home/james/Documents/research/cfd/pipe/0/U/boundaryField/inlet"
const vector axis(1, 0, 0);

            vectorField v(2.0*this->patch().Cf() ^ axis);
            v.replace(vector::X, 1.0);
            operator==(v);
//}}} end code

    this->fixedValueFvPatchField<vector>::updateCoeffs();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

