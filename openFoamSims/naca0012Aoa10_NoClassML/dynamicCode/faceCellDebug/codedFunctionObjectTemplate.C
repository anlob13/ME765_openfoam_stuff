/*---------------------------------------------------------------------------*\
  =========                 |
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

#include "codedFunctionObjectTemplate.H"
#include "volFields.H"
#include "read.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(faceCellDebugFunctionObject, 0);

addRemovableToRunTimeSelectionTable
(
    functionObject,
    faceCellDebugFunctionObject,
    dictionary
);


// * * * * * * * * * * * * * * * Global Functions  * * * * * * * * * * * * * //

extern "C"
{
    // dynamicCode:
    // SHA1 = 60cb972956e3b72c0d2b92b3b1a187aa11914d79
    //
    // unique function name that can be checked if the correct library version
    // has been loaded
    void faceCellDebug_60cb972956e3b72c0d2b92b3b1a187aa11914d79(bool load)
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


// * * * * * * * * * * * * * * * Local Functions * * * * * * * * * * * * * * //

//{{{ begin localCode

//}}} end localCode


// * * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * //

const fvMesh& faceCellDebugFunctionObject::mesh() const
{
    return refCast<const fvMesh>(obr_);
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

faceCellDebugFunctionObject::faceCellDebugFunctionObject
(
    const word& name,
    const Time& runTime,
    const dictionary& dict
)
:
    functionObjects::regionFunctionObject(name, runTime, dict)
{
    read(dict);
}


// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

faceCellDebugFunctionObject::~faceCellDebugFunctionObject()
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool faceCellDebugFunctionObject::read(const dictionary& dict)
{
    if (false)
    {
        Info<<"read faceCellDebug sha1: 60cb972956e3b72c0d2b92b3b1a187aa11914d79\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


Foam::wordList faceCellDebugFunctionObject::fields() const
{
    if (false)
    {
        Info<<"fields faceCellDebug sha1: 60cb972956e3b72c0d2b92b3b1a187aa11914d79\n";
    }

    wordList fields;
//{{{ begin code
    
//}}} end code

    return fields;
}


bool faceCellDebugFunctionObject::execute()
{
    if (false)
    {
        Info<<"execute faceCellDebug sha1: 60cb972956e3b72c0d2b92b3b1a187aa11914d79\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


bool faceCellDebugFunctionObject::write()
{
    if (false)
    {
        Info<<"write faceCellDebug sha1: 60cb972956e3b72c0d2b92b3b1a187aa11914d79\n";
    }

//{{{ begin code
    #line 30 "/home/abrongni/Documents/ME765_Proj/openFoamSims/aerofoilNACA0012Test/system/functions/faceCellDebug"

        const fvMesh& mesh = this->mesh();

        // Change this to your actual aerofoil patch name
        const word patchName = "aerofoil";
        const label patchi = mesh.boundaryMesh().findIndex(patchName);

        if (patchi == -1)
        {
            FatalErrorInFunction
                << "Patch " << patchName << " not found"
                << exit(FatalError);
        }

        const fvPatch& p = mesh.boundary()[patchi];
        const labelUList& faceCells = p.faceCells();

        // Register (or fetch) a debug field so it shows up in ParaView
        if (!mesh.foundObject<volScalarField>("faceCellDebug"))
        {
            volScalarField* debugFieldPtr = new volScalarField
            (
                IOobject
                (
                    "faceCellDebug",
                    mesh.time().timeName(mesh.time().value()),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(dimless, -1)
            );
            debugFieldPtr->store();
        }

        volScalarField& debugField =
            mesh.lookupObjectRef<volScalarField>("faceCellDebug");

        // Reset then fill with facei at each near-wall cell
        debugField.primitiveFieldRef() = -1;

        forAll(faceCells, facei)
        {
            debugField[faceCells[facei]] = facei;
        }

        // Also print a text table to the log for exact cross-checking
        Info<< "=== faceCellDebug mapping for patch " << patchName
            << " (proc " << Pstream::myProcNo() << ") ===" << endl;

        forAll(faceCells, facei)
        {
            const label celli = faceCells[facei];
            Info<< "facei " << facei
                << " -> celli " << celli
                << " @ " << mesh.C()[celli]
                << endl;
        }

        debugField.write();
    
//}}} end code

    return true;
}


bool faceCellDebugFunctionObject::end()
{
    if (false)
    {
        Info<<"end faceCellDebug sha1: 60cb972956e3b72c0d2b92b3b1a187aa11914d79\n";
    }

//{{{ begin code
    
//}}} end code

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

