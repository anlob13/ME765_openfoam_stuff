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

Description
    Template for use with codeStream.

\*---------------------------------------------------------------------------*/

#include "dictionaryEntry.H"
#include "fieldTypes.H"
#include "Ostream.H"
#include "Pstream.H"
#include "read.H"
#include "unitConversion.H"

//{{{ begin codeInclude
#line 49 "/home/abrongni/Documents/ME765_Proj/openFoamSims/periodicHill/steadyState/system/blockMeshDict/#codeStream"

        #include "pointField.H"
        #include "mathematicalConstants.H"
    
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
    void codeStream_0c0b30fdd60d7fb8bc3ec0c2c9e9482ba5eac759
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 55 "/home/abrongni/Documents/ME765_Proj/openFoamSims/periodicHill/steadyState/system/blockMeshDict/#codeStream"

        const scalar xMin = 0.0;
        const scalar xMax = 9.0;

        // Scaling between ERCOFTAC source coordinates
        // and the desired 9 x 3.035 x 4.5 geometry
        const scalar scaleGeom = 28.0;

        const label nPoints = 1000;

        const scalar dx =
            (xMax - xMin)/scalar(nPoints - 1);

        os << "(" << nl;
        os << "spline 0 1" << nl;

        pointField profile(nPoints, Zero);

        for (label i = 0; i < nPoints; ++i)
        {
            // New coordinate: 0 -> 9
            scalar x = xMin + i*dx;

            // Convert to original ERCOFTAC coordinate:
            // 0 -> 252
            scalar xOld = scaleGeom*x;

            // Original geometry is defined from crest to base
            // over x = 0 -> 54 and mirrored about the flat
            // region to give a periodic length of 252.
            if (xOld > 198.0)
            {
                xOld = 252.0 - xOld;
            }

            scalar hOld = 28.0;

            if (xOld >= 0.0 && xOld < 9.0)
            {
                hOld =
                    28.0
                  + 6.775070969851E-03*xOld*xOld
                  - 2.124527775800E-03*xOld*xOld*xOld;
            }
            else if (xOld >= 9.0 && xOld < 14.0)
            {
                hOld =
                    25.07355893131
                  + 0.9754803562315*xOld
                  - 1.016116352781E-01*xOld*xOld
                  + 1.889794677828E-03*xOld*xOld*xOld;
            }
            else if (xOld >= 14.0 && xOld < 20.0)
            {
                hOld =
                    2.579601052357E+01
                  + 8.206693007457E-01*xOld
                  - 9.055370274339E-02*xOld*xOld
                  + 1.626510569859E-03*xOld*xOld*xOld;
            }
            else if (xOld >= 20.0 && xOld < 30.0)
            {
                hOld =
                    4.046435022819E+01
                  - 1.379581654948E+00*xOld
                  + 1.945884504128E-02*xOld*xOld
                  - 2.070318932190E-04*xOld*xOld*xOld;
            }
            else if (xOld >= 30.0 && xOld < 40.0)
            {
                hOld =
                    1.792461334664E+01
                  + 8.743920332081E-01*xOld
                  - 5.567361123058E-02*xOld*xOld
                  + 6.277731764683E-04*xOld*xOld*xOld;
            }
            else if (xOld >= 40.0 && xOld < 54.0)
            {
                hOld =
                    max
                    (
                        0.0,
                        5.639011190988E+01
                      - 2.010520359035E+00*xOld
                      + 1.644919857549E-02*xOld*xOld
                      + 2.674976141766E-05*xOld*xOld*xOld
                    );
            }

            // Scale the original 28-unit hill height
            // to a crest height of 1
            const scalar y = hOld/scaleGeom;

            profile[i].x() = x;
            profile[i].y() = y;
            profile[i].z() = 0.0;
        }

        os << profile << nl;

        // Same hill profile on the opposite spanwise face
        os << "spline 4 5" << nl;

        profile.replace(2, 4.5);

        os << profile << nl;

        os << ");" << nl;
    
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

