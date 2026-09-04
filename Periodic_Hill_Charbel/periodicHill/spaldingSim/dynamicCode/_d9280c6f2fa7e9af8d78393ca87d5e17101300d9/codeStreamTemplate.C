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
#line 45 "/home/abrongni/Documents/ME765_Proj/openFoamSims/periodicHill/steadyState/system/blockMeshDict/#codeStream"

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
    void codeStream_d9280c6f2fa7e9af8d78393ca87d5e17101300d9
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 51 "/home/abrongni/Documents/ME765_Proj/openFoamSims/periodicHill/steadyState/system/blockMeshDict/#codeStream"

        // Hill profile polynomials are fitted in mm (original ERCOFTAC
        // coefficients); h is the normalisation length (crest height).
        const scalar h = 28.0;
        const scalar xMin = 0;
        const scalar xMax = 252;
        const label nPoints = 1000;
        const scalar dx = (xMax - xMin)/scalar(nPoints - 1);

        os  << "(" << nl << "spline 0 1" << nl;
        pointField profile(nPoints);

        for (label i = 0; i < nPoints; ++i)
        {
            const scalar x = xMin + i*dx;
            scalar xs = x;
            if (xs > 198) xs = 252 - xs;

            scalar y = 0;

            if (xs >= 0 && xs < 9)
            {
                y = min
                    (
                        28.0,
                        2.800000000000E+01
                      + 0.000000000000E+00*xs
                      + 6.775070969851E-03*xs*xs
                      - 2.124527775800E-03*xs*xs*xs
                    );
            }
            else if (xs >= 9 && xs < 14)
            {
                y = 2.507355893131E+01
                  + 9.754803562315E-01*xs
                  - 1.016116352781E-01*xs*xs
                  + 1.889794677828E-03*xs*xs*xs;
            }
            else if (xs >= 14 && xs < 20)
            {
                y = 2.579601052357E+01
                  + 8.206693007457E-01*xs
                  - 9.055370274339E-02*xs*xs
                  + 1.626510569859E-03*xs*xs*xs;
            }
            else if (xs >= 20 && xs < 30)
            {
                y = 4.046435022819E+01
                  - 1.379581654948E+00*xs
                  + 1.945884504128E-02*xs*xs
                  - 2.070318932190E-04*xs*xs*xs;
            }
            else if (xs >= 30 && xs < 40)
            {
                y = 1.792461334664E+01
                  + 8.743920332081E-01*xs
                  - 5.567361123058E-02*xs*xs
                  + 6.277731764683E-04*xs*xs*xs;
            }
            else if (xs >= 40 && xs < 54)
            {
                y = max
                    (
                        0.0,
                        5.639011190988E+01
                      - 2.010520359035E+00*xs
                      + 1.644919857549E-02*xs*xs
                      + 2.674976141766E-05*xs*xs*xs
                    );
            }
            else
            {
                y = 0;
            }

            // Normalise to h-units so the block spans exactly 0-9 in x
            // and the hill crest sits at y = 1.
            profile[i].x() = x/h;
            profile[i].y() = y/h;
            profile[i].z() = 0;
        }
        os << profile << nl;

        os << "spline 4 5" << nl;
        profile.replace(2, 4.5);
        os << profile << nl;

        os  << ");" << nl;
    
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

