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
#line 156 "/home/abrongni/Documents/ME765_Proj/openFoamSims/periodicHill/mlSim_classifier/system/blockMeshDict/#codeStream"

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
    void codeStream_b5bfa76a23b7d8bf5b28d6d961d9ac01b1b5d5ef
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 162 "/home/abrongni/Documents/ME765_Proj/openFoamSims/periodicHill/mlSim_classifier/system/blockMeshDict/#codeStream"


        const scalar h = 28.0;

        const scalar xMin = 0.0;
        const scalar xMax = 252.0;

        const label nPoints = 1000;

        const scalar dx =
            (xMax - xMin)/scalar(nPoints - 1);


        // ---------------------------------------------------------------
        // Wall-normal layer locations.
        //
        // eta = 0 -> hill wall
        // eta = 1 -> channel top
        //
        // Strong clustering near the wall.
        // ---------------------------------------------------------------

        const scalar eta[9] =
        {
            0.00,
            0.02,
            0.05,
            0.10,
            0.20,
            0.35,
            0.55,
            0.775,
            1.00
        };


        // ---------------------------------------------------------------
        // Generate each streamwise layer.
        // ---------------------------------------------------------------

        for (label layer = 0; layer < 9; ++layer)
        {
            pointField profile(nPoints);

            for (label i = 0; i < nPoints; ++i)
            {
                const scalar x = xMin + i*dx;

                scalar xs = x;

                if (xs > 198.0)
                {
                    xs = 252.0 - xs;
                }


                scalar y = 0.0;


                // -------------------------------------------------------
                // ERCOFTAC periodic hill
                // -------------------------------------------------------

                if (xs >= 0.0 && xs < 9.0)
                {
                    y =
                        2.800000000000E+01
                      + 0.000000000000E+00*xs
                      + 6.775070969851E-03*xs*xs
                      - 2.124527775800E-03*xs*xs*xs;

                    y = min(28.0, y);
                }
                else if (xs >= 9.0 && xs < 14.0)
                {
                    y =
                        2.507355893131E+01
                      + 9.754803562315E-01*xs
                      - 1.016116352781E-01*xs*xs
                      + 1.889794677828E-03*xs*xs*xs;
                }
                else if (xs >= 14.0 && xs < 20.0)
                {
                    y =
                        2.579601052357E+01
                      + 8.206693007457E-01*xs
                      - 9.055370274339E-02*xs*xs
                      + 1.626510569859E-03*xs*xs*xs;
                }
                else if (xs >= 20.0 && xs < 30.0)
                {
                    y =
                        4.046435022819E+01
                      - 1.379581654948E+00*xs
                      + 1.945884504128E-02*xs*xs
                      - 2.070318932190E-04*xs*xs*xs;
                }
                else if (xs >= 30.0 && xs < 40.0)
                {
                    y =
                        1.792461334664E+01
                      + 8.743920332081E-01*xs
                      - 5.567361123058E-02*xs*xs
                      + 6.277731764683E-04*xs*xs*xs;
                }
                else if (xs >= 40.0 && xs < 54.0)
                {
                    y =
                        5.639011190988E+01
                      - 2.010520359035E+00*xs
                      + 1.644919857549E-02*xs*xs
                      + 2.674976141766E-05*xs*xs*xs;

                    y = max(0.0, y);
                }
                else
                {
                    y = 0.0;
                }


                // -------------------------------------------------------
                // Convert mm -> h units.
                //
                // Wall:
                //     y_wall = y/h
                //
                // Channel top:
                //     y_top = 3.035
                //
                // Interpolate between wall and top.
                // -------------------------------------------------------

                const scalar yWall = y/h;

                const scalar yTop = 3.035;

                const scalar yLayer =
                    yWall + eta[layer]*(yTop - yWall);


                profile[i].x() = x/h;
                profile[i].y() = yLayer;
                profile[i].z() = 0.0;
            }


            // z = 0 edge

            const label v0 = 2*layer;

            const label v1 = 2*layer + 1;

            os << "spline "
               << v0 << " "
               << v1
               << nl;

            os << profile << nl;


            // z = 4.5 edge

            const label v0z = 18 + 2*layer;

            const label v1z = 18 + 2*layer + 1;

            os << "spline "
               << v0z << " "
               << v1z
               << nl;

            profile.replace(2, 4.5);

            os << profile << nl;
        }

    
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

