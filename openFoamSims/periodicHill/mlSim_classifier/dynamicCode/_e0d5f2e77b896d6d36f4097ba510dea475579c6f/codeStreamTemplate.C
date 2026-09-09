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
    void codeStream_e0d5f2e77b896d6d36f4097ba510dea475579c6f
    (
        Ostream& os,
        const dictionary& dict
    )
    {
//{{{ begin code
        #line 162 "/home/abrongni/Documents/ME765_Proj/openFoamSims/periodicHill/mlSim_classifier/system/blockMeshDict/#codeStream"


        os << "(" << nl;


        // ---------------------------------------------------------------
        // Geometry
        // ---------------------------------------------------------------

        const scalar h = 28.0;

        const scalar xMin = 0.0;
        const scalar xMax = 252.0;

        const label nPoints = 1500;

        const scalar dx =
            (xMax - xMin)/scalar(nPoints - 1);

        const scalar yTop = 3.035;


        // ---------------------------------------------------------------
        // Wall-normal locations
        //
        // eta = 0 -> wall
        // eta = 1 -> top
        //
        // These are fractions of the LOCAL NORMAL DISTANCE to the top.
        // ---------------------------------------------------------------

        const scalar eta[9] =
        {
            0.00,
            0.005,
            0.015,
            0.030,
            0.060,
            0.12,
            0.25,
            0.50,
            1.00
        };


        // ---------------------------------------------------------------
        // Periodic hill wall position.
        //
        // Input: x in mm
        // Output: y in mm
        // ---------------------------------------------------------------

        auto wallY =
        [&](const scalar xInput) -> scalar
        {
            scalar x = xInput;

            // Periodic wrapping
            while (x < 0.0)
            {
                x += 252.0;
            }

            while (x > 252.0)
            {
                x -= 252.0;
            }


            scalar xs = x;

            if (xs > 198.0)
            {
                xs = 252.0 - xs;
            }


            scalar y = 0.0;


            if (xs >= 0.0 && xs < 9.0)
            {
                y =
                    28.0
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

            return y;
        };


        // ---------------------------------------------------------------
        // Generate all curved edges.
        // ---------------------------------------------------------------

        for (label layer = 0; layer < 9; ++layer)
        {
            pointField profile(nPoints);


            for (label i = 0; i < nPoints; ++i)
            {
                const scalar xmm =
                    xMin + i*dx;


                // -------------------------------------------------------
                // Wall location
                // -------------------------------------------------------

                const scalar yWallMM =
                    wallY(xmm);


                const scalar xh =
                    xmm/h;

                const scalar yWall =
                    yWallMM/h;


                // -------------------------------------------------------
                // Compute wall slope dy/dx.
                //
                // Numerical derivative is used so the same hill
                // function defines both position and normal.
                // -------------------------------------------------------

                const scalar delta =
                    1.0e-3;


                const scalar yPlus =
                    wallY(xmm + delta);

                const scalar yMinus =
                    wallY(xmm - delta);


                const scalar dydx_mm =
                    (yPlus - yMinus)/(2.0*delta);


                // dy/dx is unchanged by the mm -> h normalization


                // -------------------------------------------------------
                // Unit normal pointing into the flow
                //
                // tangent = (1, dydx)
                //
                // normal  = (-dydx, 1)
                // -------------------------------------------------------

                const scalar magNormal =
                    sqrt(1.0 + sqr(dydx_mm));


                const scalar nx =
                    -dydx_mm/magNormal;

                const scalar ny =
                     1.0/magNormal;


                // -------------------------------------------------------
                // Distance along the normal until it intersects
                // the horizontal upper wall y = 3.035.
                //
                // yWall + d*ny = yTop
                //
                // therefore:
                //
                // d = (yTop-yWall)/ny
                // -------------------------------------------------------

                const scalar dTop =
                    (yTop - yWall)/ny;


                // -------------------------------------------------------
                // Move from wall along the NORMAL.
                //
                // eta = 0 -> wall
                // eta = 1 -> top
                // -------------------------------------------------------

                const scalar d =
                    eta[layer]*dTop;


                const scalar xLayer =
                    xh + d*nx/h;

                const scalar yLayer =
                    yWall + d*ny;


                profile[i].x() = xLayer;
                profile[i].y() = yLayer;
                profile[i].z() = 0.0;
            }


            // -----------------------------------------------------------
            // z = 0 curved edge
            // -----------------------------------------------------------

            const label v0 = 2*layer;
            const label v1 = 2*layer + 1;


            os
                << "spline "
                << v0
                << " "
                << v1
                << nl;

            os << profile << nl;


            // -----------------------------------------------------------
            // z = 4.5 curved edge
            // -----------------------------------------------------------

            const label v0z = 18 + 2*layer;
            const label v1z = 18 + 2*layer + 1;


            os
                << "spline "
                << v0z
                << " "
                << v1z
                << nl;


            profile.replace(2, 4.5);

            os << profile << nl;
        }


        os << ");" << nl;

    
//}}} end code
    }
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //

