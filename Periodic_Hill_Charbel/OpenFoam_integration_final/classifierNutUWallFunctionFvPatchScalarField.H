/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

Class
    Foam::nutUWallFunctionFvPatchScalarField

Description
    This boundary condition provides a turbulent kinematic viscosity condition
    when using wall functions, based on velocity.

Usage
    Example of the boundary condition specification:
    \verbatim
    <patchName>
    {
        type            nutUWallFunction;
        value           uniform 0;
    }
    \endverbatim

See also
    Foam::nutWallFunctionFvPatchScalarField

SourceFiles
    nutUWallFunctionFvPatchScalarField.C

\*---------------------------------------------------------------------------*/

#ifndef classifierNutUWallFunctionFvPatchScalarField_H
#define classifierNutUWallFunctionFvPatchScalarField_H

#include "nutWallFunctionFvPatchScalarField.H"
#include "classifier_wall_model.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    /*---------------------------------------------------------------------------*\
                 Class classifierNutUWallFunctionFvPatchScalarField Declaration
    \*---------------------------------------------------------------------------*/

    class classifierNutUWallFunctionFvPatchScalarField
        : public nutWallFunctionFvPatchScalarField
    {
    protected:
        // Protected data

        // nut ml Wall Model
        // autoPtr<wallModelMLP> mlModel_;

        // nut ml Wall Model Ensemble
        // autoPtr<wallModelEnsemble> ensembleModel_;
        // autoPtr<crossflowWallModel> ensembleModel_;

        // classifier ml
        autoPtr<classifierWallModel> classifierModel_;

        mutable labelList secondLayerCell_;
        mutable bool secondLayerCellValid_;

        //- Global face index of the "top" face of each wall-adjacent cell
        mutable labelList topFace_;
        mutable bool topFaceValid_;

        // Surface coord for basis
        mutable autoPtr<volScalarField> sPtr_;
        mutable bool sValid_ = false;

        // Output field management

        //- Ensure all ML diagnostic fields exist in the mesh registry
        void ensureOutputFields() const;

        //- Set one value in a registered scalar output field
        void setOutputField(
            const word &fieldName,
            const label patchi,
            const label facei,
            const scalar value) const;

        //- Set classifier mean/sigma output fields
        void setClassifierWeightFields(
            const label patchi,
            const label facei,
            const WeightVector &weightsMean,
            const WeightVector &weightsSigma) const;

        // Protected Member Functions

        //- Calculate yPLus
        virtual tmp<scalarField> yPlus(const scalarField &magUp) const;

        //- Calculate the turbulence viscosity
        virtual tmp<scalarField> nut() const;

        // Calculate second layer cell addressing for each face
        void buildSecondLayerAddressing() const;

        //- Calculate top-face addressing for each wall-adjacent cell
        void buildTopFaceAddressing() const;

        //- Interpolate a scalar volField to the top face of facei's cell
        scalar topFaceValue(
            const volScalarField &field,
            const label facei) const;

        //- Interpolate a vector volField to the top face of facei's cell
        vector topFaceValueVector(
            const volVectorField &field,
            const label facei) const;

        void ensureSurfaceCoordinate() const;

    public:
        //- Runtime type information
        TypeName("classifierNutUWallFunction");

        // Constructors

        //- Construct from patch, internal field and dictionary
        classifierNutUWallFunctionFvPatchScalarField(
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const dictionary &);

        //- Construct by mapping given
        //  classifierNutUWallFunctionFvPatchScalarField
        //  onto a new patch
        classifierNutUWallFunctionFvPatchScalarField(
            const classifierNutUWallFunctionFvPatchScalarField &,
            const fvPatch &,
            const DimensionedField<scalar, volMesh> &,
            const fieldMapper &);

        //- Disallow copy without setting internal field reference
        classifierNutUWallFunctionFvPatchScalarField(
            const classifierNutUWallFunctionFvPatchScalarField &) = delete;

        //- Copy constructor setting internal field reference
        classifierNutUWallFunctionFvPatchScalarField(
            const classifierNutUWallFunctionFvPatchScalarField &,
            const DimensionedField<scalar, volMesh> &);

        //- Construct and return a clone setting internal field reference
        virtual tmp<fvPatchScalarField> clone(
            const DimensionedField<scalar, volMesh> &iF) const
        {
            return tmp<fvPatchScalarField>(
                new classifierNutUWallFunctionFvPatchScalarField(*this, iF));
        }

        // Member Functions

        // Evaluation functions

        //- Calculate and return yPlus at the boundary
        virtual tmp<scalarField> yPlus() const;

        // I-O

        //- Write
        virtual void write(Ostream &os) const;
    };

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //
