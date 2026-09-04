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

\*---------------------------------------------------------------------------*/

#include "unifiedNutUWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"
#include <array>
#include "wallPolyPatch.H"
#include "fvmLaplacian.H"
#include "fvMatrices.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

    // * * * * * * * * * * * * * * Output field names * * * * * * * * * * * * * //

    namespace
    {
        const std::array<word, WM_N_REGIMES> regimeWeightMeanNames =
            {
                word("classifier_w_laminar_mean"),
                word("classifier_w_freestream_mean"),
                word("classifier_w_channel_mean"),
                word("classifier_w_couette_mean"),
                word("classifier_w_crossflow_mean")};

        const std::array<word, WM_N_REGIMES> regimeWeightSigmaNames =
            {
                word("classifier_w_laminar_sigma"),
                word("classifier_w_freestream_sigma"),
                word("classifier_w_channel_sigma"),
                word("classifier_w_couette_sigma"),
                word("classifier_w_crossflow_sigma")};
    }

    // * * * * * * * * * * * Protected Member Functions * * * * * * * * * * * //

    void unifiedNutUWallFunctionFvPatchScalarField::ensureOutputFields() const
    {
        const fvMesh &mesh = patch().boundaryMesh().mesh();

        auto createField =
            [&](const word &fieldName)
        {
            if (!mesh.foundObject<volScalarField>(fieldName))
            {
                volScalarField *fieldPtr =
                    new volScalarField(
                        IOobject(
                            fieldName,
                            mesh.time().name(),
                            mesh,
                            IOobject::NO_READ,
                            IOobject::AUTO_WRITE),
                        mesh,
                        dimensionedScalar(
                            fieldName,
                            dimless,
                            0.0));

                const_cast<fvMesh &>(mesh).store(fieldPtr);
            }
        };

        // ML wall-model outputs
        createField("tw_n_mean");
        createField("tw_n_sigma");

        // Classifier regime means
        for (int r = 0; r < WM_N_REGIMES; ++r)
        {
            createField(regimeWeightMeanNames[r]);
        }

        // Classifier regime sigmas
        for (int r = 0; r < WM_N_REGIMES; ++r)
        {
            createField(regimeWeightSigmaNames[r]);
        }
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    void unifiedNutUWallFunctionFvPatchScalarField::setOutputField(
        const word &fieldName,
        const label patchi,
        const label facei,
        const scalar value) const
    {
        const fvMesh &mesh = patch().boundaryMesh().mesh();

        if (!mesh.foundObject<volScalarField>(fieldName))
        {
            FatalErrorInFunction
                << "Output field " << fieldName
                << " has not been created in the mesh registry."
                << exit(FatalError);
        }

        volScalarField &field =
            const_cast<volScalarField &>(
                mesh.lookupObject<volScalarField>(fieldName));

        field.boundaryFieldRef()[patchi][facei] = value;
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    void unifiedNutUWallFunctionFvPatchScalarField::setClassifierWeightFields(
        const label patchi,
        const label facei,
        const WeightVector &weightsMean,
        const WeightVector &weightsSigma) const
    {
        ensureOutputFields();

        for (int r = 0; r < WM_N_REGIMES; ++r)
        {
            setOutputField(
                regimeWeightMeanNames[r],
                patchi,
                facei,
                weightsMean[r]);

            setOutputField(
                regimeWeightSigmaNames[r],
                patchi,
                facei,
                weightsSigma[r]);
        }
    }

    // * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

    tmp<scalarField> unifiedNutUWallFunctionFvPatchScalarField::nut() const
    {
        const label patchi = patch().index();

        // Ensure all diagnostic fields exist before evaluating wall faces
        ensureOutputFields();

        const momentumTransportModel &turbModel =
            db().lookupType<momentumTransportModel>(internalField().group());

        const tmp<scalarField> tnuw = turbModel.nu(patchi);
        const scalarField &nuw = tnuw();

        const volVectorField &Ufield = turbModel.U(); // for top layer

        const fvPatchVectorField &Uw = Ufield.boundaryField()[patchi];
        const scalarField magUp(mag(Uw.patchInternalField() - Uw));

        const tmp<volScalarField> tnuField = turbModel.nu();
        const volScalarField &nuField = tnuField();

        const volScalarField &yField = turbModel.y();
        const scalarField y(yField.boundaryField()[patchi].patchInternalField());

        // NOTE: an unconditional `db().lookupObject<volScalarField>("rho")`
        // lookup used to sit here. Removed - confirmed dead code (both its uses,
        // `rhop[facei]`, were already commented out below, so removing it
        // changes no computed output) that additionally made this class
        // uncompilable/crashing at runtime against any INCOMPRESSIBLE solver
        // (no `rho` field is ever registered there - confirmed via a real
        // end-to-end run against OpenFOAM-13's incompressibleFluid solver,
        // which aborted here with "request for volScalarField rho ... failed"
        // before this fix). If a future change genuinely needs density, look
        // it up defensively (e.g. `db().foundObject<volScalarField>("rho") ?
        // ... : dimensionedScalar("rho", dimDensity, 1.0)`), not unconditionally.
        tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
        scalarField &nutw = tnutw.ref();

        tmp<scalarField> tnutwTest(new scalarField(patch().size(), 0.0));
        scalarField &nutwTest = tnutwTest.ref();

        // compute acceleration

        const fvMesh &mesh = patch().boundaryMesh().mesh();

        // const volVectorField& U = mesh.lookupObject<volVectorField>("U");
        const fvPatchVectorField &Uw_old = Ufield.oldTime().boundaryField()[patchi];
        const scalarField magUp_old(mag(Uw_old.patchInternalField() - Uw_old));

        const scalar deltaT = mesh.time().deltaTValue();

        // Current and previous near-wall cell velocity
        const vectorField Uc_new(Uw.patchInternalField());
        const vectorField Uc_old(Uw_old.patchInternalField());

        // Simple backward-difference local acceleration
        vectorField dudt_local((Uc_new - Uc_old) / deltaT);

        // --- Build the local wall-fitted basis ---

        // Wall-normal unit vector (outward from wall into domain, per face)
        const vectorField nf(patch().nf());

        // Basis computed from geometry

        // Free-stream / reference direction
        // const vector chordDir(1, 0, 0);
        // Project chord direction onto the local wall-tangent plane at each face
        // vectorField t1(chordDir - (chordDir & nf) * nf);
        // t1 /= (mag(t1) + SMALL);
        // const vectorField t2(nf ^ t1);

        ensureSurfaceCoordinate();
        const volScalarField &s = sPtr_();

        // Streamwise direction: tangential component of grad(s)
        vectorField t1(fvc::grad(s)().boundaryField()[patchi].patchInternalField());
        t1 -= (t1 & nf) * nf; // project onto tangent plane

        const scalarField magT1(mag(t1));
        forAll(t1, facei)
        {
            if (magT1[facei] > SMALL)
            {
                t1[facei] /= magT1[facei];
            }
            else
            {
                // Degenerate face (grad(s) ~ 0, e.g. numerical flat spot) -
                // fall back to an arbitrary consistent in-plane direction
                // rather than dividing by SMALL and getting garbage.
                vector fallback(1, 0, 0);
                fallback -= (fallback & nf[facei]) * nf[facei];
                if (mag(fallback) < SMALL)
                {
                    fallback = vector(0, 1, 0);
                    fallback -= (fallback & nf[facei]) * nf[facei];
                }
                t1[facei] = fallback / mag(fallback);
            }
        }

        // Spanwise direction: complete the right-handed orthonormal triad
        vectorField t2(nf ^ t1);
        t2 /= (mag(t2) + SMALL);

        // Re-orthogonalize t1 against t2 (guards against small numerical
        // non-orthogonality after the fallback branch above)
        t1 = t2 ^ nf;
        t1 /= (mag(t1) + SMALL);

        // basis sanity check
        if (mesh.time().writeTime()) // or your own gate, so this doesn't run every timestep
        {
            if (!mesh.foundObject<volVectorField>("t1_basis"))
            {
                volVectorField *t1Field = new volVectorField(
                    IOobject(
                        "t1_basis",
                        mesh.time().timeName(mesh.time().value()),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE),
                    mesh,
                    dimensionedVector(dimless, Zero));
                t1Field->store();
            }

            if (!mesh.foundObject<volVectorField>("t2_basis"))
            {
                volVectorField *t2Field = new volVectorField(
                    IOobject(
                        "t2_basis",
                        mesh.time().timeName(mesh.time().value()),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::AUTO_WRITE),
                    mesh,
                    dimensionedVector(dimless, Zero));
                t2Field->store();
            }

            volVectorField &t1Field =
                const_cast<volVectorField &>(
                    mesh.lookupObject<volVectorField>("t1_basis"));
            volVectorField &t2Field =
                const_cast<volVectorField &>(
                    mesh.lookupObject<volVectorField>("t2_basis"));

            t1Field.boundaryFieldRef()[patchi] == t1;
            t2Field.boundaryFieldRef()[patchi] == t2;
        }
        // Project acceleration onto the basis

        scalarField dudt_n(dudt_local & (-nf));   // wall-normal
        scalarField dudt_stream(dudt_local & t1); // streamwise
        scalarField dudt_span(dudt_local & t2);   // spanwise

        // scalarField dudt_n(dudt_local.component(vector::Y));       // wall-normal
        // scalarField dudt_stream(dudt_local.component(vector::X));  // streamwise
        // scalarField dudt_span(dudt_local.component(vector::Z));    // spanwise

        // end compute acceleration

        // compute kinetic energy (RANS)

        scalarField meanKE(0.5 * magSqr(Uc_new)); // 0.5 * (U·U), per cell

        // const volScalarField& k = mesh.lookupObject<volScalarField>("k"); //RANS only
        // scalarField turbKE(k.boundaryField()[patchi].patchInternalField());

        // For LES
        scalarField turbKE(patch().size(), 0.0);

        // Subgrid scale TKE
        tmp<volScalarField> tkSgs = turbModel.k(); // valid for all LES models, RANS models too
        const volScalarField &kSgs = tkSgs();
        scalarField kSgsPatch(kSgs.boundaryField()[patchi].patchInternalField());

        // Resolved TKE
        scalarField kResolvedPatch(patch().size(), 0.0);

        if (turbModel.type() != "RAS" && mesh.foundObject<volSymmTensorField>("UPrime2Mean")) // deactivate for RANS
        {
            const volSymmTensorField &UPrime2Mean =
                mesh.lookupObject<volSymmTensorField>("UPrime2Mean");

            const symmTensorField UPrime2MeanPatch(
                UPrime2Mean.boundaryField()[patchi].patchInternalField());

            kResolvedPatch = 0.5 * tr(UPrime2MeanPatch);
        }

        turbKE = kSgsPatch + kResolvedPatch; // SGS + Resolved

        // end compute kinetic energy

        const scalarField yPlus(this->yPlus(magUp));

        const scalar Aplus = 26.0; // For LES wall function, not used in RANS

        static Foam::OFstream wmlog("nutwmlog.log");

        forAll(yPlus, facei)
        {
            if (yPlus[facei] > yPlusLam_)
            {
                nutwTest[facei] = nuw[facei] * (yPlus[facei] * kappa_ / log(E_ * yPlus[facei]) - 1);
            }
        }

        // wallModelMLP nutMlModel(this->db().time().constant()/"wallModel");

        // if (!secondLayerCellValid_)
        //{
        //     buildSecondLayerAddressing();
        // }

        if (!topFaceValid_)
        {
            buildTopFaceAddressing();
        }

        std::cout << "nutU wall function loop" << std::endl; // debug

        forAll(nutw, facei)
        {
            // Top face values

            const scalar yT = topFaceValue(yField, facei);
            const vector UT = topFaceValueVector(Ufield, facei);
            const scalar magUpT = mag(UT - Uw[facei]);
            const scalar nuwT = topFaceValue(nuField, facei);

            // Streamwise magnitudes

            const scalar magU_stream = mag(Uc_new[facei] & t1[facei]);
            const scalar magUT_stream = mag(UT & t1[facei]);

            // spanwise magnitudes
            const scalar magU_span = mag(Uc_new[facei] & t2[facei]);
            const scalar magUT_span = mag(UT & t2[facei]);

            // Streamwise Reynolds

            const scalar centerRe = y[facei] * sqrt(magU_stream * magU_stream + magU_span * magU_span) / nuw[facei];
            const scalar topRe = yT * sqrt(magU_stream * magU_stream + magU_span * magU_span) / nuwT;

            // Angles
            const scalar angle_center_top = std::remainder(atan2(UT & t2[facei], UT & t1[facei]) - atan2(Uc_new[facei] & t2[facei], Uc_new[facei] & t1[facei]), 2 * M_PI);
            // const scalar angle_center_top_prev = std::remainder(atan2(topFaceValueVector(Ufield.oldTime(), facei) & t2[facei], topFaceValueVector(Ufield.oldTime(), facei) & t1[facei]) - atan2(Uc_old[facei] & t2[facei], Uc_old[facei] & t1[facei]), 2 * M_PI);
            // const scalar dtheta_dt = ((angle_center_top - angle_center_top_prev)/deltaT)*pow(y[facei],2)/nuw[facei];
            const scalar dtheta_dt = (std::remainder((atan2(Uc_new[facei] & t2[facei], Uc_new[facei] & t1[facei]) - atan2(Uc_old[facei] & t2[facei], Uc_old[facei] & t1[facei])), 2 * M_PI) / deltaT) * sqr(y[facei]) / nuw[facei];

            // Inputs

            const scalar dudt_stream_input = dudt_stream[facei] * pow(y[facei], 3) / pow(nuw[facei], 2);
            // const scalar dudt_n_input = dudt_n[facei]*pow(y[facei], 3)/pow(nuw[facei], 2);
            // const scalar dudt_span_input = dudt_span[facei]*pow(y[facei], 3)/pow(nuw[facei], 2);

            std::vector<scalar> inputs = {centerRe, topRe, dudt_stream_input, angle_center_top, dtheta_dt, turbKE[facei] / max(SMALL, meanKE[facei])};

            unifiedWallModel::Prediction values = classifierModel_->predict(inputs);

            const scalar tw_n_prediction = max(0, values.twMean);
            // values.twSigma already includes the (C,Z)-calibration-driven
            // epistemic term (see unified_wall_model.C's predict()) on top of
            // the NN aleatoric+epistemic and classifier terms it already had -
            // no separate combination needed here, unifiedWallModel::predict()
            // already did it. Written into a registered volScalarField below
            // (previously computed here and only ever printed to a debug log,
            // never exposed as an actual field - see writeTwField()).
            const scalar tw_n_sigma = values.twSigma;
            const scalar nutw_sigma =
                (values.twMean > 0)
                    ? nuw[facei] * (magU_stream / max(SMALL, magUp[facei])) * tw_n_sigma
                    : 0;

            // std::cout << "utau flag" << nuw[facei] *tw_n_prediction * magU_stream / y[facei] << std::endl; //debug
            // Info << "Reached utau" << endl;
            const scalar uTau = sqrt(nuw[facei] * tw_n_prediction * magU_stream / max(SMALL, y[facei]));

            // Store ML outputs
            setOutputField(
                "tw_n_mean",
                patchi,
                facei,
                uTau*uTau);

            setOutputField(
                "tw_n_sigma",
                patchi,
                facei,
                tw_n_sigma);

            setClassifierWeightFields(
                patchi,
                facei,
                values.weightsMean,
                values.weightsSigma);


            // const scalar uTau = ensembleModel_->predict(inputs)[0]*magUp[facei]*nuw[facei]/(y[facei]*rhop[facei]);
            const scalar yPlusCalc = uTau * y[facei] / max(SMALL, nuw[facei]);

            if (turbModel.type() == "RAS") // only for RANS
            {
                nutw[facei] = max(0, (uTau * uTau * y[facei]) / max(SMALL, magUp[facei]) - nuw[facei]); // RANS
            }
            else
            {
                nutw[facei] = kappa_ * uTau * y[facei] * (1 - std::exp(-yPlusCalc / Aplus)) * (1 - std::exp(-yPlusCalc / Aplus)); // LES
            }

            //const scalar expected = nutwTest[facei];
            //const scalar absError = mag(nutw[facei] - expected);

            //const scalar percentError =(mag(expected) > SMALL)? 100.0 * absError / mag(expected): 0.0; // or -1, GREAT, NaN, etc. depending on how you want to flag it


            // LOGS
            // wmlog
            //  << "facei: " << facei
            //  << ", y: " << y[facei]
            //  << ", yT: " << yT
            //  << ", t1: " << t1[facei]
            //  << ", magUp: " << magUp[facei]
            //  << ", magU_stream: " << magU_stream
            //  << ", Uc_new: " << Uc_new[facei]
            //  << ", magUT: " << magUpT
            //  << ", magUT_stream: " << magUT_stream
            //  << ", UT: " << UT
            //  << ", magUp_old: " << magUp_old[facei]
            //  << ", Uc_old: " << Uc_old[facei]
            //  << ", nuw: " << nuw[facei]
            //  << ", nuwT: " << nuwT
            //  <<  ", rho: " << rhop[facei]
            //  << ", yPlus: " << yPlus[facei]
            //  << ", yPlusCalc: " << yPlusCalc
            //  << ", dudt_stream: " << dudt_stream[facei]
            //  << ", dudt_n: " << dudt_n[facei]
            //  << ", dudt_span: " << dudt_span[facei]
            //  << ", meanKE: " << meanKE[facei]
            //  << ", turbKE: " << turbKE[facei]
            //  << ", //INPUTS:"
            //  << " centerRe: " << centerRe
            //  << ", topRe: " << topRe
            //  << ", denormalized dudt_stream: " << dudt_stream_input
            //  << ", denormalized dudt_n: " << dudt_n_input
            //  << ", denormalized dudt_span: " << dudt_span_input
            //  << ", angle_center_top: " << angle_center_top
            //  << ", dtheta_dt: " << dtheta_dt
            //  << ", turbKE/meanKE: " << turbKE[facei]/max(SMALL,meanKE[facei])
            //  << Foam::nl; //debug
            //  << "predict: " << nutw[facei]
            //  << ", expected: " << nutwTest[facei]
            //  << ", error: " << mag(nutw[facei] - nutwTest[facei])
            //  << ", percentError: " << percentError
            //  << ", sigma: " << nutw_sigma << Foam::nl; //debug
        }

        return tnutw;
    }

    tmp<scalarField> unifiedNutUWallFunctionFvPatchScalarField::yPlus(
        const scalarField &magUp) const
    {
        const label patchi = patch().index();

        const momentumTransportModel &turbModel =
            db().lookupType<momentumTransportModel>(internalField().group());

        const scalarField &y = turbModel.yb()[patchi];
        const tmp<scalarField> tnuw = turbModel.nu(patchi);
        const scalarField &nuw = tnuw();

        tmp<scalarField> tyPlus(new scalarField(patch().size(), 0.0));
        scalarField &yPlus = tyPlus.ref();

        forAll(yPlus, facei)
        {
            const scalar Re = magUp[facei] * y[facei] / nuw[facei];
            const scalar ryPlusLam = 1 / yPlusLam_;

            // std::cout  << "facei: " << facei << "Re: " << Re << std::endl; //debug

            int iter = 0;
            scalar yp = yPlusLam_;
            scalar yPlusLast = yp;

            do
            {
                yPlusLast = yp;
                if (yp > yPlusLam_)
                {
                    yp = (kappa_ * Re + yp) / (1 + log(E_ * yp));
                }
                else
                {
                    yp = sqrt(Re);
                }
            } while (mag(ryPlusLam * (yp - yPlusLast)) > 0.0001 && ++iter < 20);

            yPlus[facei] = yp;
        }

        return tyPlus;
    }

    void unifiedNutUWallFunctionFvPatchScalarField::buildSecondLayerAddressing() const
    {
        const fvMesh &mesh = patch().boundaryMesh().mesh();
        const labelUList &faceCells = patch().faceCells();
        const labelListList &cellCells = mesh.cellCells();
        const vectorField &Cc = mesh.C();
        const vectorField nf(patch().nf());

        // const momentumTransportModel& turbModel =
        // db().lookupType<momentumTransportModel>(internalField().group());

        // const volScalarField& yField = turbModel.y();

        secondLayerCell_.setSize(patch().size(), -1);

        forAll(faceCells, facei)
        {
            label ownCell = faceCells[facei];
            const labelList &neighbours = cellCells[ownCell];

            scalar bestAlignment = GREAT;

            label bestCell = -1;

            forAll(neighbours, ni)
            {
                label neiCell = neighbours[ni];
                vector d = Cc[neiCell] - Cc[ownCell];
                scalar alignment = (d & nf[facei]) / (mag(d) + SMALL);

                if (alignment < bestAlignment)
                {
                    bestAlignment = alignment;
                    bestCell = neiCell;
                }
            }

            secondLayerCell_[facei] = bestCell;
        }

        secondLayerCellValid_ = true;
    }

    void unifiedNutUWallFunctionFvPatchScalarField::buildTopFaceAddressing() const
    {
        const fvMesh &mesh = patch().boundaryMesh().mesh();
        const labelUList &faceCells = patch().faceCells();
        const cellList &cells = mesh.cells();
        const vectorField &Sf = mesh.faceAreas();
        const vectorField nf(patch().nf());

        topFace_.setSize(patch().size(), -1);

        forAll(faceCells, facei)
        {
            const label ownCell = faceCells[facei];
            const labelList &cFaces = cells[ownCell];
            const label wallFacei = patch().start() + facei;

            scalar bestAlignment = -GREAT;
            label bestFace = -1;

            forAll(cFaces, fi)
            {
                const label globalFacei = cFaces[fi];

                if (globalFacei == wallFacei)
                {
                    continue;
                }

                vector faceN = Sf[globalFacei] / (mag(Sf[globalFacei]) + SMALL);

                // orient outward from ownCell, consistent with patch nf
                if (mesh.owner()[globalFacei] != ownCell)
                {
                    faceN = -faceN;
                }

                const scalar alignment = faceN & (-nf[facei]);

                if (alignment > bestAlignment)
                {
                    bestAlignment = alignment;
                    bestFace = globalFacei;
                }
            }

            topFace_[facei] = bestFace;
        }

        topFaceValid_ = true;
    }

    scalar unifiedNutUWallFunctionFvPatchScalarField::topFaceValue(
        const volScalarField &field,
        const label facei) const
    {
        const fvMesh &mesh = patch().boundaryMesh().mesh();
        const label fFace = topFace_[facei];

        if (fFace < mesh.nInternalFaces())
        {
            const label own = mesh.owner()[fFace];
            const label nei = mesh.neighbour()[fFace];
            const scalar w = mesh.weights()[fFace];
            return w * field[own] + (1 - w) * field[nei];
        }
        else
        {
            const label patchi = mesh.boundaryMesh().whichPatch(fFace);
            const label pFacei = fFace - mesh.boundaryMesh()[patchi].start();
            return field.boundaryField()[patchi][pFacei];
        }
    }

    vector unifiedNutUWallFunctionFvPatchScalarField::topFaceValueVector(
        const volVectorField &field,
        const label facei) const
    {
        const fvMesh &mesh = patch().boundaryMesh().mesh();
        const label fFace = topFace_[facei];

        if (fFace < mesh.nInternalFaces())
        {
            const label own = mesh.owner()[fFace];
            const label nei = mesh.neighbour()[fFace];
            const scalar w = mesh.weights()[fFace];
            return w * field[own] + (1 - w) * field[nei];
        }
        else
        {
            const label patchi = mesh.boundaryMesh().whichPatch(fFace);
            const label pFacei = fFace - mesh.boundaryMesh()[patchi].start();
            return field.boundaryField()[patchi][pFacei];
        }
    }

    void unifiedNutUWallFunctionFvPatchScalarField::ensureSurfaceCoordinate() const
    {
        if (sValid_)
        {
            return;
        }

        const fvMesh &mesh = patch().boundaryMesh().mesh();

        if (mesh.foundObject<volScalarField>("chordwisePotential"))
        {
            sPtr_.reset(
                new volScalarField(
                    mesh.lookupObject<volScalarField>("chordwisePotential")));
            sValid_ = true;
            return;
        }

        // Build s from scratch: Laplace equation, s = x on farfield,
        // zeroGradient on wall(s).
        sPtr_.reset(
            new volScalarField(
                IOobject(
                    "chordwisePotential",
                    mesh.time().timeName(mesh.time().value()),
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE),
                mesh,
                dimensionedScalar(dimless, 0),
                "zeroGradient"));

        volScalarField &s = sPtr_();

        // Now overwrite the farfield patch(es) with fixedValue = local x
        // Loop over all patches; treat anything that's a wall as
        // zeroGradient (already the default), everything else (farfield,
        // inlet, outlet, symmetry-ish far boundary) as fixedValue = x.
        forAll(mesh.boundary(), patchi)
        {
            const polyPatch &pp = mesh.boundaryMesh()[patchi];

            const bool isWall = isA<wallPolyPatch>(pp);
            const bool isCoupled = pp.coupled(); // processor, cyclic, cyclicAMI
            const bool isConstraint = polyPatch::constraintType(pp.type());
            // empty, symmetry, wedge, etc.

            if (!isWall && !isCoupled && !isConstraint)
            {
                // genuine external farfield/inlet/outlet patch
                const vectorField &Cp = mesh.boundary()[patchi].Cf();
                scalarField xVals(Cp.component(vector::X));

                s.boundaryFieldRef().set(
                    patchi,
                    fvPatchField<scalar>::New(
                        "fixedValue",
                        mesh.boundary()[patchi],
                        s));

                fixedValueFvPatchScalarField &fvp =
                    refCast<fixedValueFvPatchScalarField>(
                        s.boundaryFieldRef()[patchi]);

                fvp == xVals;
            }
            // else: leave at whatever the field's default constructor gave it
            // (zeroGradient) - correct behavior for wall AND for coupled/constraint
            // patches, since those are handled by their own patch field types
            // regardless of what we assign here.
        }

        dictionary solverDict;
        solverDict.add("solver", "GAMG");
        solverDict.add("tolerance", 1e-8);
        solverDict.add("relTol", 0.01);
        solverDict.add("smoother", "GaussSeidel");

        // Solve Laplace's equation for s
        fvScalarMatrix sEqn(fvm::laplacian(s));
        sEqn.solve(solverDict);

        s.correctBoundaryConditions();
        s.write();

        sValid_ = true;
    }

    // * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

    unifiedNutUWallFunctionFvPatchScalarField::unifiedNutUWallFunctionFvPatchScalarField(
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const dictionary &dict)
        : nutWallFunctionFvPatchScalarField(p, iF, dict),
          // mlModel_(new wallModelMLP(this->db().time().constant()/"wallModel")),
          // ensembleModel_(new crossflowWallModel(this->db().time().constant()/"wallModel_pytorch")),
          // utauMlModel_(new wallModelMLP(this->db().time().constant()/"utauWallModel")),
          //  Directory renamed from "wallModel_classifier" to "wallModel_unified" -
          //  this now holds the single unified-network export (wmLayers_openfoam/
          //  unified/ + cz_uncertainty.txt), not the old classifier+4-regime-model
          //  layout. classifierModel_ itself is left named as-is (private, internal,
          //  still an autoPtr<unifiedWallModel> with the same external interface).
          classifierModel_(new unifiedWallModel(this->db().time().constant() / "wallModel_full" / "wallModel_unified")),
          secondLayerCell_(),
          secondLayerCellValid_(false),
          topFace_(),
          topFaceValid_(false)

    {
    }

    unifiedNutUWallFunctionFvPatchScalarField::unifiedNutUWallFunctionFvPatchScalarField(
        const unifiedNutUWallFunctionFvPatchScalarField &ptf,
        const fvPatch &p,
        const DimensionedField<scalar, volMesh> &iF,
        const fieldMapper &mapper)
        : nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper),
          classifierModel_(new unifiedWallModel(this->db().time().constant() / "wallModel_full" / "wallModel_unified")),
          secondLayerCell_(),
          secondLayerCellValid_(false),
          topFace_(),
          topFaceValid_(false)
    {
    }

    unifiedNutUWallFunctionFvPatchScalarField::unifiedNutUWallFunctionFvPatchScalarField(
        const unifiedNutUWallFunctionFvPatchScalarField &sawfpsf,
        const DimensionedField<scalar, volMesh> &iF)
        : nutWallFunctionFvPatchScalarField(sawfpsf, iF),
          classifierModel_(new unifiedWallModel(this->db().time().constant() / "wallModel_full" / "wallModel_unified")),
          secondLayerCell_(),
          secondLayerCellValid_(false),
          topFace_(),
          topFaceValid_(false)
    {
    }

    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    tmp<scalarField> unifiedNutUWallFunctionFvPatchScalarField::yPlus() const
    {
        const label patchi = patch().index();

        const momentumTransportModel &turbModel =
            db().lookupType<momentumTransportModel>(internalField().group());

        const fvPatchVectorField &Uw = turbModel.U().boundaryField()[patchi];
        const scalarField magUp(mag(Uw.patchInternalField() - Uw));

        return yPlus(magUp);
    }

    void unifiedNutUWallFunctionFvPatchScalarField::write(Ostream &os) const
    {
        fvPatchField<scalar>::write(os);
        writeLocalEntries(os);
        writeEntry(os, "value", *this);
    }

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

    makePatchTypeField(
        fvPatchScalarField,
        unifiedNutUWallFunctionFvPatchScalarField);

    // * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
