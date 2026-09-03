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

#include "mlNutUWallFunctionFvPatchScalarField.H"
#include "momentumTransportModel.H"
#include "fieldMapper.H"
#include "volFields.H"
#include "addToRunTimeSelectionTable.H"
#include "OFstream.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

namespace
{
    void writeTwField
    (
        const fvMesh& mesh,
        const label patchi,
        const label facei,
        const scalar value,
        const word& fieldName
    )
    {
        if (!mesh.foundObject<volScalarField>(fieldName))
        {
            volScalarField* fieldPtr = new volScalarField
            (
                IOobject
                (
                    fieldName,
#ifdef OPENFOAM
                    mesh.time().timeName(),
#else
                    mesh.time().name(),
#endif
                    mesh,
                    IOobject::NO_READ,
                    IOobject::AUTO_WRITE
                ),
                mesh,
                dimensionedScalar(fieldName, dimless, 0.0)
            );
            fieldPtr->store();
        }

        volScalarField& twField =
            const_cast<volScalarField&>
            (
                mesh.lookupObject<volScalarField>(fieldName)
            );
        twField.boundaryFieldRef()[patchi][facei] = value;
    }

    const std::array<word, WM_N_REGIMES> regimeWeightMeanNames =
    {
        word("classifier_w_laminar_mean"),
        word("classifier_w_freestream_mean"),
        word("classifier_w_channel_mean"),
        word("classifier_w_couette_mean"),
        word("classifier_w_crossflow_mean")
    };
    const std::array<word, WM_N_REGIMES> regimeWeightSigmaNames =
    {
        word("classifier_w_laminar_sigma"),
        word("classifier_w_freestream_sigma"),
        word("classifier_w_channel_sigma"),
        word("classifier_w_couette_sigma"),
        word("classifier_w_crossflow_sigma")
    };

    void writeClassifierWeightFields
    (
        const fvMesh& mesh,
        const label patchi,
        const label facei,
        const WeightVector& weightsMean,
        const WeightVector& weightsSigma
    )
    {
        for (int r = 0; r < WM_N_REGIMES; ++r)
        {
            writeTwField(mesh, patchi, facei, weightsMean[r], regimeWeightMeanNames[r]);
            writeTwField(mesh, patchi, facei, weightsSigma[r], regimeWeightSigmaNames[r]);
        }
    }
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

tmp<scalarField> mlNutUWallFunctionFvPatchScalarField::nut() const
{
    const label patchi = patch().index();

    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    const volVectorField& Ufield = turbModel.U(); //for top layer

    const fvPatchVectorField& Uw = Ufield.boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    const tmp<volScalarField> tnuField = turbModel.nu();
    const volScalarField& nuField = tnuField();

    const volScalarField& yField = turbModel.y();
    const scalarField y(yField.boundaryField()[patchi].patchInternalField());

    tmp<scalarField> tnutw(new scalarField(patch().size(), 0.0));
    scalarField& nutw = tnutw.ref();

    tmp<scalarField> tnutwTest(new scalarField(patch().size(), 0.0));
    scalarField& nutwTest = tnutwTest.ref();

    //compute acceleration

    const fvMesh& mesh = patch().boundaryMesh().mesh();

    //const volVectorField& U = mesh.lookupObject<volVectorField>("U");
    const fvPatchVectorField& Uw_old = Ufield.oldTime().boundaryField()[patchi];
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

    // Directions computed from flow
    //vectorField t1(Uc_new - (Uc_new & nf) * nf);
    //t1 /= (mag(t1) + SMALL);
    //const vectorField t2(nf ^ t1);

    // Basis computed from geometry

    // Free-stream / reference direction
    const vector chordDir(1, 0, 0);  

    // Project chord direction onto the local wall-tangent plane at each face
    vectorField t1(chordDir - (chordDir & nf) * nf);
    t1 /= (mag(t1) + SMALL);

    const vectorField t2(nf ^ t1);

    //Project acceleration onto the basis 

    scalarField dudt_n(dudt_local & (-nf));       // wall-normal
    scalarField dudt_stream(dudt_local & t1);  // streamwise
    scalarField dudt_span(dudt_local & t2);    // spanwise

    //end compute acceleration

    //compute kinetic energy (RANS)

    scalarField meanKE(0.5 * magSqr(Uc_new));   // 0.5 * (U·U), per cell

    //For LES
    scalarField turbKE(patch().size(),0.0);

	//Subgrid scale TKE 
	tmp<volScalarField> tkSgs = turbModel.k();  // valid for all LES models, RANS models too
	const volScalarField& kSgs = tkSgs();
	scalarField kSgsPatch(kSgs.boundaryField()[patchi].patchInternalField());

	//Resolved TKE
	scalarField kResolvedPatch(patch().size(), 0.0);

	if (turbModel.type() != "RAS" && mesh.foundObject<volSymmTensorField>("UPrime2Mean")) //deactivate for RANS
	{
		const volSymmTensorField& UPrime2Mean =
			mesh.lookupObject<volSymmTensorField>("UPrime2Mean");

		const symmTensorField UPrime2MeanPatch
		(
			UPrime2Mean.boundaryField()[patchi].patchInternalField()
		);

		kResolvedPatch = 0.5*tr(UPrime2MeanPatch);
	}

    turbKE = kSgsPatch + kResolvedPatch; //SGS + Resolved

    //end compute kinetic energy

    const scalarField yPlus(this->yPlus(magUp));

    //const scalar Aplus = 26.0; //For LES wall function, not used in RANS

    static Foam::OFstream wmlog("nutwmlog.log");

    forAll(yPlus, facei)
    {
        if (yPlus[facei] > yPlusLam_)
        {
           nutwTest[facei] = nuw[facei]*(yPlus[facei]*kappa_/log(E_*yPlus[facei]) - 1);
        }
    }

    //wallModelMLP nutMlModel(this->db().time().constant()/"wallModel");

    //if (!secondLayerCellValid_)
    //{
    //    buildSecondLayerAddressing();
    //}

    if (!topFaceValid_)
    {
        buildTopFaceAddressing();
    }

    std::cout << "nutU wall function loop" << std::endl; //debug

    forAll(nutw, facei)
    {
        //Top face values

        const scalar yT = topFaceValue(yField, facei);
        const vector UT = topFaceValueVector(Ufield, facei);
        const scalar magUpT = mag(UT - Uw[facei]);
        const scalar nuwT = topFaceValue(nuField, facei);

        //Streamwise magnitudes

        const scalar magU_stream = mag(Uc_new[facei] & t1[facei]);
        const scalar magUT_stream = mag(UT & t1[facei]);

        //spanwise magnitudes
        const scalar magU_span = mag(Uc_new[facei] & t2[facei]);
        const scalar magUT_span = mag(UT & t2[facei]);

        //Streamwise Reynolds

        const scalar centerRe = y[facei]*sqrt(magU_stream*magU_stream+magU_span*magU_span)/nuw[facei];
        const scalar topRe = yT*sqrt(magU_stream*magU_stream+magU_span*magU_span)/nuwT;

        //Angles
        const scalar angle_center_top = std::remainder(atan2(UT & t2[facei], UT & t1[facei]) - atan2(Uc_new[facei] & t2[facei], Uc_new[facei] & t1[facei]), 2 * M_PI);
        const scalar dtheta_dt = (std::remainder((atan2(Uc_new[facei] & t2[facei], Uc_new[facei] & t1[facei])-atan2(Uc_old[facei] & t2[facei], Uc_old[facei] & t1[facei])),2*M_PI)/deltaT)*sqr(y[facei])/nuw[facei];
        
        //Inputs

        const scalar dudt_stream_input = dudt_stream[facei]*pow(y[facei], 3)/pow(nuw[facei], 2);

        std::vector<scalar> inputs = {centerRe, topRe, dudt_stream_input, angle_center_top, dtheta_dt, turbKE[facei]/max(SMALL,meanKE[facei])};

        unifiedWallModel::Prediction values = classifierModel_->predict(inputs);

        const scalar tw_n_prediction = max(0, values.twMean);

        const scalar tw_n_sigma = values.twSigma;
        const scalar nutw_sigma =
            (values.twMean > 0)
          ? nuw[facei] * (magU_stream / max(SMALL, magUp[facei])) * tw_n_sigma
          : 0;

        writeTwField(mesh, patchi, facei, tw_n_sigma, "tw_n_sigma");
        writeTwField(mesh, patchi, facei, tw_n_prediction, "tw_n_mean");

        writeClassifierWeightFields(mesh, patchi, facei, values.weightsMean, values.weightsSigma);



        //std::cout << "utau flag" << nuw[facei] *tw_n_prediction * magU_stream / y[facei] << std::endl; //debug
        //Info << "Reached utau" << endl;
        const scalar uTau = sqrt(nuw[facei] *tw_n_prediction * magU_stream / max(SMALL, y[facei]));

        //const scalar uTau = ensembleModel_->predict(inputs)[0]*magUp[facei]*nuw[facei]/(y[facei]*rhop[facei]);
        const scalar yPlusCalc = uTau*y[facei]/max(SMALL, nuw[facei]);

        //LOGS
        /*
        wmlog 
        << "facei: " << facei 
        << ", y: " << y[facei] 
        << ", yT: " << yT 
        << ", t1: " << t1[facei]
        << ", magUp: " << magUp[facei] 
        << ", magUT: " << magUpT
        << ", yPlusCalc: " << yPlusCalc 
        << ", dudt_stream: " << dudt_stream[facei]
        << ", dudt_n: " << dudt_n[facei]
        << ", dudt_span: " << dudt_span[facei]
        << ", meanKE: " << meanKE[facei]
        << ", turbKE: " << turbKE[facei]
        << ", //INPUTS:"
        << " centerRe: " << centerRe
        << ", topRe: " << topRe
        << ", denormalized dudt_stream: " << dudt_stream_input
        << ", angle_center_top: " << angle_center_top
        << ", dtheta_dt: " << dtheta_dt
        << ", turbKE/meanKE: " << turbKE[facei]/max(SMALL,meanKE[facei])
        << Foam::nl; //debug
        */

        if (turbModel.type() != "RAS") //deactivate for RANS
	    {
            nutw[facei] = max(0, (uTau*uTau*y[facei])/max(SMALL, magUp[facei]) - nuw[facei]); //RANS
        }
        else 
        {
            nutw[facei] = kappa_*uTau*y[facei]*(1-std::exp(-yPlusCalc/Aplus))*(1-std::exp(-yPlusCalc/Aplus)); //LES
        }
            //Info << "Reached nutw" << endl;
        

        const scalar expected = nutwTest[facei];
        const scalar absError = mag(nutw[facei] - expected);

        const scalar percentError =
        (mag(expected) > SMALL)
        ? 100.0*absError/mag(expected)
        : 0.0;   // or -1, GREAT, NaN, etc. depending on how you want to flag it

        /*
        wmlog 
        << "predict: " << nutw[facei] 
        << ", expected: " << nutwTest[facei] 
        << ", error: " << mag(nutw[facei] - nutwTest[facei]) 
        << ", percentError: " << percentError
        << ", sigma: " << nutw_sigma << Foam::nl; //debug
        */
    }

    return tnutw;
}


tmp<scalarField> mlNutUWallFunctionFvPatchScalarField::yPlus
(
    const scalarField& magUp
) const
{
    const label patchi = patch().index();

    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const scalarField& y = turbModel.yb()[patchi];
    const tmp<scalarField> tnuw = turbModel.nu(patchi);
    const scalarField& nuw = tnuw();

    tmp<scalarField> tyPlus(new scalarField(patch().size(), 0.0));
    scalarField& yPlus = tyPlus.ref();

    forAll(yPlus, facei)
    {
        const scalar Re = magUp[facei]*y[facei]/nuw[facei];
        const scalar ryPlusLam = 1/yPlusLam_;

        //std::cout  << "facei: " << facei << "Re: " << Re << std::endl; //debug

        int iter = 0;
        scalar yp = yPlusLam_;
        scalar yPlusLast = yp;

        do
        {
            yPlusLast = yp;
            if (yp > yPlusLam_)
            {
                yp = (kappa_*Re + yp)/(1 + log(E_*yp));
            }
            else
            {
                yp = sqrt(Re);
            }
        } while(mag(ryPlusLam*(yp - yPlusLast)) > 0.0001 && ++iter < 20);

        yPlus[facei] = yp;
    }

    return tyPlus;
}

void mlNutUWallFunctionFvPatchScalarField::buildSecondLayerAddressing() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const labelUList& faceCells = patch().faceCells();
    const labelListList& cellCells = mesh.cellCells();
    const vectorField& Cc = mesh.C();
    const vectorField nf(patch().nf());

    secondLayerCell_.setSize(patch().size(), -1);

    forAll(faceCells, facei)
    {
        label ownCell = faceCells[facei];
        const labelList& neighbours = cellCells[ownCell];

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

void mlNutUWallFunctionFvPatchScalarField::buildTopFaceAddressing() const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const labelUList& faceCells = patch().faceCells();
    const cellList& cells = mesh.cells();
    const vectorField& Sf = mesh.faceAreas();
    const vectorField nf(patch().nf());

    topFace_.setSize(patch().size(), -1);

    forAll(faceCells, facei)
    {
        const label ownCell = faceCells[facei];
        const labelList& cFaces = cells[ownCell];
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

scalar mlNutUWallFunctionFvPatchScalarField::topFaceValue
(
    const volScalarField& field,
    const label facei
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const label fFace = topFace_[facei];

    if (fFace < mesh.nInternalFaces())
    {
        const label own = mesh.owner()[fFace];
        const label nei = mesh.neighbour()[fFace];
        const scalar w = mesh.weights()[fFace];
        return w*field[own] + (1 - w)*field[nei];
    }
    else
    {
        const label patchi = mesh.boundaryMesh().whichPatch(fFace);
        const label pFacei = fFace - mesh.boundaryMesh()[patchi].start();
        return field.boundaryField()[patchi][pFacei];
    }
}

vector mlNutUWallFunctionFvPatchScalarField::topFaceValueVector
(
    const volVectorField& field,
    const label facei
) const
{
    const fvMesh& mesh = patch().boundaryMesh().mesh();
    const label fFace = topFace_[facei];

    if (fFace < mesh.nInternalFaces())
    {
        const label own = mesh.owner()[fFace];
        const label nei = mesh.neighbour()[fFace];
        const scalar w = mesh.weights()[fFace];
        return w*field[own] + (1 - w)*field[nei];
    }
    else
    {
        const label patchi = mesh.boundaryMesh().whichPatch(fFace);
        const label pFacei = fFace - mesh.boundaryMesh()[patchi].start();
        return field.boundaryField()[patchi][pFacei];
    }
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

mlNutUWallFunctionFvPatchScalarField::mlNutUWallFunctionFvPatchScalarField
(
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const dictionary& dict
)
:
    nutWallFunctionFvPatchScalarField(p, iF, dict),
    classifierModel_(new unifiedWallModel(this->db().time().constant()/"wallModel_classifier")),
    secondLayerCell_(),
    secondLayerCellValid_(false),
    topFace_(),
    topFaceValid_(false) 
{}

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


mlNutUWallFunctionFvPatchScalarField::mlNutUWallFunctionFvPatchScalarField
(
    const mlNutUWallFunctionFvPatchScalarField& ptf,
    const fvPatch& p,
    const DimensionedField<scalar, volMesh>& iF,
    const fieldMapper& mapper
)
:
    nutWallFunctionFvPatchScalarField(ptf, p, iF, mapper)
{}


mlNutUWallFunctionFvPatchScalarField::mlNutUWallFunctionFvPatchScalarField
(
    const mlNutUWallFunctionFvPatchScalarField& sawfpsf,
    const DimensionedField<scalar, volMesh>& iF
)
:
    nutWallFunctionFvPatchScalarField(sawfpsf, iF)
{}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<scalarField> mlNutUWallFunctionFvPatchScalarField::yPlus() const
{
    const label patchi = patch().index();


    const momentumTransportModel& turbModel =
        db().lookupType<momentumTransportModel>(internalField().group());

    const fvPatchVectorField& Uw = turbModel.U().boundaryField()[patchi];
    const scalarField magUp(mag(Uw.patchInternalField() - Uw));

    return yPlus(magUp);
}


void mlNutUWallFunctionFvPatchScalarField::write(Ostream& os) const
{
    fvPatchField<scalar>::write(os);
    writeLocalEntries(os);
    writeEntry(os, "value", *this);
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

makePatchTypeField
(
    fvPatchScalarField,
    mlNutUWallFunctionFvPatchScalarField
);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// ************************************************************************* //
