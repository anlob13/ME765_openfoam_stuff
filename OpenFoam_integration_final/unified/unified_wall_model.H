#ifndef unified_wall_model_H
#define unified_wall_model_H

#include <array>
#include <string>
#include <utility>
#include <vector>
#include "scalar.H"
#include "fileName.H"
#include "unified_ensemble.H"

namespace Foam
{

// Regime label order - fixed everywhere in this folder (C++ and the Python
// training/export side alike): ../neural_network_v4/collect_training_data.py's
// own LABEL_NAMES list (unchanged from the original classifier's own order).
enum RegimeLabel
{
    REGIME_LAMINAR = 0,
    REGIME_FREESTREAM = 1,
    REGIME_CHANNEL = 2,
    REGIME_COUETTE = 3,
    REGIME_CROSSFLOW = 4
};

// The full "evaluate the flow" model - SUPERSEDES the original two-stage
// design (a separately-trained classifier blending 4 separately-trained
// per-regime regression ensembles via combine()/combineCircular(), see this
// project's OpenFoam_integration_final/unified_wall_model.C for that version).
//
// Here, ONE shared unifiedEnsemble (see unified_ensemble.H) already produces
// (weightsMean, weightsCov, twMean, twVar, angleMean, angleVar) directly, in a
// single forward pass per member - there is no separate per-regime blend left
// to do here. This class's only remaining job is adding the (C,Z)-calibration
// epistemic uncertainty term on top of the ensemble's own tw variance (see
// tau2Channel_/tau2Couette_ below) - everything else is a thin pass-through.
//
// Same external interface (class name, Prediction fields, constructor
// signature) as the design this replaces, so mlNutUWallFunctionFvPatchScalarField.C
// needed no changes beyond which subdirectory name it points at.
class unifiedWallModel
{
private:

    unifiedEnsemble ensemble_;

    // (C,Z)-calibration-driven epistemic uncertainty on tw - see
    // Uncertainty_propagation/README.md for the full derivation. Loaded once
    // at construction from wmLayersDir/"cz_uncertainty.txt" (optional -
    // missing file defaults both to 0).
    //
    // UNLIKE the two-stage design this replaces (which added tau2 to each
    // regime's OWN sigma2 before blending), there is only one tw output here
    // - so the (C,Z) term is blended by the classification weights instead,
    // exactly mirroring how those weights would have selected a regime's tau2
    // in the old design:
    //   tau2_effective = w_freestream*tau2Channel_ + w_channel*tau2Channel_
    //                   + w_crossflow*tau2Channel_ + w_couette*tau2Couette_
    //                   (+ w_laminar*0 - laminar has no ODT/(C,Z) dependence)
    // This still uses the classification output, but only as a soft indicator
    // of which calibration uncertainty applies here - NOT to gate or reweight
    // the regression prediction itself, preserving this design's core
    // decoupling property (see unified_ensemble.H's own module comment).
    scalar tau2Channel_;
    scalar tau2Couette_;

public:

    struct Prediction
    {
        scalar twMean;
        scalar twSigma;
        scalar angleMean;
        scalar angleSigma;
        WeightVector weightsMean;   // label order: see RegimeLabel above
        WeightVector weightsSigma;  // marginal per-regime weight uncertainty
        WeightCovMatrix weightsCov; // full covariance, in case a caller wants
                                    // more than the marginals
    };

    // dir should contain wmLayers_openfoam/unified/ (a unifiedEnsemble
    // directory - member_0/../member_4/, input_metadata.txt, output_metadata.txt)
    // and optionally wmLayers_openfoam/cz_uncertainty.txt.
    explicit unifiedWallModel(const fileName& wmLayersDir);

    // rawInputs: (Re_center, Re_top, ax_n, gamma12, lambda_roc, tke_n), raw
    // physical values, in that column order - unchanged from the original
    // design. The input transform (log1p/arcsinh) happens internally.
    Prediction predict(const std::vector<scalar>& rawInputs) const;
};

}

#endif
