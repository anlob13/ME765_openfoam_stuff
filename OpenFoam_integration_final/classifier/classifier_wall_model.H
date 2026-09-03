#ifndef classifier_wall_model_H
#define classifier_wall_model_H

#include <array>
#include <string>
#include <utility>
#include <vector>
#include "scalar.H"
#include "fileName.H"
#include "regime_wall_model.H"
#include "classifier_ensemble.H"

namespace Foam
{

// Regime label order - fixed everywhere in this folder (C++ and the Python
// training/export side alike): classifier_building_block/train_classifier.py's
// own LABEL_NAMES list.
enum RegimeLabel
{
    REGIME_LAMINAR = 0,
    REGIME_FREESTREAM = 1,
    REGIME_CHANNEL = 2,
    REGIME_COUETTE = 3,
    REGIME_CROSSFLOW = 4
};

// The full "evaluate the flow" model: classifies which of the 5 regimes a
// local point resembles (with genuine, simplex-summing-to-1 uncertainty on
// that classification), runs all 4 DISTINCT trained building-block models on
// the same input, and blends them into one (tw, tw_direction_angle) + UQ
// prediction, propagating BOTH each block's own uncertainty AND the
// classifier's uncertainty about which block even applies.
//
// "Freestream" has no separately-trained model (per the project's own
// framing: freestream is channel flow's own low-turbulence outer region, not
// a distinct physical regime with its own simulation) - its classifier label
// aliases the channel regimeWallModel's prediction directly. See
// classifier_building_block/README.md section 1 for the full reasoning.
//
// Combination formula (derived from the law of total variance conditioned on
// the random weight vector W; verified against a 5,000,000-sample Monte Carlo
// simulation in classifier_building_block/verify_combination.py before being
// trusted here - same formula, Python side):
//     E[Y]   = sum_i w_bar_i * mu_i
//     Var(Y) = sum_i (Var(W_i) + w_bar_i^2) * sigma_i^2      <- aleatoric,
//                                                                inflated by
//                                                                weight
//                                                                uncertainty
//            + mu^T * Sigma_W * mu                            <- "model
//                                                                selection"
//                                                                uncertainty,
//                                                                large exactly
//                                                                when the
//                                                                classifier is
//                                                                unsure AND the
//                                                                candidate
//                                                                blocks
//                                                                disagree
class classifierWallModel
{
private:

    classifierEnsemble classifier_;
    regimeWallModel laminarModel_;
    regimeWallModel channelModel_;
    regimeWallModel couetteModel_;
    regimeWallModel crossflowModel_;

    // (C,Z)-calibration-driven epistemic uncertainty on tw, ADDED to each
    // regime's own sigma2 (laminarPred.twVar etc.) before combine() runs -
    // see Uncertainty_propagation/README.md for the full derivation and
    // ../Uncertainty_propagation/anova_reducer.py for how these two numbers
    // were produced. Physically: "how much would this regime's trained
    // model's mean prediction have differed had ODT's (C,Z) eddy-model
    // constants been drawn from a different, still-plausible point in the
    // calibration's own posterior, instead of the single best-fit point
    // actually used to generate its training data" - a source of
    // uncertainty the NN ensembles cannot see on their own, since each one
    // only ever saw ONE (C,Z)-generated dataset. Loaded once at
    // construction from wmLayersDir/"cz_uncertainty.txt" (optional - MISSING
    // file defaults both to 0, i.e. behaves exactly like every
    // wmLayers_openfoam/ export that predates this feature - same
    // backward-compatible pattern angle_type.txt already established).
    //
    // Regime mapping (applied in predict(), not stored per-regime here -
    // laminar gets NEITHER value, it has no ODT/(C,Z) dependence at all,
    // being a closed-form analytical solution):
    //   channel, freestream (aliases channel), crossflow (shares channel's
    //     calibrated C,Z - confirmed via each building block's own
    //     generator.py) -> tau2Channel_
    //   couette                                                -> tau2Couette_
    scalar tau2Channel_;
    scalar tau2Couette_;

    static std::pair<scalar, scalar> combine
    (
        const WeightVector& weightsMean,
        const WeightCovMatrix& weightsCov,
        const std::array<scalar, WM_N_REGIMES>& mu,
        const std::array<scalar, WM_N_REGIMES>& sigma2
    );

    // tw_direction_angle-only: circular-resultant-vector combination, NOT
    // the linear combine() above - see classifier_building_block's own
    // blend_predict.py module docstring and verify_circular_combination.py
    // for the full derivation and Monte Carlo verification. A real, measured
    // bug (not just a theoretical risk): the linear formula gave ~pi errors
    // on real Couette data at extreme pressure differentials, whenever two
    // significant-weight regimes' angle predictions straddled the +-pi
    // wraparound. Deliberately takes only weightsMean, not weightsCov -
    // unlike combine()'s own Var(Y) (a second-moment quantity that genuinely
    // needs the weight covariance via the law of total variance), the
    // resultant vector is a first-moment quantity and by linearity of
    // expectation only needs the MEAN weights, regardless of Var(W) - proven
    // both analytically and empirically (see verify_circular_combination.py).
    static std::pair<scalar, scalar> combineCircular
    (
        const WeightVector& weightsMean,
        const std::array<scalar, WM_N_REGIMES>& mu,
        const std::array<scalar, WM_N_REGIMES>& sigma2
    );

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

    // dir should contain wmLayers_openfoam/{laminar,channel,couette,crossflow}/
    // (each a regimeWallModel directory - tw/, tw_direction_angle/,
    // input_metadata.txt) and wmLayers_openfoam/classifier/ (a
    // classifierEnsemble directory - member_0/../member_4/, input_metadata.txt).
    explicit classifierWallModel(const fileName& wmLayersDir);

    // rawInputs: (Re_center, Re_top, ax_n, gamma12, lambda_roc, tke_n), raw
    // physical values, in that column order - the same 6 inputs every
    // building block and the classifier itself already use. Each model's own
    // input transform (log1p/arcsinh) happens internally; feed raw values only.
    Prediction predict(const std::vector<scalar>& rawInputs) const;
};

}

#endif
