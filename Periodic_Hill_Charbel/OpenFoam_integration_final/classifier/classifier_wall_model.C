#include <cmath>
#include <fstream>
#include <sstream>
#include <tuple>
#include <utility>

#include "classifier_wall_model.H"

namespace Foam
{

namespace
{
    // Optional file, two "key value" lines (tau2_channel, tau2_couette) -
    // see classifier_wall_model.H's own comment on tau2Channel_/tau2Couette_
    // for what these mean and why missing = (0,0) is the correct,
    // backward-compatible default (every wmLayers_openfoam/ export that
    // predates this feature has no such file). Deliberately file-local
    // (matches regime_wall_model.C's own readAngleType/readKeyValueFile
    // pattern of small anonymous-namespace metadata readers, not shared
    // via the header).
    std::pair<scalar, scalar> readCZUncertainty(const fileName& wmLayersDir)
    {
        const fileName path = wmLayersDir/"cz_uncertainty.txt";
        std::ifstream f(path);
        if (!f)
        {
            return {scalar(0.0), scalar(0.0)};
        }

        scalar tau2Channel = 0.0;
        scalar tau2Couette = 0.0;
        std::string line;
        while (std::getline(f, line))
        {
            std::istringstream iss(line);
            std::string key;
            double value;
            if (!(iss >> key >> value))
            {
                continue;
            }
            if (key == "tau2_channel")
            {
                tau2Channel = value;
            }
            else if (key == "tau2_couette")
            {
                tau2Couette = value;
            }
        }
        return {tau2Channel, tau2Couette};
    }
}

classifierWallModel::classifierWallModel(const fileName& wmLayersDir)
:
    classifier_(wmLayersDir/"classifier"),
    laminarModel_(wmLayersDir/"laminar"),
    channelModel_(wmLayersDir/"channel"),
    couetteModel_(wmLayersDir/"couette"),
    crossflowModel_(wmLayersDir/"crossflow"),
    tau2Channel_(0.0),
    tau2Couette_(0.0)
{
    std::tie(tau2Channel_, tau2Couette_) = readCZUncertainty(wmLayersDir);
}

std::pair<scalar, scalar> classifierWallModel::combine
(
    const WeightVector& weightsMean,
    const WeightCovMatrix& weightsCov,
    const std::array<scalar, WM_N_REGIMES>& mu,
    const std::array<scalar, WM_N_REGIMES>& sigma2
)
{
    scalar combinedMean = 0.0;
    for (int i = 0; i < WM_N_REGIMES; ++i)
    {
        combinedMean += weightsMean[i] * mu[i];
    }

    // term A: sum_i (Var(W_i) + w_i^2) * sigma_i^2
    scalar termA = 0.0;
    for (int i = 0; i < WM_N_REGIMES; ++i)
    {
        const scalar varWi = weightsCov[i][i];
        termA += (varWi + weightsMean[i]*weightsMean[i]) * sigma2[i];
    }

    // term B: mu^T * Sigma_W * mu (quadratic form - "model selection"
    // uncertainty, large exactly when the classifier is unsure AND the
    // candidate blocks disagree)
    scalar termB = 0.0;
    for (int i = 0; i < WM_N_REGIMES; ++i)
    {
        for (int j = 0; j < WM_N_REGIMES; ++j)
        {
            termB += mu[i] * mu[j] * weightsCov[i][j];
        }
    }

    const scalar combinedVar = termA + termB;
    return {combinedMean, combinedVar};
}

classifierWallModel::Prediction classifierWallModel::predict
(
    const std::vector<scalar>& rawInputs
) const
{
    // Not a structured binding (auto [a, b] = ...) deliberately - that's a
    // C++17 feature, and this project (matching OpenFOAM's own build
    // standard) targets C++14. GCC accepts it anyway as an extension, with a
    // warning, but relying on a non-standard extension here would be fragile.
    const std::pair<WeightVector, WeightCovMatrix> classifierResult = classifier_.predict(rawInputs);
    const WeightVector& weightsMean = classifierResult.first;
    const WeightCovMatrix& weightsCov = classifierResult.second;

    const regimeWallModel::Prediction laminarPred = laminarModel_.predict(rawInputs);
    const regimeWallModel::Prediction channelPred = channelModel_.predict(rawInputs);
    const regimeWallModel::Prediction couettePred = couetteModel_.predict(rawInputs);
    const regimeWallModel::Prediction crossflowPred = crossflowModel_.predict(rawInputs);
    // freestream has no separately-trained model - aliases channel's
    // prediction directly (see the .H file's class comment).
    const regimeWallModel::Prediction& freestreamPred = channelPred;

    std::array<scalar, WM_N_REGIMES> muTw{};
    std::array<scalar, WM_N_REGIMES> sigma2Tw{};
    std::array<scalar, WM_N_REGIMES> muAngle{};
    std::array<scalar, WM_N_REGIMES> sigma2Angle{};

    muTw[REGIME_LAMINAR] = laminarPred.twMean;      sigma2Tw[REGIME_LAMINAR] = laminarPred.twVar;
    muTw[REGIME_FREESTREAM] = freestreamPred.twMean; sigma2Tw[REGIME_FREESTREAM] = freestreamPred.twVar;
    muTw[REGIME_CHANNEL] = channelPred.twMean;       sigma2Tw[REGIME_CHANNEL] = channelPred.twVar;
    muTw[REGIME_COUETTE] = couettePred.twMean;       sigma2Tw[REGIME_COUETTE] = couettePred.twVar;
    muTw[REGIME_CROSSFLOW] = crossflowPred.twMean;   sigma2Tw[REGIME_CROSSFLOW] = crossflowPred.twVar;

    // Fold in (C,Z)-calibration-driven epistemic uncertainty (see
    // unified_wall_model.H's own comment on tau2Channel_/tau2Couette_) -
    // ADDED onto each regime's own NN-predicted variance, BEFORE combine()
    // runs, so the existing, already-verified combine() formula needs no
    // structural change at all: the (Var(W_i)+w_i^2) weighting it already
    // applies to sigma2Tw[i] is exactly the mathematically correct weight
    // for this new term too (re-derived the same way as the original
    // formula - see Uncertainty_propagation/README.md). tw ONLY - angle
    // uncertainty (sigma2Angle) is deliberately untouched, matching the
    // scope this analysis was requested for.
    sigma2Tw[REGIME_FREESTREAM] += tau2Channel_;
    sigma2Tw[REGIME_CHANNEL] += tau2Channel_;
    sigma2Tw[REGIME_CROSSFLOW] += tau2Channel_;
    sigma2Tw[REGIME_COUETTE] += tau2Couette_;
    // REGIME_LAMINAR: no addition - closed-form analytical solution, no
    // ODT/(C,Z) dependence.

    muAngle[REGIME_LAMINAR] = laminarPred.angleMean;      sigma2Angle[REGIME_LAMINAR] = laminarPred.angleVar;
    muAngle[REGIME_FREESTREAM] = freestreamPred.angleMean; sigma2Angle[REGIME_FREESTREAM] = freestreamPred.angleVar;
    muAngle[REGIME_CHANNEL] = channelPred.angleMean;       sigma2Angle[REGIME_CHANNEL] = channelPred.angleVar;
    muAngle[REGIME_COUETTE] = couettePred.angleMean;       sigma2Angle[REGIME_COUETTE] = couettePred.angleVar;
    muAngle[REGIME_CROSSFLOW] = crossflowPred.angleMean;   sigma2Angle[REGIME_CROSSFLOW] = crossflowPred.angleVar;

    Prediction result;
    result.weightsMean = weightsMean;
    result.weightsCov = weightsCov;
    for (int i = 0; i < WM_N_REGIMES; ++i)
    {
        result.weightsSigma[i] = std::sqrt(weightsCov[i][i]);
    }

    scalar twVar, angleVar;
    std::tie(result.twMean, twVar) = combine(weightsMean, weightsCov, muTw, sigma2Tw);
    std::tie(result.angleMean, angleVar) = combineCircular(weightsMean, muAngle, sigma2Angle);
    result.twSigma = std::sqrt(twVar);
    result.angleSigma = std::sqrt(angleVar);

    return result;
}

std::pair<scalar, scalar> classifierWallModel::combineCircular
(
    const WeightVector& weightsMean,
    const std::array<scalar, WM_N_REGIMES>& mu,
    const std::array<scalar, WM_N_REGIMES>& sigma2
)
{
    scalar Zx = 0.0, Zy = 0.0;
    for (int i = 0; i < WM_N_REGIMES; ++i)
    {
        const scalar sigma_i = std::sqrt(sigma2[i]);
        const scalar R_i = std::exp(-sigma_i * sigma_i / 2.0);
        Zx += weightsMean[i] * R_i * std::cos(mu[i]);
        Zy += weightsMean[i] * R_i * std::sin(mu[i]);
    }

    const scalar combinedMean = std::atan2(Zy, Zx);
    const scalar RBar = std::min(std::sqrt(Zx*Zx + Zy*Zy), scalar(1.0 - 1e-12));
    const scalar combinedVar = -2.0 * std::log(RBar);

    return {combinedMean, combinedVar};
}

}
