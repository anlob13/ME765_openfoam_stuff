#include <cmath>
#include <fstream>
#include <sstream>
#include <map>
#include <string>

#include "unified_ensemble.H"

namespace Foam
{

namespace
{
    std::map<std::string, std::vector<double>> readKeyValueFile(const fileName& path)
    {
        std::ifstream f(path);
        if (!f)
        {
            FatalErrorInFunction << "Cannot open " << path << exit(FatalError);
        }

        std::map<std::string, std::vector<double>> result;
        std::string line;
        while (std::getline(f, line))
        {
            std::istringstream iss(line);
            std::string key;
            iss >> key;
            if (key.empty()) continue;

            std::vector<double> values;
            double v;
            while (iss >> v) values.push_back(v);
            result[key] = values;
        }
        return result;
    }

    std::vector<int> toIntVector(const std::vector<double>& v)
    {
        std::vector<int> result;
        result.reserve(v.size());
        for (double x : v) result.push_back(static_cast<int>(x));
        return result;
    }

    inline scalar softplus(scalar x)
    {
        return std::max<scalar>(x, 0) + std::log1p(std::exp(-std::fabs(x)));
    }
}

unifiedEnsemble::unifiedEnsemble(const fileName& dir)
:
    members_{{
        wallModelMLP(dir/"member_0"),
        wallModelMLP(dir/"member_1"),
        wallModelMLP(dir/"member_2"),
        wallModelMLP(dir/"member_3"),
        wallModelMLP(dir/"member_4")
    }},
    log1pCols_(),
    arcsinhCols_(),
    twMean_(0), twScale_(1),
    sinMean_(0), sinScale_(1),
    cosMean_(0), cosScale_(1)
{
    const auto inMeta = readKeyValueFile(dir/"input_metadata.txt");
    log1pCols_ = toIntVector(inMeta.at("log1p_cols"));
    arcsinhCols_ = toIntVector(inMeta.at("arcsinh_cols"));

    const auto outMeta = readKeyValueFile(dir/"output_metadata.txt");
    twMean_ = outMeta.at("tw_mean").at(0);
    twScale_ = outMeta.at("tw_scale").at(0);
    sinMean_ = outMeta.at("sin_mean").at(0);
    sinScale_ = outMeta.at("sin_scale").at(0);
    cosMean_ = outMeta.at("cos_mean").at(0);
    cosScale_ = outMeta.at("cos_scale").at(0);
}

unifiedEnsemble::Prediction unifiedEnsemble::predict
(
    const std::vector<scalar>& rawInputs
) const
{
    std::vector<scalar> x = rawInputs;
    for (int c : log1pCols_)   x[c] = std::log1p(x[c]);
    for (int c : arcsinhCols_) x[c] = std::asinh(x[c]);

    std::array<WeightVector, WM_N_MEMBERS> memberAlphaMeans;
    std::array<WeightCovMatrix, WM_N_MEMBERS> memberAlphaCovs;
    std::array<scalar, WM_N_MEMBERS> twMus, twVars;
    std::array<scalar, WM_N_MEMBERS> sinMus, sinVars;
    std::array<scalar, WM_N_MEMBERS> cosMus, cosVars;

    for (int k = 0; k < WM_N_MEMBERS; ++k)
    {
        // Each member: Normalization -> Dense(tanh) x4 -> Dense(linear, 11
        // outputs), exactly as loaded from member_k/layer_metadata.txt. See
        // unified_ensemble.H for the fixed output slot order.
        std::vector<scalar> raw = members_[k].predict(x);

        WeightVector alpha;
        scalar S = 0.0;
        for (int i = 0; i < WM_N_REGIMES; ++i)
        {
            alpha[i] = softplus(raw[i]) + 1.0;
            S += alpha[i];
        }
        for (int i = 0; i < WM_N_REGIMES; ++i)
        {
            memberAlphaMeans[k][i] = alpha[i] / S;
        }
        for (int i = 0; i < WM_N_REGIMES; ++i)
        {
            for (int j = 0; j < WM_N_REGIMES; ++j)
            {
                const scalar diagTerm = (i == j) ? S : 0.0;
                memberAlphaCovs[k][i][j] = alpha[i] * (diagTerm - alpha[j]) / (S * S * (S + 1.0));
            }
        }

        twMus[k] = raw[5];             twVars[k] = std::exp(raw[6]);
        sinMus[k] = raw[7];             sinVars[k] = std::exp(raw[8]);
        cosMus[k] = raw[9];             cosVars[k] = std::exp(raw[10]);
    }

    // ---- classification: Dirichlet ensemble combine (identical math to the
    // original classifier_ensemble.C) ----
    WeightVector combinedMean{};
    for (int i = 0; i < WM_N_REGIMES; ++i)
    {
        scalar sum = 0.0;
        for (int k = 0; k < WM_N_MEMBERS; ++k) sum += memberAlphaMeans[k][i];
        combinedMean[i] = sum / WM_N_MEMBERS;
    }
    WeightCovMatrix combinedCov{};
    for (int i = 0; i < WM_N_REGIMES; ++i)
    {
        for (int j = 0; j < WM_N_REGIMES; ++j)
        {
            scalar aleatoric = 0.0;
            for (int k = 0; k < WM_N_MEMBERS; ++k) aleatoric += memberAlphaCovs[k][i][j];
            aleatoric /= WM_N_MEMBERS;

            scalar epistemic = 0.0;
            for (int k = 0; k < WM_N_MEMBERS; ++k)
            {
                epistemic += (memberAlphaMeans[k][i] - combinedMean[i]) * (memberAlphaMeans[k][j] - combinedMean[j]);
            }
            epistemic /= WM_N_MEMBERS;

            combinedCov[i][j] = aleatoric + epistemic;
        }
    }

    // ---- tw: deep-ensemble moment matching + sinh back-transform (identical
    // math to the original regime_wall_model.C's wallModelEnsemble) ----
    scalar twMuZ = 0.0, twVarMean = 0.0;
    for (int k = 0; k < WM_N_MEMBERS; ++k) { twMuZ += twMus[k]; twVarMean += twVars[k]; }
    twMuZ /= WM_N_MEMBERS; twVarMean /= WM_N_MEMBERS;
    scalar twDisagreement = 0.0;
    for (int k = 0; k < WM_N_MEMBERS; ++k) { const scalar d = twMus[k] - twMuZ; twDisagreement += d*d; }
    twDisagreement /= WM_N_MEMBERS;
    const scalar twVarZ = twVarMean + twDisagreement;
    const scalar twMuLin = twMuZ * twScale_ + twMean_;
    const scalar twVarLin = twVarZ * twScale_ * twScale_;
    const scalar twMeanOut = std::sinh(twMuLin);
    const scalar twSigmaOut = std::cosh(twMuLin) * std::sqrt(twVarLin);

    // ---- angle: independent sin/cos deep-ensemble moment matching, then
    // atan2 combine (identical math to the original regime_wall_model.C's
    // angleModel circular branch) ----
    auto ensembleCombine1D = [](const std::array<scalar, WM_N_MEMBERS>& mus,
                                 const std::array<scalar, WM_N_MEMBERS>& vars,
                                 scalar yMean, scalar yScale, scalar& outMean, scalar& outVar)
    {
        scalar muZ = 0.0, varMean = 0.0;
        for (int k = 0; k < WM_N_MEMBERS; ++k) { muZ += mus[k]; varMean += vars[k]; }
        muZ /= WM_N_MEMBERS; varMean /= WM_N_MEMBERS;
        scalar disagreement = 0.0;
        for (int k = 0; k < WM_N_MEMBERS; ++k) { const scalar d = mus[k] - muZ; disagreement += d*d; }
        disagreement /= WM_N_MEMBERS;
        const scalar varZ = varMean + disagreement;
        outMean = muZ * yScale + yMean;
        outVar = varZ * yScale * yScale;
    };

    scalar meanSin, varSin, meanCos, varCos;
    ensembleCombine1D(sinMus, sinVars, sinMean_, sinScale_, meanSin, varSin);
    ensembleCombine1D(cosMus, cosVars, cosMean_, cosScale_, meanCos, varCos);

    const scalar angleMeanOut = std::atan2(meanSin, meanCos);
    const scalar r2 = std::max(meanSin*meanSin + meanCos*meanCos, scalar(1e-12));
    const scalar angleVarOut = (meanCos*meanCos*varSin + meanSin*meanSin*varCos) / (r2*r2);

    Prediction result;
    result.weightsMean = combinedMean;
    result.weightsCov = combinedCov;
    result.twMean = twMeanOut;
    result.twVar = twSigmaOut * twSigmaOut;
    result.angleMean = angleMeanOut;
    result.angleVar = angleVarOut;
    return result;
}

}
