#include <cmath>
#include <fstream>
#include <sstream>
#include <map>
#include <string>

#include "classifier_ensemble.H"

namespace Foam
{

namespace
{
    // Same minimal "key value value..." reader every building block's own
    // wall_model_ensemble_combine.C already uses for ensemble_metadata.txt/
    // input_metadata.txt - avoids hardcoding column indices as C++ constants.
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
        // Numerically stable softplus: log(1+exp(x)) overflows for large x,
        // so use the standard identity softplus(x) = max(x,0) + log1p(exp(-|x|)).
        return std::max<scalar>(x, 0) + std::log1p(std::exp(-std::fabs(x)));
    }
}

classifierEnsemble::classifierEnsemble(const fileName& dir)
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
    sqrtCols_()
{
    const auto meta = readKeyValueFile(dir/"input_metadata.txt");
    log1pCols_ = toIntVector(meta.at("log1p_cols"));
    arcsinhCols_ = toIntVector(meta.at("arcsinh_cols"));
    // Optional - see the .H file comment for why a missing key is not an error.
    if (meta.count("sqrt_cols"))
    {
        sqrtCols_ = toIntVector(meta.at("sqrt_cols"));
    }
}

std::pair<WeightVector, WeightCovMatrix> classifierEnsemble::predict
(
    const std::vector<scalar>& rawInputs
) const
{
    std::vector<scalar> x = rawInputs;
    for (int c : log1pCols_)   x[c] = std::log1p(x[c]);
    for (int c : arcsinhCols_) x[c] = std::asinh(x[c]);
    for (int c : sqrtCols_)    x[c] = std::sqrt(std::max<scalar>(x[c], 0.0));

    std::array<WeightVector, WM_N_MEMBERS> memberMeans;
    std::array<WeightCovMatrix, WM_N_MEMBERS> memberCovs;

    for (int k = 0; k < WM_N_MEMBERS; ++k)
    {
        // Each member's predict() runs Normalization -> Dense(tanh) x4 ->
        // Dense(linear, 5 outputs) exactly as loaded from
        // member_k/layer_metadata.txt - the RAW evidence logits (see the .H
        // file comment for why softplus isn't baked into the exported layers).
        std::vector<scalar> rawEvidence = members_[k].predict(x);

        WeightVector alpha;
        scalar S = 0.0;
        for (int i = 0; i < WM_N_REGIMES; ++i)
        {
            alpha[i] = softplus(rawEvidence[i]) + 1.0;
            S += alpha[i];
        }

        for (int i = 0; i < WM_N_REGIMES; ++i)
        {
            memberMeans[k][i] = alpha[i] / S;
        }

        for (int i = 0; i < WM_N_REGIMES; ++i)
        {
            for (int j = 0; j < WM_N_REGIMES; ++j)
            {
                const scalar diagTerm = (i == j) ? S : 0.0;
                memberCovs[k][i][j] = alpha[i] * (diagTerm - alpha[j]) / (S * S * (S + 1.0));
            }
        }
    }

    WeightVector combinedMean{};
    for (int i = 0; i < WM_N_REGIMES; ++i)
    {
        scalar sum = 0.0;
        for (int k = 0; k < WM_N_MEMBERS; ++k) sum += memberMeans[k][i];
        combinedMean[i] = sum / WM_N_MEMBERS;
    }

    WeightCovMatrix combinedCov{};
    for (int i = 0; i < WM_N_REGIMES; ++i)
    {
        for (int j = 0; j < WM_N_REGIMES; ++j)
        {
            // aleatoric: mean of the members' own covariances
            scalar aleatoric = 0.0;
            for (int k = 0; k < WM_N_MEMBERS; ++k) aleatoric += memberCovs[k][i][j];
            aleatoric /= WM_N_MEMBERS;

            // epistemic: population covariance (divide by M, not M-1, matching
            // np.var()'s default ddof=0 used when this formula was validated)
            // of the members' own mean vectors
            scalar epistemic = 0.0;
            for (int k = 0; k < WM_N_MEMBERS; ++k)
            {
                epistemic += (memberMeans[k][i] - combinedMean[i]) * (memberMeans[k][j] - combinedMean[j]);
            }
            epistemic /= WM_N_MEMBERS;

            combinedCov[i][j] = aleatoric + epistemic;
        }
    }

    return {combinedMean, combinedCov};
}

}
