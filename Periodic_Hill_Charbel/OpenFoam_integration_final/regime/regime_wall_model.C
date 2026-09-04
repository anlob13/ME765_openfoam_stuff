// wallModelEnsemble / angleModel / regimeWallModel - identical logic to every
// individual building block's own wall_model_ensemble_combine.C for the tw
// path (moment-matching combination, sinh/cosh vs. linear output back-
// transform decided per-target from ensemble_metadata.txt, not hardcoded).
// angleModel is new here (see regime_wall_model.H's own docstring) - adds
// support for couette's circular sin/cos tw_direction_angle representation
// alongside the original plain direct-angle regression every other regime
// still uses, auto-selected from an explicit angle_type.txt marker.

#include <cmath>
#include <fstream>
#include <sstream>
#include <map>
#include <string>
#include <tuple>
#include <utility>

#include "regime_wall_model.H"

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

    bool fileExists(const fileName& path)
    {
        std::ifstream f(path);
        return static_cast<bool>(f);
    }

    // Reads angle_type.txt's single word ("plain" or "circular"). Missing
    // file -> "plain", matching every pre-existing wmLayers_openfoam/ export
    // (see regime_wall_model.H's own docstring for why this default is safe).
    std::string readAngleType(const fileName& dir)
    {
        const fileName path = dir/"angle_type.txt";
        if (!fileExists(path))
        {
            return "plain";
        }
        std::ifstream f(path);
        std::string word;
        f >> word;
        if (word != "plain" && word != "circular")
        {
            FatalErrorInFunction << "angle_type.txt at " << path
                << " must contain 'plain' or 'circular', got '" << word << "'"
                << exit(FatalError);
        }
        return word;
    }
}

// ---- wallModelEnsemble ----

wallModelEnsemble::wallModelEnsemble(const fileName& dir)
:
    members_{{
        wallModelMLP(dir/"member_0"),
        wallModelMLP(dir/"member_1"),
        wallModelMLP(dir/"member_2"),
        wallModelMLP(dir/"member_3"),
        wallModelMLP(dir/"member_4")
    }},
    yMean_(0),
    yScale_(1),
    applyOutputArcsinh_(true)
{
    const auto meta = readKeyValueFile(dir/"ensemble_metadata.txt");
    yMean_ = meta.at("y_mean").at(0);
    yScale_ = meta.at("y_scale").at(0);
    applyOutputArcsinh_ = meta.at("apply_target_arcsinh").at(0) != 0.0;
}

std::pair<scalar, scalar> wallModelEnsemble::predict
(
    const std::vector<scalar>& x
) const
{
    std::array<scalar, WM_N_ENSEMBLE_MEMBERS> mus;
    std::array<scalar, WM_N_ENSEMBLE_MEMBERS> vars;

    for (int k = 0; k < WM_N_ENSEMBLE_MEMBERS; ++k)
    {
        std::vector<scalar> out = members_[k].predict(x);
        mus[k] = out[0];
        vars[k] = std::exp(out[1]);
    }

    scalar muZ = 0.0;
    scalar varMean = 0.0;
    for (int k = 0; k < WM_N_ENSEMBLE_MEMBERS; ++k) { muZ += mus[k]; varMean += vars[k]; }
    muZ /= WM_N_ENSEMBLE_MEMBERS;
    varMean /= WM_N_ENSEMBLE_MEMBERS;

    scalar muDisagreement = 0.0;
    for (int k = 0; k < WM_N_ENSEMBLE_MEMBERS; ++k)
    {
        const scalar d = mus[k] - muZ;
        muDisagreement += d*d;
    }
    muDisagreement /= WM_N_ENSEMBLE_MEMBERS;

    const scalar varZ = varMean + muDisagreement;

    const scalar muLin = muZ*yScale_ + yMean_;
    const scalar varLin = varZ*yScale_*yScale_;

    if (applyOutputArcsinh_)
    {
        const scalar prediction = std::sinh(muLin);
        const scalar sigma = std::cosh(muLin)*std::sqrt(varLin);
        return {prediction, sigma};
    }
    else
    {
        return {muLin, std::sqrt(varLin)};
    }
}


// ---- angleModel ----

angleModel::angleModel(const fileName& dir)
:
    circular_(readAngleType(dir) == "circular"),
    plainEnsemble_(),
    sinEnsemble_(),
    cosEnsemble_()
{
    if (circular_)
    {
        sinEnsemble_.reset(new wallModelEnsemble(dir/"sin"));
        cosEnsemble_.reset(new wallModelEnsemble(dir/"cos"));
    }
    else
    {
        plainEnsemble_.reset(new wallModelEnsemble(dir));
    }
}

std::pair<scalar, scalar> angleModel::predict(const std::vector<scalar>& x) const
{
    if (!circular_)
    {
        return plainEnsemble_->predict(x);
    }

    scalar meanSin, sigmaSin, meanCos, sigmaCos;
    std::tie(meanSin, sigmaSin) = sinEnsemble_->predict(x);
    std::tie(meanCos, sigmaCos) = cosEnsemble_->predict(x);

    const scalar varSin = sigmaSin*sigmaSin;
    const scalar varCos = sigmaCos*sigmaCos;

    const scalar angleMean = std::atan2(meanSin, meanCos);
    const scalar r2 = std::max(meanSin*meanSin + meanCos*meanCos, scalar(1e-12));
    const scalar angleVar = (meanCos*meanCos*varSin + meanSin*meanSin*varCos) / (r2*r2);

    return {angleMean, std::sqrt(angleVar)};
}


// ---- regimeWallModel ----

regimeWallModel::regimeWallModel(const fileName& dir)
:
    log1pCols_(),
    arcsinhCols_(),
    twEnsemble_(dir/"tw"),
    angleModel_(dir/"tw_direction_angle")
{
    const auto meta = readKeyValueFile(dir/"input_metadata.txt");
    log1pCols_ = toIntVector(meta.at("log1p_cols"));
    arcsinhCols_ = toIntVector(meta.at("arcsinh_cols"));
}

regimeWallModel::Prediction regimeWallModel::predict
(
    const std::vector<scalar>& rawInputs
) const
{
    std::vector<scalar> x = rawInputs;

    for (int c : log1pCols_)   x[c] = std::log1p(x[c]);
    for (int c : arcsinhCols_) x[c] = std::asinh(x[c]);

    Prediction result;
    scalar twSigma, angleSigma;
    std::tie(result.twMean, twSigma) = twEnsemble_.predict(x);
    std::tie(result.angleMean, angleSigma) = angleModel_.predict(x);
    result.twVar = twSigma * twSigma;
    result.angleVar = angleSigma * angleSigma;
    return result;
}

}
