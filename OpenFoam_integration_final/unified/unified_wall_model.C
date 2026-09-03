#include <cmath>
#include <fstream>
#include <sstream>
#include <tuple>
#include <utility>

#include "unified_wall_model.H"

namespace Foam
{

namespace
{
    // Optional file, two "key value" lines (tau2_channel, tau2_couette) - see
    // unified_wall_model.H's own comment on tau2Channel_/tau2Couette_ for what
    // these mean and why missing = (0,0) is the correct, backward-compatible
    // default. Unchanged from the original unified_wall_model.C.
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

unifiedWallModel::unifiedWallModel(const fileName& wmLayersDir)
:
    ensemble_(wmLayersDir/"unified"),
    tau2Channel_(0.0),
    tau2Couette_(0.0)
{
    std::tie(tau2Channel_, tau2Couette_) = readCZUncertainty(wmLayersDir);
}

unifiedWallModel::Prediction unifiedWallModel::predict
(
    const std::vector<scalar>& rawInputs
) const
{
    const unifiedEnsemble::Prediction pred = ensemble_.predict(rawInputs);

    // Fold in (C,Z)-calibration-driven epistemic uncertainty (see
    // unified_wall_model.H's own comment on tau2Channel_/tau2Couette_) -
    // blended by the classification weights, since there is only one tw
    // output left to add it to (unlike the old per-regime design, which added
    // it to each regime's own sigma2 before its own separate blend).
    const scalar tau2Effective =
        pred.weightsMean[REGIME_FREESTREAM] * tau2Channel_
      + pred.weightsMean[REGIME_CHANNEL] * tau2Channel_
      + pred.weightsMean[REGIME_CROSSFLOW] * tau2Channel_
      + pred.weightsMean[REGIME_COUETTE] * tau2Couette_;
      // REGIME_LAMINAR: no addition - closed-form analytical solution, no
      // ODT/(C,Z) dependence.

    Prediction result;
    result.weightsMean = pred.weightsMean;
    result.weightsCov = pred.weightsCov;
    for (int i = 0; i < WM_N_REGIMES; ++i)
    {
        result.weightsSigma[i] = std::sqrt(pred.weightsCov[i][i]);
    }

    result.twMean = pred.twMean;
    result.twSigma = std::sqrt(pred.twVar + tau2Effective);
    result.angleMean = pred.angleMean;
    result.angleSigma = std::sqrt(pred.angleVar);

    return result;
}

}
