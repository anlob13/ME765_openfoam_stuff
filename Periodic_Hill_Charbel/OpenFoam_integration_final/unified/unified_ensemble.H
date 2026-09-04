#ifndef unified_ensemble_H
#define unified_ensemble_H

#include <array>
#include <utility>
#include <vector>
#include "scalar.H"
#include "fileName.H"
#include "mlForwardPass_full.H"

namespace Foam
{

constexpr int WM_N_MEMBERS = 5;
constexpr int WM_N_REGIMES = 5;   // laminar, freestream, channel, couette, crossflow
// Raw per-member output width: see ../neural_network_v4/export_unified_for_
// openfoam.py's module docstring for the fixed slot order this depends on.
constexpr int WM_N_RAW_OUTPUTS = WM_N_REGIMES + 2 + 2 + 2;  // cls(5) + tw(2) + sin(2) + cos(2)

using WeightVector = std::array<scalar, WM_N_REGIMES>;
using WeightCovMatrix = std::array<std::array<scalar, WM_N_REGIMES>, WM_N_REGIMES>;

// unifiedEnsemble replaces classifier_ensemble.C + regime_wall_model.C's
// wallModelEnsemble/angleModel entirely: instead of 5 independently-trained
// classifier networks blending 4 SEPARATELY-trained per-regime regression
// ensembles at inference time, this is ONE shared-trunk, 4-head network (see
// ../neural_network_v4/train_unified_v4.py's module docstring for the full
// design rationale - classification is now a genuinely auxiliary output that
// does not gate or blend the regression prediction at all).
//
// Each of the WM_N_MEMBERS members is an unmodified wallModelMLP producing
// WM_N_RAW_OUTPUTS raw values in the fixed order documented in
// export_unified_for_openfoam.py:
//     [0:5]  raw classification evidence logits (softplus+1 -> alpha here, not
//            inside the shared MLP runner - same established pattern as the
//            original classifier_ensemble.C)
//     [5:7]  tw (mu, log_var) on arcsinh(tw)
//     [7:9]  sin(tw_direction_angle) (mu, log_var)
//     [9:11] cos(tw_direction_angle) (mu, log_var)
//
// Combination formulas (identical to the ones this replaces, just computed
// from ONE ensemble's raw outputs instead of two SEPARATE ensembles' outputs):
//   - classification: per-member Dirichlet moments + aleatoric/epistemic
//     covariance combine - byte-for-byte the same math as the old
//     classifier_ensemble.C (see that file's own header comment, preserved
//     here for reference).
//   - tw / sin / cos: per-target deep-ensemble moment matching (aleatoric =
//     mean of member variances, epistemic = variance of member means),
//     exactly regime_wall_model.C's wallModelEnsemble - just run once instead
//     of once per regime.
//   - angle: atan2(mean_sin, mean_cos) + delta-method variance propagation,
//     exactly regime_wall_model.C's angleModel (circular branch) - always
//     circular here (no "plain" variant), since the whole point of unifying
//     across regimes needs one angle representation that works everywhere,
//     including couette's near-+-pi wraparound.
class unifiedEnsemble
{
private:

    std::array<wallModelMLP, WM_N_MEMBERS> members_;
    std::vector<int> log1pCols_;
    std::vector<int> arcsinhCols_;

    scalar twMean_, twScale_;
    scalar sinMean_, sinScale_;
    scalar cosMean_, cosScale_;

public:

    struct Prediction
    {
        WeightVector weightsMean;
        WeightCovMatrix weightsCov;
        scalar twMean;
        scalar twVar;
        scalar angleMean;
        scalar angleVar;
    };

    // dir should contain member_0/ .. member_4/ (each loadable by wallModelMLP
    // unchanged), input_metadata.txt (log1p_cols, arcsinh_cols), and
    // output_metadata.txt (tw_mean/scale, sin_mean/scale, cos_mean/scale) -
    // written by ../neural_network_v4/export_unified_for_openfoam.py.
    explicit unifiedEnsemble(const fileName& dir);

    // rawInputs: (Re_center, Re_top, ax_n, gamma12, lambda_roc, tke_n), raw
    // physical values - the log1p/arcsinh input transform happens internally.
    Prediction predict(const std::vector<scalar>& rawInputs) const;
};

}

#endif
