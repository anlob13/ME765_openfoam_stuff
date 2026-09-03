#ifndef classifier_ensemble_H
#define classifier_ensemble_H

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

using WeightVector = std::array<scalar, WM_N_REGIMES>;
using WeightCovMatrix = std::array<std::array<scalar, WM_N_REGIMES>, WM_N_REGIMES>;

// The building-block classifier: 5 independently-trained evidential (Dirichlet)
// networks, each an unmodified wallModelMLP run through to a RAW linear output
// (mlForwardPass_full.C only implements tanh/relu/linear - no softplus - so the
// softplus(x)+1 evidence->alpha step happens here, not inside the shared MLP
// runner, exactly like every regression building block already keeps its own
// non-Layer-representable transforms in its own wrapper .C file instead of
// touching mlForwardPass_full.C).
//
// Per-member Dirichlet moments (alpha = softplus(rawEvidence)+1, S = sum(alpha)):
//     mean_k       = alpha_k / S                         (sums to exactly 1)
//     Cov(k,l)     = (S*delta_kl - alpha_l)*alpha_k / (S^2*(S+1))
//                    (diagonal reduces to mean_k*(1-mean_k)/(S+1))
//
// Combined across the 5 members via this project's usual aleatoric+epistemic
// moment-matching (Lakshminarayanan, Pritzel & Blundell 2017), generalized from
// a scalar variance to a full 5x5 covariance matrix (weights are correlated
// with each other by the simplex constraint, so a plain per-class variance
// would throw away real information a full covariance keeps):
//     combined mean       = mean_k(member_k mean)
//     combined covariance = mean_k(member_k covariance)      <- aleatoric
//                          + Cov_k(member_k mean vectors)     <- epistemic
// Both the per-member Dirichlet formulas and this combination were checked
// against a direct Monte Carlo simulation before being trusted for anything
// real - see classifier_building_block/smoke_test_classifier.py and
// verify_combination.py in this project (same formulas, Python side).
class classifierEnsemble
{
private:

    std::array<wallModelMLP, WM_N_MEMBERS> members_;
    std::vector<int> log1pCols_;
    std::vector<int> arcsinhCols_;
    // sqrt(clamp(x,0)) transform, added for this classifier's own tke_n
    // input (v1_improved's own change - verified via excess kurtosis to fit
    // better than log1p for this feature; see classifier_building_block/
    // train_classifier_v1_improved.py's own module docstring). Empty vector
    // (no-op) if input_metadata.txt has no "sqrt_cols" line at all - keeps
    // this class working unchanged against the ORIGINAL v1 classifier's
    // 2-transform-category exports too, not just this 3-category one.
    std::vector<int> sqrtCols_;

public:

    // dir should contain member_0/ .. member_4/ (each loadable by wallModelMLP
    // unchanged) and input_metadata.txt (log1p_cols, arcsinh_cols, and
    // optionally sqrt_cols) - written by export_classifier_for_openfoam.py.
    explicit classifierEnsemble(const fileName& dir);

    // rawInputs: (Re_center, Re_top, ax_n, gamma12, lambda_roc, tke_n), raw
    // physical values - the log1p/arcsinh input transform happens internally,
    // same convention as every regression building block.
    // Returns (weightsMean, weightsCov) - weightsMean sums to 1 to floating
    // point precision, in label order [laminar, freestream, channel, couette,
    // crossflow].
    std::pair<WeightVector, WeightCovMatrix> predict(const std::vector<scalar>& rawInputs) const;
};

}

#endif
