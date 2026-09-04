#ifndef regime_wall_model_H
#define regime_wall_model_H

#include <array>
#include <memory>
#include <utility>
#include <vector>
#include "scalar.H"
#include "fileName.H"
#include "mlForwardPass_full.H"

namespace Foam
{

constexpr int WM_N_ENSEMBLE_MEMBERS = 5;

// One target's (tw, or the "sin"/"cos" half of a circular tw_direction_angle)
// 5-member ensemble - unchanged wallModelMLP members, moment-matching
// combination, and a target-specific output back-transform. Identical class
// (by design) to every individual building block's own wallModelEnsemble.
class wallModelEnsemble
{
private:

    std::array<wallModelMLP, WM_N_ENSEMBLE_MEMBERS> members_;
    scalar yMean_;
    scalar yScale_;
    bool applyOutputArcsinh_;

public:

    explicit wallModelEnsemble(const fileName& dir);

    std::pair<scalar, scalar> predict(const std::vector<scalar>& x) const;
};

// Wraps EITHER a plain wallModelEnsemble (the original, direct-angle-
// regression representation - still correct and used for laminar/channel/
// crossflow, which have no known wraparound issue) OR a pair of independent
// sin(angle)/cos(angle) wallModelEnsemble instances combined via atan2 (the
// fix couette's own tw_direction_angle needed - see optimizer-couette/
// couette neural network/retrain_tw_direction_angle_circular.py's module
// docstring for the full derivation of why direct-angle regression breaks
// down at couette's own extreme, near-+-pi pressure-differential conditions,
// and this formula:
//     angle_mean = atan2(mean_sin, mean_cos)
//     Var(angle) = (mean_cos^2*var_sin + mean_sin^2*var_cos)
//                  / (mean_sin^2+mean_cos^2)^2
// - the standard delta-method for atan2 of two independent Gaussian
// estimates, verified against Monte Carlo before being trusted in the Python
// training pipeline this mirrors).
//
// Which representation a given directory uses is read from an explicit
// "angle_type.txt" marker file (contents: "plain" or "circular") written by
// the corresponding export script - not inferred from which subdirectories
// happen to exist, to avoid silently guessing wrong. If angle_type.txt is
// ABSENT (every wmLayers_openfoam/ export produced before this couette fix
// was integrated here), defaults to "plain" - preserves old deployments'
// exact existing behavior with no action needed.
class angleModel
{
private:

    bool circular_;
    std::unique_ptr<wallModelEnsemble> plainEnsemble_;
    std::unique_ptr<wallModelEnsemble> sinEnsemble_;
    std::unique_ptr<wallModelEnsemble> cosEnsemble_;

public:

    explicit angleModel(const fileName& dir);

    // Returns (angleMean, angleSigma) - same external interface regardless
    // of which internal representation this instance actually uses.
    std::pair<scalar, scalar> predict(const std::vector<scalar>& x) const;
};

// One building block's full (tw, tw_direction_angle) model - a
// wallModelEnsemble for tw plus an angleModel (plain or circular, auto-
// detected) for tw_direction_angle, plus the shared log1p/arcsinh input
// transform. One of these is instantiated per DISTINCT trained model
// (laminar, channel, couette, crossflow - NOT one per classifier label,
// since "freestream" has no separately-trained model and aliases channel's
// - see unified_wall_model.H).
class regimeWallModel
{
private:

    std::vector<int> log1pCols_;
    std::vector<int> arcsinhCols_;
    wallModelEnsemble twEnsemble_;
    angleModel angleModel_;

public:

    struct Prediction
    {
        scalar twMean;
        scalar twVar;
        scalar angleMean;
        scalar angleVar;
    };

    // dir should contain tw/ (a wallModelEnsemble directory) and
    // tw_direction_angle/ (either a plain wallModelEnsemble directory, or a
    // circular directory containing sin/ and cos/ subfolders plus
    // angle_type.txt) plus input_metadata.txt (log1p_cols, arcsinh_cols) -
    // the same wmLayers_openfoam/<regime>/ layout every individual building
    // block's own OpenFOAM integration already produces.
    explicit regimeWallModel(const fileName& dir);

    // rawInputs: (Re_center, Re_top, ax_n, gamma12, lambda_roc, tke_n), raw
    // physical values. Returns variances (not sigmas) directly, since
    // unifiedWallModel needs variances for the blend combination formula
    // and would otherwise have to immediately square the sigma back.
    Prediction predict(const std::vector<scalar>& rawInputs) const;
};

}

#endif
