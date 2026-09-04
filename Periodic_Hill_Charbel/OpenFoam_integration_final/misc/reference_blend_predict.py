"""Python reference for THIS deployment - same blend formula (linear for tw,
circular-resultant-vector for tw_direction_angle) as classifier_building_
block/blend_predict.py and '../V1 classifier updated Re selective data/
reference_blend_predict.py', copied here (not imported) because MODEL_PATHS
below points at THIS package's own sources - not identical to that prior
version's paths (couette and crossflow both changed, see comments below).
COMBINATION LOGIC (BuildingBlockModel.predict, _combine, _combine_circular)
is reused verbatim, unchanged.

Used by generate_verification_rows.py as the Python side of the C++
verification (test_main.C).
"""
import os
import torch

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(THIS_DIR)

LABEL_NAMES = ["laminar", "freestream", "channel", "couette", "crossflow"]
MODEL_ALIAS = {"freestream": "channel"}

# laminar/channel: unchanged from every prior integration in this project.
# couette: CHANGED - now the backflow-filtered retrain ('../updated couette
#   no backflow building block/'), not 'couette building block updated Re/'.
# crossflow: CHANGED - now the full-scale (all 4 Re_tau) Aug-11 retrain via
#   its canonical location ('lukas_cross_flow_download/cross_flow_building_
#   block/'), not the older reduced-scale (2-of-4 Re_tau) 'crossflow
#   building block updated Re/' - see export_regimes_for_openfoam.py's own
#   module docstring for why this is a deliberate correction, not an
#   arbitrary choice.
MODEL_PATHS = {
    "laminar": (
        os.path.join(PROJECT_ROOT, "wall_model_ai", "laminar_building_block", "tw_wall_model_deployable.pt"),
        os.path.join(PROJECT_ROOT, "wall_model_ai", "laminar_building_block", "tw_direction_angle_wall_model_deployable.pt"),
    ),
    "channel": (
        os.path.join(PROJECT_ROOT, "channel building block updated Re", "test_NEW_RE_tw_wall_model_deployable.pt"),
        os.path.join(PROJECT_ROOT, "channel building block updated Re", "test_NEW_RE_tw_direction_angle_wall_model_deployable.pt"),
    ),
    "couette": (
        os.path.join(PROJECT_ROOT, "updated couette no backflow building block", "tw_wall_model_deployable.pt"),
        os.path.join(PROJECT_ROOT, "updated couette no backflow building block", "tw_direction_angle_wall_model_deployable.pt"),
    ),
    "crossflow": (
        os.path.join(PROJECT_ROOT, "lukas_cross_flow_download", "cross_flow_building_block", "tw_wall_model_deployable.pt"),
        os.path.join(PROJECT_ROOT, "lukas_cross_flow_download", "cross_flow_building_block", "tw_direction_angle_wall_model_deployable.pt"),
    ),
}


class BuildingBlockModel:
    def __init__(self, tw_path, angle_path):
        self.tw_model = torch.jit.load(tw_path)
        self.angle_model = torch.jit.load(angle_path)
        self.tw_model.eval()
        self.angle_model.eval()

    def predict(self, x):
        with torch.no_grad():
            tw_mean, tw_sigma = self.tw_model(x)
            angle_mean, angle_sigma = self.angle_model(x)
        return tw_mean, tw_sigma**2, angle_mean, angle_sigma**2


class BuildingBlockClassifier:
    def __init__(self, deployable_path):
        self.model = torch.jit.load(deployable_path)
        self.model.eval()

    def predict(self, x):
        with torch.no_grad():
            mean_w, cov_w = self.model(x)
        return mean_w, cov_w


def read_cz_uncertainty(wmlayers_dir):
    """Mirrors unified_wall_model.C's own readCZUncertainty() EXACTLY (same
    file, same "key value" line format, same missing-file->(0,0) default) -
    single source of truth (wmLayers_openfoam/cz_uncertainty.txt) read by
    BOTH the Python reference and the compiled C++, so they cannot silently
    drift apart the way two independently-hardcoded copies of the same two
    numbers could. See unified_wall_model.H's own comment on
    tau2Channel_/tau2Couette_ for what these mean."""
    path = os.path.join(wmlayers_dir, "cz_uncertainty.txt")
    tau2_channel, tau2_couette = 0.0, 0.0
    if not os.path.exists(path):
        return tau2_channel, tau2_couette
    with open(path) as f:
        for line in f:
            parts = line.split()
            if len(parts) != 2:
                continue
            key, value = parts
            if key == "tau2_channel":
                tau2_channel = float(value)
            elif key == "tau2_couette":
                tau2_couette = float(value)
    return tau2_channel, tau2_couette


class CrossflowWallModel:
    def __init__(self, classifier_path=None):
        if classifier_path is None:
            classifier_path = os.path.join(
                PROJECT_ROOT, "V1 classifier no backflow couette", "classifier_no_backflow_deployable.pt"
            )

        self.distinct_models = {
            name: BuildingBlockModel(tw_path, angle_path)
            for name, (tw_path, angle_path) in MODEL_PATHS.items()
        }
        self.classifier = BuildingBlockClassifier(classifier_path)

        tau2_channel, tau2_couette = read_cz_uncertainty(os.path.join(THIS_DIR, "wmLayers_openfoam"))
        # Same regime mapping as unified_wall_model.C: laminar gets neither
        # (closed-form, no ODT/(C,Z) dependence); freestream/crossflow reuse
        # channel's value (freestream aliases channel's prediction outright;
        # crossflow shares channel's own calibrated C,Z - confirmed via each
        # building block's own generator.py).
        self.tau2_by_label = {
            "laminar": 0.0,
            "freestream": tau2_channel,
            "channel": tau2_channel,
            "couette": tau2_couette,
            "crossflow": tau2_channel,
        }

    def _model_for_label(self, label):
        resolved = MODEL_ALIAS.get(label, label)
        return self.distinct_models[resolved]

    def blend_predict(self, raw_inputs):
        x = torch.as_tensor(raw_inputs, dtype=torch.float32)
        if x.ndim == 1:
            x = x.unsqueeze(0)

        w_mean, w_cov = self.classifier.predict(x)

        block_means_tw, block_vars_tw = [], []
        block_means_angle, block_vars_angle = [], []
        for label in LABEL_NAMES:
            model = self._model_for_label(label)
            tw_mean, tw_var, angle_mean, angle_var = model.predict(x)
            # (C,Z)-calibration-driven epistemic uncertainty, added onto
            # this regime's own NN-predicted variance BEFORE _combine() runs
            # - see unified_wall_model.C's identical comment. tw only, NOT
            # angle_var (out of scope, matching what was actually asked for).
            tw_var = tw_var + self.tau2_by_label[label]
            block_means_tw.append(tw_mean)
            block_vars_tw.append(tw_var)
            block_means_angle.append(angle_mean)
            block_vars_angle.append(angle_var)

        mu_tw = torch.stack(block_means_tw, dim=1)
        sigma2_tw = torch.stack(block_vars_tw, dim=1)
        mu_angle = torch.stack(block_means_angle, dim=1)
        sigma2_angle = torch.stack(block_vars_angle, dim=1)

        tw_mean, tw_var = self._combine(w_mean, w_cov, mu_tw, sigma2_tw)
        angle_mean, angle_var = self._combine_circular(w_mean, mu_angle, sigma2_angle)

        return {
            "tw_mean": tw_mean.numpy(),
            "tw_sigma": torch.sqrt(tw_var).numpy(),
            "angle_mean": angle_mean.numpy(),
            "angle_sigma": torch.sqrt(angle_var).numpy(),
            "weights_mean": w_mean.numpy(),
            "weights_sigma": torch.sqrt(torch.diagonal(w_cov, dim1=1, dim2=2)).numpy(),
            "weights_cov": w_cov.numpy(),
        }

    @staticmethod
    def _combine(w_mean, w_cov, mu, sigma2):
        combined_mean = torch.sum(w_mean * mu, dim=1)
        var_w = torch.diagonal(w_cov, dim1=1, dim2=2)
        term_a = torch.sum((var_w + w_mean**2) * sigma2, dim=1)
        mu_row = mu.unsqueeze(1)
        term_b = torch.bmm(torch.bmm(mu_row, w_cov), mu_row.transpose(1, 2)).squeeze(-1).squeeze(-1)
        combined_var = term_a + term_b
        return combined_mean, combined_var

    @staticmethod
    def _combine_circular(w_mean, mu, sigma2):
        sigma = torch.sqrt(sigma2)
        R_i = torch.exp(-sigma**2 / 2.0)
        Zx_i = R_i * torch.cos(mu)
        Zy_i = R_i * torch.sin(mu)
        Zx = torch.sum(w_mean * Zx_i, dim=1)
        Zy = torch.sum(w_mean * Zy_i, dim=1)
        combined_mean = torch.atan2(Zy, Zx)
        R_bar = torch.clamp(torch.sqrt(Zx**2 + Zy**2), max=1.0 - 1e-12)
        combined_var = -2.0 * torch.log(R_bar)
        return combined_mean, combined_var
