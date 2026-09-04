"""Assembles wmLayers_openfoam/{laminar,channel,couette,crossflow}/ for this
final deployment by COPYING each regime's already-built, already-verified
export verbatim from its own canonical training directory - no re-deriving
from .pt files here, since every source below already ran its own
export_ensemble_for_openfoam.py-equivalent and the .bin/layer_metadata.txt
format is a plain file copy away from being reused (same "copy verbatim"
pattern '../V1 classifier updated Re selective data/export_regimes_for_
openfoam.py' already used for laminar specifically - extended here to all 4
regimes since none of them changed architecture, only weights).

Sources - each one deliberately picked, not just "whatever was lying around":

  laminar  <- '../V1 classifier updated Re selective data/wmLayers_openfoam/
              laminar' - unretrained (w=0 always, old/new Re formulas
              coincide exactly, see that directory's own README §7).
  channel  <- same directory's wmLayers_openfoam/channel - unretrained this
              round (last retrained for the NEW Re formula, no further
              change needed).
  couette  <- '../updated couette no backflow building block/wmLayers_
              openfoam/couette' - THE retrain this package exists for.
              couette_data_aug_12 filtered by purge.py (15 of 56 nominal
              conditions removed for genuine backflow), fixing
              tw_direction_angle's near-+-pi wraparound blowup (MAE
              0.278->0.040, see that directory's train_couette.log).
  crossflow <- 'lukas_cross_flow_download/cross_flow_building_block/
              wmLayers_openfoam/crossflow' - the FULL-SCALE (all 4 Re_tau)
              Aug-11 retrain, confirmed byte-identical to
              'crossflow_building_block_aug_11/'s own copy (md5sum). NOTE:
              this is a deliberate correction relative to '../V1 classifier
              updated Re selective data/', whose own wmLayers_openfoam/
              crossflow was exported from the OLDER, reduced-scale (2-of-4
              Re_tau) 'crossflow building block updated Re/' - a real,
              flagged-but-not-fixed inconsistency in that directory (its own
              README §9) between what the classifier was trained to
              recognize as "crossflow" (full-scale Aug-11 data, via
              collect_training_data.py) and which regression model actually
              got deployed for it (the old reduced-scale one). This package
              uses the full-scale model for both, closing that gap.

classifier <- '../V1 classifier no backflow couette/wmLayers_openfoam/
              classifier' - copied by a separate step in build_final_
              package.py (that directory's classifier training must finish
              first; this script only handles the 4 regression building
              blocks, which have no such dependency).
"""
import os
import shutil

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(THIS_DIR)
OUT_ROOT = os.path.join(THIS_DIR, "wmLayers_openfoam")

SOURCES = {
    "laminar": os.path.join(PROJECT_ROOT, "V1 classifier updated Re selective data", "wmLayers_openfoam", "laminar"),
    "channel": os.path.join(PROJECT_ROOT, "V1 classifier updated Re selective data", "wmLayers_openfoam", "channel"),
    "couette": os.path.join(PROJECT_ROOT, "updated couette no backflow building block", "wmLayers_openfoam", "couette"),
    "crossflow": os.path.join(PROJECT_ROOT, "lukas_cross_flow_download", "cross_flow_building_block", "wmLayers_openfoam", "crossflow"),
}


def main():
    for name, src in SOURCES.items():
        if not os.path.isdir(src):
            raise FileNotFoundError(f"{name}: source not found: {src}")
        dst = os.path.join(OUT_ROOT, name)
        if os.path.isdir(dst):
            shutil.rmtree(dst)
        shutil.copytree(src, dst)
        angle_type_path = os.path.join(dst, "tw_direction_angle", "angle_type.txt")
        angle_type = open(angle_type_path).read().strip() if os.path.exists(angle_type_path) else "MISSING"
        print(f"  {name}: copied from {src} (angle_type={angle_type})")
    print(f"\nAll 4 regression regimes copied into {OUT_ROOT}")
    print("NOTE: wmLayers_openfoam/classifier/ is NOT written by this script - "
          "see build_final_package.py / README.md step for that piece.")


if __name__ == "__main__":
    main()
