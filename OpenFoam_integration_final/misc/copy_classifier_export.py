"""Copies wmLayers_openfoam/classifier/ verbatim from '../V1 classifier no
backflow couette/wmLayers_openfoam/classifier' - the 5th and last piece of
this package's wmLayers_openfoam/ tree (the other 4 come from export_
regimes_for_openfoam.py). Kept as its own tiny script, not folded into that
one, because the classifier's own training/export in that directory is the
long-running step this package depends on finishing first - running this
separately makes that dependency explicit rather than silently baked into a
bigger script.

Run '../V1 classifier no backflow couette/train_classifier.py' then
'../V1 classifier no backflow couette/export_classifier_for_openfoam.py'
first; this script only copies their already-written output.
"""
import os
import shutil

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(THIS_DIR)
SRC = os.path.join(PROJECT_ROOT, "V1 classifier no backflow couette", "wmLayers_openfoam", "classifier")
DST = os.path.join(THIS_DIR, "wmLayers_openfoam", "classifier")


def main():
    if not os.path.isdir(SRC):
        raise FileNotFoundError(
            f"{SRC} does not exist yet - run this directory's own train_classifier.py "
            "and export_classifier_for_openfoam.py first"
        )
    if os.path.isdir(DST):
        shutil.rmtree(DST)
    shutil.copytree(SRC, DST)
    print(f"Copied {SRC} -> {DST}")


if __name__ == "__main__":
    main()
