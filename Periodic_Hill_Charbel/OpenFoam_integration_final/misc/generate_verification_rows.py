"""Picks real rows spanning all 5 regimes from '../V1 classifier no backflow
couette/'s own training data (classifier_training_data_no_backflow.npy) and
runs them through reference_blend_predict.py's CrossflowWallModel - the
Python reference the compiled C++ unifiedWallModel gets checked against in
test_main.C. Same role as '../V1 classifier updated Re selective data/
generate_verification_rows.py', adapted to this package's own data/model
paths (see reference_blend_predict.py's own docstring for what changed).
"""
import os
import numpy as np

THIS_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECT_ROOT = os.path.dirname(THIS_DIR)

from reference_blend_predict import CrossflowWallModel, LABEL_NAMES  # noqa: E402

data = np.load(os.path.join(PROJECT_ROOT, "V1 classifier no backflow couette", "classifier_training_data_no_backflow.npy"))

rng = np.random.default_rng(7)
rows = []
for label in range(5):
    idx = np.where(data[:, 8] == label)[0]
    pick = rng.choice(idx)
    rows.append(data[pick])
rows = np.array(rows)

X = rows[:, :6]

model = CrossflowWallModel()
result = model.blend_predict(X)

with open(os.path.join(THIS_DIR, "verification_rows.txt"), "w") as f:
    f.write("# Re_center Re_top ax_n gamma12 lambda_roc tke_n | "
            "tw_mean tw_sigma angle_mean angle_sigma | "
            "w_laminar w_freestream w_channel w_couette w_crossflow (true_label)\n")
    for i in range(len(X)):
        w = result["weights_mean"][i]
        f.write(
            " ".join(f"{v:.17g}" for v in X[i]) + " | "
            f"{result['tw_mean'][i]:.17g} {result['tw_sigma'][i]:.17g} "
            f"{result['angle_mean'][i]:.17g} {result['angle_sigma'][i]:.17g} | "
            + " ".join(f"{v:.17g}" for v in w)
            + f" ({LABEL_NAMES[int(rows[i, 8])]})\n"
        )

print("Wrote verification_rows.txt")
print("\nReference (Python) predictions:")
for i in range(len(X)):
    true_label = LABEL_NAMES[int(rows[i, 8])]
    print(f"  row {i} (true={true_label}): tw={result['tw_mean'][i]:.6f}+-{result['tw_sigma'][i]:.6f}  "
          f"angle={result['angle_mean'][i]:.6f}+-{result['angle_sigma'][i]:.6f}  "
          f"weights={np.round(result['weights_mean'][i], 4)}")
