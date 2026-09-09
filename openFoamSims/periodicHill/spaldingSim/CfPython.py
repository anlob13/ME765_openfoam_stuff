#!/usr/bin/env python3

"""
Periodic-hill wall-shear / skin-friction analysis.

Reads:
    constant/polyMesh/points
    constant/polyMesh/faces
    constant/polyMesh/boundary
    <time>/wallShearStressMean

Calculates:
    - wall-face centres
    - wall-face normals
    - wall tangent direction
    - tangential wall shear
    - streamwise-averaged wall shear
    - skin-friction coefficient Cf
    - separation / reattachment locations
    - CSV and PNG output

Important:
OpenFOAM wallShearStress is a KINEMATIC quantity:
    tau / rho

Therefore:
    Cf = (tau_w/rho) / (0.5 Ub^2)

The sign convention is adjusted so that conventional Cf is:
    positive  -> attached flow
    negative  -> separated flow

For the supplied periodic-hill case:
    h  = 0.028 m
    Ub = 1.0 m/s

The streamwise origin is automatically shifted so that the
minimum wall x-coordinate corresponds to x/h = 0.
"""

import re
from pathlib import Path

import numpy as np
import matplotlib

# Non-interactive backend: important for SSH / cluster runs.
matplotlib.use("Agg")

import matplotlib.pyplot as plt


# ============================================================================
# USER SETTINGS
# ============================================================================

UB = 1.0          # Bulk/reference velocity [m/s]
H = 0.028         # Hill height / reference height [m]

PATCH_NAME = "hills"

CASES = [
    {
        "name": "Spalding",
        "dir": ".",
        "time": "0.5",
    },
]

# ---------------------------------------------------------------------------
# Coordinate origin
#
# None:
#     Automatically use the minimum x-coordinate of the hills patch.
#
# A number, e.g. 0.028:
#     Explicitly use that physical x-coordinate as x/h = 0.
#
# For your current case, automatic mode is recommended initially.
# ---------------------------------------------------------------------------
X_ORIGIN = None

# ---------------------------------------------------------------------------
# OpenFOAM wallShearStress sign convention.
#
# For the supplied lower-wall field, the raw OpenFOAM tangential stress
# points opposite to the conventional positive skin-friction direction.
#
# Therefore use:
#     Cf = -tau_t / (0.5 Ub^2)
#
# Set this to +1.0 if you want the raw OpenFOAM sign instead.
# ---------------------------------------------------------------------------
CF_SIGN = -1.0

# Tolerance for grouping faces at the same streamwise x location.
# Increase slightly if numerical roundoff causes duplicate x stations.
X_TOL = 1e-10

# Minimum allowed face area.
AREA_TOL = 1e-15

# Output names.
PNG_OUTPUT = "periodic_hill_Cf.png"
CSV_OUTPUT = "periodic_hill_Cf.csv"


# ============================================================================
# BASIC FILE READING
# ============================================================================

def read_file(path):
    path = Path(path)

    with path.open("r", errors="replace") as f:
        return f.read()


def remove_comments(text):
    # C-style block comments
    text = re.sub(r"/\*.*?\*/", "", text, flags=re.S)

    # C++ single-line comments
    text = re.sub(r"//.*", "", text)

    return text


# ============================================================================
# READ OPENFOAM points
# ============================================================================

def read_points(path):
    text = remove_comments(read_file(path))

    # Find the first integer followed by "(".
    match = re.search(r"(\d+)\s*\(", text)

    if not match:
        raise RuntimeError(
            f"Could not find point count in {path}"
        )

    n_points = int(match.group(1))

    data = text[match.end():]

    tokens = re.findall(
        r"\(\s*"
        r"([-+0-9.eE]+)\s+"
        r"([-+0-9.eE]+)\s+"
        r"([-+0-9.eE]+)"
        r"\s*\)",
        data,
    )

    if len(tokens) < n_points:
        raise RuntimeError(
            f"{path}: expected {n_points} points, "
            f"found {len(tokens)}"
        )

    points = np.array(
        [
            [float(x), float(y), float(z)]
            for x, y, z in tokens[:n_points]
        ],
        dtype=float,
    )

    return points


# ============================================================================
# READ OPENFOAM faces
# ============================================================================

def read_faces(path):
    text = remove_comments(read_file(path))

    # First integer followed by "(" is the number of faces.
    match = re.search(r"(\d+)\s*\(", text)

    if not match:
        raise RuntimeError(
            f"Could not find face count in {path}"
        )

    n_faces = int(match.group(1))

    data = text[match.end():]

    # Standard OpenFOAM format:
    #
    # 4(1 2 3 4)
    #
    # or
    #
    # 4
    # (
    # 1 2 3 4
    # )
    #
    # We explicitly capture one face at a time.

    pattern = re.compile(
        r"(\d+)\s*\(\s*([0-9\s]+?)\s*\)"
    )

    matches = pattern.findall(data)

    if len(matches) < n_faces:
        raise RuntimeError(
            f"{path}: expected {n_faces} faces, "
            f"found {len(matches)}"
        )

    faces = []

    for n_str, body in matches[:n_faces]:

        n_expected = int(n_str)

        ids = np.fromstring(
            body,
            sep=" ",
            dtype=int,
        )

        if len(ids) != n_expected:
            raise RuntimeError(
                f"Invalid face in {path}: "
                f"expected {n_expected} vertices, "
                f"found {len(ids)}"
            )

        if len(ids) < 3:
            raise RuntimeError(
                f"Degenerate face in {path}: {body}"
            )

        faces.append(ids)

    return faces


# ============================================================================
# READ PATCH INFORMATION
# ============================================================================

def read_patch_info(path, patch_name):
    text = remove_comments(read_file(path))

    pattern = (
        rf"\b{re.escape(patch_name)}\s*"
        rf"\{{(.*?)\}}"
    )

    match = re.search(
        pattern,
        text,
        flags=re.S,
    )

    if not match:
        raise RuntimeError(
            f"Patch '{patch_name}' not found in {path}"
        )

    block = match.group(1)

    start_match = re.search(
        r"\bstartFace\s+(\d+)\s*;",
        block,
    )

    n_match = re.search(
        r"\bnFaces\s+(\d+)\s*;",
        block,
    )

    if not start_match or not n_match:
        raise RuntimeError(
            f"Could not read startFace/nFaces "
            f"for patch '{patch_name}'"
        )

    start_face = int(start_match.group(1))
    n_faces = int(n_match.group(1))

    return start_face, n_faces


# ============================================================================
# READ wallShearStressMean
# ============================================================================

def read_wall_shear(path, patch_name):
    text = remove_comments(read_file(path))

    pattern = (
        rf"\b{re.escape(patch_name)}\s*"
        rf"\{{(.*?)\}}"
    )

    match = re.search(
        pattern,
        text,
        flags=re.S,
    )

    if not match:
        raise RuntimeError(
            f"Patch '{patch_name}' not found in {path}"
        )

    block = match.group(1)

    match_list = re.search(
        r"nonuniform\s+List<vector>\s+(\d+)\s*\(",
        block,
        flags=re.S,
    )

    if not match_list:
        raise RuntimeError(
            f"Could not find nonuniform vector list "
            f"for patch '{patch_name}' in {path}"
        )

    n = int(match_list.group(1))

    data = block[match_list.end():]

    tokens = re.findall(
        r"\(\s*"
        r"([-+0-9.eE]+)\s+"
        r"([-+0-9.eE]+)\s+"
        r"([-+0-9.eE]+)"
        r"\s*\)",
        data,
    )

    if len(tokens) < n:
        raise RuntimeError(
            f"{path}: expected {n} wall-shear vectors, "
            f"found {len(tokens)}"
        )

    tau = np.array(
        [
            [float(a), float(b), float(c)]
            for a, b, c in tokens[:n]
        ],
        dtype=float,
    )

    return tau


# ============================================================================
# FACE GEOMETRY
# ============================================================================

def face_geometry(points, face):
    vertices = points[face]

    # Geometric centre.
    centre = np.mean(vertices, axis=0)

    # Polygon area vector by triangulation.
    area_vector = np.zeros(3)

    p0 = vertices[0]

    for i in range(1, len(vertices) - 1):

        p1 = vertices[i]
        p2 = vertices[i + 1]

        area_vector += (
            0.5 * np.cross(
                p1 - p0,
                p2 - p0,
            )
        )

    area = np.linalg.norm(area_vector)

    if area < AREA_TOL:
        raise RuntimeError(
            "Encountered zero-area wall face."
        )

    normal = area_vector / area

    # The hills patch is the lower wall.
    #
    # We explicitly choose the wall-normal direction so that:
    #
    #     ny < 0
    #
    # This makes the tangent below point approximately +x.
    if normal[1] > 0.0:
        normal *= -1.0

    # Tangent in x-y plane:
    #
    # n = (nx, ny, 0)
    #
    # t = (-ny, nx, 0)
    #
    tangent = np.array(
        [
            -normal[1],
            normal[0],
            0.0,
        ]
    )

    tangent_magnitude = np.linalg.norm(tangent)

    if tangent_magnitude < AREA_TOL:
        raise RuntimeError(
            "Unable to construct wall tangent."
        )

    tangent /= tangent_magnitude

    # Force streamwise tangent direction to be +x.
    if tangent[0] < 0.0:
        tangent *= -1.0

    return centre, normal, tangent, area


# ============================================================================
# GROUP VALUES BY STREAMWISE LOCATION
# ============================================================================

def average_by_x(x, values, weights, tol):
    """
    Area-weighted average of values for wall faces having the same x.

    Returns:
        x_out
        value_out
        area_out
        count_out
    """

    order = np.argsort(x)

    x_sorted = x[order]
    v_sorted = values[order]
    w_sorted = weights[order]

    x_out = []
    v_out = []
    area_out = []
    count_out = []

    start = 0

    for i in range(1, len(x_sorted)):

        if abs(x_sorted[i] - x_sorted[start]) > tol:

            x_group = x_sorted[start:i]
            v_group = v_sorted[start:i]
            w_group = w_sorted[start:i]

            total_area = np.sum(w_group)

            x_out.append(
                np.average(
                    x_group,
                    weights=w_group,
                )
            )

            v_out.append(
                np.average(
                    v_group,
                    weights=w_group,
                )
            )

            area_out.append(total_area)
            count_out.append(i - start)

            start = i

    # Last group.
    x_group = x_sorted[start:]
    v_group = v_sorted[start:]
    w_group = w_sorted[start:]

    total_area = np.sum(w_group)

    x_out.append(
        np.average(
            x_group,
            weights=w_group,
        )
    )

    v_out.append(
        np.average(
            v_group,
            weights=w_group,
        )
    )

    area_out.append(total_area)
    count_out.append(
        len(x_sorted) - start
    )

    return (
        np.array(x_out),
        np.array(v_out),
        np.array(area_out),
        np.array(count_out),
    )


# ============================================================================
# ZERO CROSSINGS
# ============================================================================

def find_zero_crossings(x, y):
    roots = []

    for i in range(len(y) - 1):

        y1 = y[i]
        y2 = y[i + 1]

        if y1 == 0.0:

            roots.append(x[i])

        elif y1 * y2 < 0.0:

            root = (
                x[i]
                - y1
                * (x[i + 1] - x[i])
                / (y2 - y1)
            )

            roots.append(root)

    return roots


# ============================================================================
# PROCESS ONE CASE
# ============================================================================

def process_case(case):

    case_dir = Path(case["dir"])

    mesh_dir = (
        case_dir
        / "constant"
        / "polyMesh"
    )

    field_file = (
        case_dir
        / str(case["time"])
        / "wallShearStressMean"
    )

    points_file = mesh_dir / "points"
    faces_file = mesh_dir / "faces"
    boundary_file = mesh_dir / "boundary"

    required_files = [
        points_file,
        faces_file,
        boundary_file,
        field_file,
    ]

    for path in required_files:

        if not path.exists():

            raise FileNotFoundError(
                f"Missing file:\n{path}"
            )

    print()
    print("=" * 78)
    print(case["name"])
    print("=" * 78)

    print(f"Reading:")
    print(f"  points : {points_file}")
    print(f"  faces  : {faces_file}")
    print(f"  boundary: {boundary_file}")
    print(f"  shear  : {field_file}")

    # ------------------------------------------------------------------------
    # Read mesh and field.
    # ------------------------------------------------------------------------

    points = read_points(points_file)
    faces = read_faces(faces_file)

    start_face, n_patch_faces = read_patch_info(
        boundary_file,
        PATCH_NAME,
    )

    tau = read_wall_shear(
        field_file,
        PATCH_NAME,
    )

    print()
    print("RAW DATA")
    print("-" * 78)

    print(
        f"points                  = {len(points)}"
    )

    print(
        f"mesh faces              = {len(faces)}"
    )

    print(
        f"hills startFace         = {start_face}"
    )

    print(
        f"hills nFaces            = {n_patch_faces}"
    )

    print(
        f"wallShearStress vectors = {len(tau)}"
    )

    if len(tau) != n_patch_faces:

        raise RuntimeError(
            "\nMismatch between boundary patch and "
            "wallShearStressMean:\n"
            f"  hills patch = {n_patch_faces}\n"
            f"  shear field = {len(tau)}"
        )

    # ------------------------------------------------------------------------
    # Build geometry.
    # ------------------------------------------------------------------------

    centres = np.zeros(
        (n_patch_faces, 3)
    )

    normals = np.zeros(
        (n_patch_faces, 3)
    )

    tangents = np.zeros(
        (n_patch_faces, 3)
    )

    areas = np.zeros(
        n_patch_faces
    )

    for i in range(n_patch_faces):

        face_index = start_face + i

        if face_index >= len(faces):

            raise RuntimeError(
                f"Patch face index {face_index} "
                f"is outside faces array."
            )

        (
            centre,
            normal,
            tangent,
            area,
        ) = face_geometry(
            points,
            faces[face_index],
        )

        centres[i] = centre
        normals[i] = normal
        tangents[i] = tangent
        areas[i] = area

    # ------------------------------------------------------------------------
    # Geometry diagnostics.
    # ------------------------------------------------------------------------

    x_raw = centres[:, 0]
    y_raw = centres[:, 1]
    z_raw = centres[:, 2]

    print()
    print("WALL GEOMETRY")
    print("-" * 78)

    print(
        f"x range = "
        f"{x_raw.min():.10e} ... "
        f"{x_raw.max():.10e} m"
    )

    print(
        f"y range = "
        f"{y_raw.min():.10e} ... "
        f"{y_raw.max():.10e} m"
    )

    print(
        f"z range = "
        f"{z_raw.min():.10e} ... "
        f"{z_raw.max():.10e} m"
    )

    print(
        f"face area range = "
        f"{areas.min():.10e} ... "
        f"{areas.max():.10e} m^2"
    )

    # ------------------------------------------------------------------------
    # Raw shear diagnostics.
    #
    # The supplied field has dimensions [0 2 -2 0 0 0 0], i.e. tau/rho.
    # ------------------------------------------------------------------------

    tau_x = tau[:, 0]
    tau_y = tau[:, 1]
    tau_z = tau[:, 2]

    print()
    print("RAW wallShearStressMean")
    print("-" * 78)

    print(
        f"tau_x/rho min = "
        f"{tau_x.min(): .8e}"
    )

    print(
        f"tau_x/rho max = "
        f"{tau_x.max(): .8e}"
    )

    print(
        f"tau_x/rho mean = "
        f"{tau_x.mean(): .8e}"
    )

    print(
        f"tau_y/rho min = "
        f"{tau_y.min(): .8e}"
    )

    print(
        f"tau_y/rho max = "
        f"{tau_y.max(): .8e}"
    )

    print(
        f"tau_z/rho min = "
        f"{tau_z.min(): .8e}"
    )

    print(
        f"tau_z/rho max = "
        f"{tau_z.max(): .8e}"
    )

    # ------------------------------------------------------------------------
    # Project shear onto wall tangent.
    #
    # Since wallShearStress should already be tangential, this is also a
    # useful consistency check.
    # ------------------------------------------------------------------------

    tau_t = np.einsum(
        "ij,ij->i",
        tau,
        tangents,
    )

    tau_normal = np.einsum(
        "ij,ij->i",
        tau,
        normals,
    )

    print()
    print("SHEAR PROJECTION DIAGNOSTICS")
    print("-" * 78)

    print(
        f"tau_t/rho min = "
        f"{tau_t.min(): .8e}"
    )

    print(
        f"tau_t/rho max = "
        f"{tau_t.max(): .8e}"
    )

    print(
        f"tau_t/rho mean = "
        f"{tau_t.mean(): .8e}"
    )

    print(
        f"|tau_normal/rho| max = "
        f"{np.max(np.abs(tau_normal)): .8e}"
    )

    tau_mag = np.linalg.norm(tau, axis=1)

    normal_fraction = (
        np.abs(tau_normal)
        / np.maximum(tau_mag, 1e-30)
    )

    print(
        f"|tau_normal| / |tau| max = "
        f"{normal_fraction.max(): .8e}"
    )

    # If wallShearStress is implemented correctly, the normal component
    # should be very small compared to the tangential component.

    # ------------------------------------------------------------------------
    # Group by actual wall x position.
    # ------------------------------------------------------------------------

    (
        x_unique,
        tau_t_mean,
        area_x,
        counts,
    ) = average_by_x(
        x_raw,
        tau_t,
        areas,
        X_TOL,
    )

    print()
    print("STREAMWISE GROUPING")
    print("-" * 78)

    print(
        f"x stations = {len(x_unique)}"
    )

    print(
        f"faces per x station = "
        f"{counts.min():.0f} ... "
        f"{counts.max():.0f}"
    )

    if np.any(counts != counts[0]):

        print(
            "WARNING: number of wall faces per x station "
            "is not constant."
        )

    # ------------------------------------------------------------------------
    # Streamwise coordinate normalization.
    # ------------------------------------------------------------------------

    if X_ORIGIN is None:

        x0 = np.min(x_unique)

        print()
        print(
            "X ORIGIN: automatic"
        )

        print(
            f"Using x0 = min(x) = "
            f"{x0:.10e} m"
        )

    else:

        x0 = float(X_ORIGIN)

        print()
        print(
            "X ORIGIN: user supplied"
        )

        print(
            f"Using x0 = "
            f"{x0:.10e} m"
        )

    xh = (
        x_unique - x0
    ) / H

    print(
        f"x/h range = "
        f"{xh.min():.8f} ... "
        f"{xh.max():.8f}"
    )

    # ------------------------------------------------------------------------
    # Skin-friction coefficient.
    #
    # tau_t is kinematic wall shear stress:
    #
    #       tau/rho
    #
    # Conventional periodic-hill Cf:
    #
    #       Cf = - (tau_t/rho) / (0.5 Ub^2)
    #
    # for this OpenFOAM sign convention.
    # ------------------------------------------------------------------------

    denominator = 0.5 * UB**2

    Cf = (
        CF_SIGN
        * tau_t_mean
        / denominator
    )

    print()
    print("SKIN FRICTION")
    print("-" * 78)

    print(
        f"Ub = {UB:.8f} m/s"
    )

    print(
        f"H  = {H:.8f} m"
    )

    print(
        f"0.5*Ub^2 = "
        f"{denominator:.8e} m^2/s^2"
    )

    print(
        f"Cf min = "
        f"{Cf.min(): .8e}"
    )

    print(
        f"Cf max = "
        f"{Cf.max(): .8e}"
    )

    print(
        f"Cf mean = "
        f"{Cf.mean(): .8e}"
    )

    # ------------------------------------------------------------------------
    # Print selected stations.
    # ------------------------------------------------------------------------

    print()
    print("SAMPLE Cf VALUES")
    print("-" * 78)

    sample_indices = np.linspace(
        0,
        len(xh) - 1,
        min(20, len(xh)),
        dtype=int,
    )

    print(
        f"{'x/h':>12s}"
        f"{'tau_t/rho':>18s}"
        f"{'Cf':>18s}"
        f"{'faces':>10s}"
    )

    for i in sample_indices:

        print(
            f"{xh[i]:12.6f}"
            f"{tau_t_mean[i]:18.8e}"
            f"{Cf[i]:18.8e}"
            f"{counts[i]:10.0f}"
        )

    # ------------------------------------------------------------------------
    # Separation / reattachment.
    # ------------------------------------------------------------------------

    roots = find_zero_crossings(
        xh,
        Cf,
    )

    print()
    print("=" * 78)
    print("SEPARATION / REATTACHMENT")
    print("=" * 78)

    if len(roots) == 0:

        print(
            "No Cf zero crossings detected."
        )

    else:

        for root in roots:

            print(
                f"Cf = 0 at x/h = "
                f"{root:.8f}"
            )

    # ------------------------------------------------------------------------
    # Return results.
    # ------------------------------------------------------------------------

    return {
        "name": case["name"],
        "x": x_unique,
        "xh": xh,
        "tau": tau_t_mean,
        "Cf": Cf,
        "centres": centres,
        "normals": normals,
        "tangents": tangents,
        "areas": areas,
        "rawTau": tau,
        "tauNormal": tau_normal,
    }


# ============================================================================
# MAIN
# ============================================================================

def main():

    results = []

    for case in CASES:

        result = process_case(case)

        results.append(result)

    # ------------------------------------------------------------------------
    # Plot Cf.
    # ------------------------------------------------------------------------

    plt.figure(
        figsize=(10, 5.5)
    )

    for result in results:

        plt.plot(
            result["xh"],
            result["Cf"],
            linewidth=1.8,
            label=result["name"],
        )

    plt.axhline(
        0.0,
        linestyle="--",
        linewidth=1.0,
    )

    plt.xlabel(
        r"$x/h$"
    )

    plt.ylabel(
        r"$C_f$"
    )

    plt.xlim(
        0.0,
        max(
            9.0,
            max(r["xh"].max() for r in results),
        ),
    )

    plt.title(
        "Periodic-hill skin friction"
    )

    plt.grid(
        True,
        alpha=0.3,
    )

    plt.legend()

    plt.tight_layout()

    plt.savefig(
        PNG_OUTPUT,
        dpi=250,
        bbox_inches="tight",
    )

    plt.close()

    print()
    print(
        f"Saved plot: {PNG_OUTPUT}"
    )

    # ------------------------------------------------------------------------
    # Save CSV.
    # ------------------------------------------------------------------------

    # Use the first case as the reference x grid.
    xref = results[0]["xh"]

    columns = [
        xref
    ]

    headers = [
        "x_over_h"
    ]

    for result in results:

        cf = result["Cf"]

        if (
            len(result["xh"]) == len(xref)
            and np.allclose(
                result["xh"],
                xref,
            )
        ):

            cf_out = cf

        else:

            cf_out = np.interp(
                xref,
                result["xh"],
                cf,
                left=np.nan,
                right=np.nan,
            )

        columns.append(
            cf_out
        )

        headers.append(
            "Cf_" + result["name"]
        )

    data = np.column_stack(
        columns
    )

    np.savetxt(
        CSV_OUTPUT,
        data,
        delimiter=",",
        header=",".join(headers),
        comments="",
    )

    print(
        f"Saved CSV: {CSV_OUTPUT}"
    )


# ============================================================================
# ENTRY POINT
# ============================================================================

if __name__ == "__main__":
    main()