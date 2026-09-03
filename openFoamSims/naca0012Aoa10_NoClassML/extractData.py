#!/usr/bin/env python3
"""
Extract aerodynamic data from ONE 3D incompressible OpenFOAM NACA 0012 case.

Usage:
    python3 extract_naca3d_incompressible.py naca0012Aoa10_ML

Outputs:
    Cp.csv
        x/c, z/c, Cp

    Cf.csv
        x/c, z/c, Cf

    wallShearStress.csv
        Full 3D wall-shear data

    raw_surface.csv
        Full 3D surface pressure data

    summary.csv
        Integrated forces and coefficients

    surface_coefficients.png
        Cp and Cf distributions

Coordinate convention:
    x = chord
    y = span
    z = airfoil thickness

Incompressible OpenFOAM:
    p               = kinematic pressure [m2/s2]
    wallShearStress = kinematic shear [m2/s2]

Physical quantities:
    p_phys = rho*p
    tau_w  = rho*wallShearStress

Cp = (p - p_inf)/(0.5*U_inf^2)
"""

import os
import sys
import re
import csv
import numpy as np
import matplotlib.pyplot as plt
import pyvista as pv

# ======================================================================
# CONFIGURATION
# ======================================================================

U_INF = 43.824
RHO = 1.225
P_INF = 0.0

CHORD = 1.0
REFERENCE_AREA = 0.2

AOA = None

AIRFOIL_PATCH_NAMES = [
    "aerofoil",
    "airfoil",
    "wing",
]

WALL_SHEAR_IS_KINEMATIC = True

N_X_STATIONS = 500
DUPLICATE_ROUND_DIGITS = 8

OUTPUT_DIR = "comparison_results"

# ======================================================================
# COMMAND LINE
# ======================================================================

def get_case_path():
    if len(sys.argv) != 2:
        print(
            "\nUsage:\n"
            "    python3 extract_naca3d_incompressible.py "
            "<simulation_folder>\n\n"
            "Example:\n"
            "    python3 extract_naca3d_incompressible.py "
            "naca0012Aoa10_ML\n"
        )
        sys.exit(1)

    case_path = sys.argv[1]

    if not os.path.isdir(case_path):
        raise RuntimeError(
            f"Simulation folder does not exist:\n    {case_path}"
        )

    return os.path.abspath(case_path)


# ======================================================================
# AOA
# ======================================================================

def infer_aoa(case_path):
    name = os.path.basename(os.path.normpath(case_path))
    patterns = [
        r"Aoa(-?\d+(?:\.\d+)?)",
        r"AoA(-?\d+(?:\.\d+)?)",
        r"aoa(-?\d+(?:\.\d+)?)",
    ]

    for pattern in patterns:
        match = re.search(pattern, name)
        if match:
            return float(match.group(1))

    return None


# ======================================================================
# BASIC
# ======================================================================

def q_inf():
    return 0.5 * RHO * U_INF**2


def q_inf_kinematic():
    return 0.5 * U_INF**2


def ensure_foam_stub(case_path):
    name = os.path.basename(os.path.normpath(case_path))
    foam = os.path.join(case_path, f"{name}.foam")

    if not os.path.exists(foam):
        open(foam, "w").close()

    return foam


# ======================================================================
# READ OPENFOAM
# ======================================================================

def read_case(case_path):
    foam = ensure_foam_stub(case_path)

    print(f"\nReading:\n    {foam}")

    reader = pv.OpenFOAMReader(foam)

    if not reader.time_values:
        raise RuntimeError("No OpenFOAM time values found.")

    print("\nAvailable times:")
    print(reader.time_values)

    final_time = reader.time_values[-1]
    print(f"\nUsing final time: {final_time}")

    reader.set_active_time_value(final_time)
    reader.cell_to_point_creation = True

    mesh = reader.read()

    internal = mesh["internalMesh"] if "internalMesh" in mesh.keys() else None
    boundary = mesh["boundary"] if "boundary" in mesh.keys() else None

    return internal, boundary, final_time


# ======================================================================
# FIND AIRFOIL
# ======================================================================

def get_airfoil_patch(boundary):
    if boundary is None:
        raise RuntimeError("No boundary mesh found.")

    print("\nBoundary patches:")

    for name in boundary.keys():
        try:
            patch = boundary[name]
            print(
                f"  {name:25s} "
                f"points={patch.n_points:8d} "
                f"cells={patch.n_cells:8d} "
                f"bounds={patch.bounds}"
            )
        except Exception:
            pass

    for name in AIRFOIL_PATCH_NAMES:
        if name in boundary.keys():
            print(f"\nUsing airfoil patch: {name}")
            return boundary[name]

    raise RuntimeError("Could not find the airfoil patch.")


# ======================================================================
# SURFACE GEOMETRY
# ======================================================================

def get_surface_geometry(patch):
    print("\nExtracting surface...")

    surf = patch.extract_surface().triangulate()

    points = np.asarray(surf.points, dtype=float)
    faces = np.asarray(surf.faces, dtype=np.int64)

    if len(faces) % 4 != 0:
        raise RuntimeError("Unexpected surface connectivity.")

    tri = faces.reshape(-1, 4)

    if not np.all(tri[:, 0] == 3):
        raise RuntimeError("Surface is not fully triangulated.")

    r0 = points[tri[:, 1]]
    r1 = points[tri[:, 2]]
    r2 = points[tri[:, 3]]

    Sf = 0.5 * np.cross(r1 - r0, r2 - r0)
    areas = np.linalg.norm(Sf, axis=1)
    centres = (r0 + r1 + r2) / 3.0

    good = areas > 1e-14

    if not np.any(good):
        raise RuntimeError("No valid surface faces found.")

    # Orient normals outward in the x-z plane.
    cx = np.mean(centres[good, 0])
    cz = np.mean(centres[good, 2])

    radial = np.column_stack([
        centres[:, 0] - cx,
        np.zeros(len(centres)),
        centres[:, 2] - cz,
    ])

    radial_mag = np.linalg.norm(radial, axis=1)
    valid = good & (radial_mag > 1e-14)

    alignment = (
        np.sum(Sf[valid] * radial[valid], axis=1)
        / (areas[valid] * radial_mag[valid])
    )

    if np.mean(alignment) < 0:
        print("  Reversing surface orientation.")
        Sf *= -1.0

    normals = np.zeros_like(Sf)
    normals[good] = Sf[good] / areas[good, None]

    print(f"  Surface triangles = {len(areas)}")
    print(f"  Surface area      = {np.sum(areas):.8e} m^2")
    print(
        f"  x range           = "
        f"{centres[:,0].min():.8e} to {centres[:,0].max():.8e}"
    )
    print(
        f"  y range           = "
        f"{centres[:,1].min():.8e} to {centres[:,1].max():.8e}"
    )
    print(
        f"  z range           = "
        f"{centres[:,2].min():.8e} to {centres[:,2].max():.8e}"
    )

    return surf, centres, Sf, normals, areas


# ======================================================================
# PRESSURE
# ======================================================================

def get_cell_pressure(poly):
    if "p" in poly.cell_data:
        return np.asarray(poly.cell_data["p"], dtype=float)

    if "p" in poly.point_data:
        converted = poly.point_data_to_cell_data()
        return np.asarray(converted.cell_data["p"], dtype=float)

    raise RuntimeError("Pressure field 'p' not found.")


# ======================================================================
# WALL SHEAR
# ======================================================================

def get_wall_shear(poly):
    if "wallShearStress" not in poly.cell_data:
        print("\nWARNING: wallShearStress not found.")
        return None

    tau = np.asarray(poly.cell_data["wallShearStress"], dtype=float)

    if tau.ndim != 2 or tau.shape[1] != 3:
        raise RuntimeError(
            "wallShearStress must be a 3-component vector field."
        )

    if WALL_SHEAR_IS_KINEMATIC:
        print("  Converting wallShearStress from kinematic to physical.")
        tau = RHO * tau

    return tau


# ======================================================================
# FORCE INTEGRATION
# ======================================================================

def compute_forces(p, tau, centres, Sf, normals, areas, aoa):
    pressure_phys = RHO * (p - P_INF)

    Fp = -np.sum(pressure_phys[:, None] * Sf, axis=0)

    if tau is None:
        Fv = np.zeros(3)
    else:
        tau_n = np.sum(tau * normals, axis=1)
        tau_t = tau - tau_n[:, None] * normals
        Fv = -np.sum(tau_t * areas[:, None], axis=0)

    Ft = Fp + Fv

    alpha = np.deg2rad(aoa)

    drag_dir = np.array([
        np.cos(alpha),
        0.0,
        np.sin(alpha),
    ])

    lift_dir = np.array([
        -np.sin(alpha),
        0.0,
        np.cos(alpha),
    ])

    Dp = np.dot(Fp, drag_dir)
    Lp = np.dot(Fp, lift_dir)

    Dv = np.dot(Fv, drag_dir)
    Lv = np.dot(Fv, lift_dir)

    D = np.dot(Ft, drag_dir)
    L = np.dot(Ft, lift_dir)

    denominator = q_inf() * REFERENCE_AREA

    result = {
        "Cl": L / denominator,
        "Cd": D / denominator,
        "Cl_pressure": Lp / denominator,
        "Cd_pressure": Dp / denominator,
        "Cl_viscous": Lv / denominator,
        "Cd_viscous": Dv / denominator,
        "Fx": Ft[0],
        "Fy": Ft[1],
        "Fz": Ft[2],
        "Fx_pressure": Fp[0],
        "Fy_pressure": Fp[1],
        "Fz_pressure": Fp[2],
        "Fx_viscous": Fv[0],
        "Fy_viscous": Fv[1],
        "Fz_viscous": Fv[2],
        "surface_area": np.sum(areas),
        "span": np.ptp(centres[:, 1]),
    }

    print("\nIntegrated forces")
    print("-" * 60)
    print(f"Pressure: Fx={Fp[0]:.8e} Fy={Fp[1]:.8e} Fz={Fp[2]:.8e}")
    print(f"Viscous : Fx={Fv[0]:.8e} Fy={Fv[1]:.8e} Fz={Fv[2]:.8e}")
    print(f"Total   : Fx={Ft[0]:.8e} Fy={Ft[1]:.8e} Fz={Ft[2]:.8e}")
    print()
    print(f"Cl      = {result['Cl']:.8e}")
    print(f"Cd      = {result['Cd']:.8e}")
    print(f"Cl,p    = {result['Cl_pressure']:.8e}")
    print(f"Cd,p    = {result['Cd_pressure']:.8e}")
    print(f"Cl,f    = {result['Cl_viscous']:.8e}")
    print(f"Cd,f    = {result['Cd_viscous']:.8e}")

    return result


# ======================================================================
# X AVERAGING PRESERVING UPPER / LOWER SURFACES
# ======================================================================

def average_by_x_side(x, z, values, weights, n_stations):
    xmin = np.min(x)
    xmax = np.max(x)

    bins = np.linspace(xmin, xmax, n_stations + 1)

    x_out = []
    z_out = []
    value_out = []

    for side in ("upper", "lower"):
        side_mask = z >= 0.0 if side == "upper" else z < 0.0

        for i in range(n_stations):
            if i == n_stations - 1:
                mask = (
                    side_mask
                    & (x >= bins[i])
                    & (x <= bins[i + 1])
                )
            else:
                mask = (
                    side_mask
                    & (x >= bins[i])
                    & (x < bins[i + 1])
                )

            if not np.any(mask):
                continue

            w = weights[mask]

            if np.sum(w) <= 0:
                continue

            x_mean = np.sum(x[mask] * w) / np.sum(w)
            z_mean = np.sum(z[mask] * w) / np.sum(w)
            value_mean = np.sum(values[mask] * w) / np.sum(w)

            x_out.append(x_mean)
            z_out.append(z_mean)
            value_out.append(value_mean)

    x_out = np.asarray(x_out)
    z_out = np.asarray(z_out)
    value_out = np.asarray(value_out)

    order = np.lexsort([
        z_out,
        x_out,
    ])

    return x_out[order], z_out[order], value_out[order]


# ======================================================================
# SURFACE COEFFICIENTS
# ======================================================================

def extract_surface_coefficients(p, tau, centres, normals, areas):
    x = centres[:, 0]
    y = centres[:, 1]
    z = centres[:, 2]

    # Incompressible Cp.
    Cp = (p - P_INF) / q_inf_kinematic()

    x_cp, z_cp, Cp_avg = average_by_x_side(
        x, z, Cp, areas, N_X_STATIONS
    )

    result = {
        "Cp_x": x_cp,
        "Cp_z": z_cp,
        "Cp": Cp_avg,
    }

    if tau is None:
        return result

    tau_n = np.sum(tau * normals, axis=1)
    tau_t = tau - tau_n[:, None] * normals
    tau_mag = np.linalg.norm(tau_t, axis=1)

    # tau is physical here, so use physical dynamic pressure.
    Cf = tau_mag / q_inf()

    x_cf, z_cf, Cf_avg = average_by_x_side(
        x, z, Cf, areas, N_X_STATIONS
    )

    result.update({
        "Cf_x": x_cf,
        "Cf_z": z_cf,
        "Cf": Cf_avg,
        "Cf_raw": Cf,
        "tau_raw": tau_t,
        "x_raw": x,
        "y_raw": y,
        "z_raw": z,
        "areas": areas,
        "tau_mag": tau_mag,
        "Cf_area_average": np.sum(Cf * areas) / np.sum(areas),
        "Cf_min": np.min(Cf),
        "Cf_max": np.max(Cf),
    })

    return result


# ======================================================================
# SAVE CP
# ======================================================================

def save_cp(cp, outdir):
    path = os.path.join(outdir, "Cp.csv")

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "x/c",
            "z/c",
            "Cp",
        ])

        for x, z, value in zip(
            cp["Cp_x"],
            cp["Cp_z"],
            cp["Cp"],
        ):
            writer.writerow([
                x / CHORD,
                z / CHORD,
                value,
            ])

    print(f"\nSaved: {path}")


# ======================================================================
# SAVE CF
# ======================================================================

def save_cf(cf, outdir):
    path = os.path.join(outdir, "Cf.csv")

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "x/c",
            "z/c",
            "Cf",
        ])

        for x, z, value in zip(
            cf["Cf_x"],
            cf["Cf_z"],
            cf["Cf"],
        ):
            writer.writerow([
                x / CHORD,
                z / CHORD,
                value,
            ])

    print(f"Saved: {path}")


# ======================================================================
# SAVE RAW 3D WALL DATA
# ======================================================================

def save_raw_wall(cf, outdir):
    path = os.path.join(outdir, "wallShearStress.csv")

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "x/c",
            "y/c",
            "z/c",
            "tau_w_Pa",
            "tau_x_Pa",
            "tau_y_Pa",
            "tau_z_Pa",
            "Cf",
            "surface_area_m2",
        ])

        order = np.lexsort([
            cf["z_raw"],
            cf["y_raw"],
            cf["x_raw"],
        ])

        for i in order:
            tau = cf["tau_raw"][i]

            writer.writerow([
                cf["x_raw"][i] / CHORD,
                cf["y_raw"][i] / CHORD,
                cf["z_raw"][i] / CHORD,
                np.linalg.norm(tau),
                tau[0],
                tau[1],
                tau[2],
                cf["Cf_raw"][i],
                cf["areas"][i],
            ])

    print(f"Saved: {path}")


# ======================================================================
# SAVE RAW SURFACE DATA
# ======================================================================

def save_raw_surface(p, centres, areas, tau, outdir):
    path = os.path.join(outdir, "raw_surface.csv")

    Cp = (p - P_INF) / q_inf_kinematic()

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "x/c",
            "y/c",
            "z/c",
            "p_kinematic_m2_s2",
            "p_physical_Pa",
            "Cp",
            "surface_area_m2",
        ])

        order = np.lexsort([
            centres[:, 2],
            centres[:, 1],
            centres[:, 0],
        ])

        for i in order:
            writer.writerow([
                centres[i, 0] / CHORD,
                centres[i, 1] / CHORD,
                centres[i, 2] / CHORD,
                p[i],
                RHO * p[i],
                Cp[i],
                areas[i],
            ])

    print(f"Saved: {path}")


# ======================================================================
# PLOT
# ======================================================================

def make_plot(cp, cf, outdir):
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    # Cp
    upper = cp["Cp_z"] >= 0.0
    lower = cp["Cp_z"] < 0.0

    upper_order = np.argsort(cp["Cp_x"][upper])
    lower_order = np.argsort(cp["Cp_x"][lower])

    axes[0].plot(
        cp["Cp_x"][upper][upper_order] / CHORD,
        cp["Cp"][upper][upper_order],
        label="Upper",
    )

    axes[0].plot(
        cp["Cp_x"][lower][lower_order] / CHORD,
        cp["Cp"][lower][lower_order],
        "--",
        label="Lower",
    )

    axes[0].invert_yaxis()
    axes[0].set(
        xlabel="x/c",
        ylabel="$C_p$",
        title="Spanwise-averaged $C_p$",
    )
    axes[0].grid(alpha=0.3)
    axes[0].legend()

    # Cf
    if cf is not None:
        upper = cf["Cf_z"] >= 0.0
        lower = cf["Cf_z"] < 0.0

        upper_order = np.argsort(cf["Cf_x"][upper])
        lower_order = np.argsort(cf["Cf_x"][lower])

        axes[1].plot(
            cf["Cf_x"][upper][upper_order] / CHORD,
            cf["Cf"][upper][upper_order],
            label="Upper",
        )

        axes[1].plot(
            cf["Cf_x"][lower][lower_order] / CHORD,
            cf["Cf"][lower][lower_order],
            "--",
            label="Lower",
        )

        axes[1].axhline(0.0, linewidth=0.8)
        axes[1].set(
            xlabel="x/c",
            ylabel="$C_f$",
            title="Spanwise-averaged $C_f$",
        )
        axes[1].grid(alpha=0.3)
        axes[1].legend()
    else:
        axes[1].set_visible(False)

    fig.tight_layout()

    path = os.path.join(outdir, "surface_coefficients.png")
    fig.savefig(path, dpi=200, bbox_inches="tight")
    plt.close(fig)

    print(f"Saved: {path}")


# ======================================================================
# SUMMARY
# ======================================================================

def save_summary(forces, cf, time, aoa, outdir):
    path = os.path.join(outdir, "summary.csv")

    with open(path, "w", newline="") as f:
        writer = csv.writer(f)

        writer.writerow([
            "time",
            "AoA_deg",
            "Cl",
            "Cd",
            "Cl_pressure",
            "Cd_pressure",
            "Cl_viscous",
            "Cd_viscous",
            "Fx",
            "Fy",
            "Fz",
            "Fx_pressure",
            "Fy_pressure",
            "Fz_pressure",
            "Fx_viscous",
            "Fy_viscous",
            "Fz_viscous",
            "surface_area",
            "span",
            "reference_area",
            "Cf_area_average",
            "Cf_min",
            "Cf_max",
        ])

        writer.writerow([
            time,
            aoa,
            forces["Cl"],
            forces["Cd"],
            forces["Cl_pressure"],
            forces["Cd_pressure"],
            forces["Cl_viscous"],
            forces["Cd_viscous"],
            forces["Fx"],
            forces["Fy"],
            forces["Fz"],
            forces["Fx_pressure"],
            forces["Fy_pressure"],
            forces["Fz_pressure"],
            forces["Fx_viscous"],
            forces["Fy_viscous"],
            forces["Fz_viscous"],
            forces["surface_area"],
            forces["span"],
            REFERENCE_AREA,
            cf["Cf_area_average"] if cf is not None else "",
            cf["Cf_min"] if cf is not None else "",
            cf["Cf_max"] if cf is not None else "",
        ])

    print(f"Saved: {path}")


# ======================================================================
# MAIN
# ======================================================================

def main():
    global AOA

    case_path = get_case_path()

    if AOA is None:
        AOA = infer_aoa(case_path)

    if AOA is None:
        raise RuntimeError(
            "Could not infer AoA from case name.\n"
            "Use a name such as naca0012Aoa10_ML or set AOA manually."
        )

    outdir = os.path.join(case_path, OUTPUT_DIR)
    os.makedirs(outdir, exist_ok=True)

    print("\n" + "=" * 70)
    print("NACA 0012 3D INCOMPRESSIBLE POST-PROCESSING")
    print("=" * 70)

    print(f"\nCase:\n    {case_path}")
    print(f"AoA:\n    {AOA} deg")

    print("\nIncompressible formulation:")
    print("    OpenFOAM p               = kinematic [m2/s2]")
    print("    Physical pressure        = rho*p")
    print("    OpenFOAM wallShearStress = kinematic [m2/s2]")
    print("    Physical wall shear      = rho*wallShearStress")

    print("\nFreestream:")
    print(f"    U_inf = {U_INF:.6f} m/s")
    print(f"    rho   = {RHO:.6f} kg/m^3")
    print(f"    p_inf = {P_INF:.6f} m2/s2")
    print(f"    q_inf = {q_inf():.6f} Pa")
    print(f"    q_inf/rho = {q_inf_kinematic():.6f} m2/s2")

    print("\nReference area:")
    print(f"    Aref = {REFERENCE_AREA:.6f} m2")

    # Read case.
    internal, boundary, time = read_case(case_path)

    # Airfoil.
    patch = get_airfoil_patch(boundary)

    # Geometry.
    surf, centres, Sf, normals, areas = get_surface_geometry(patch)

    # Fields.
    p = get_cell_pressure(surf)
    tau = get_wall_shear(surf)

    # Check sizes.
    n_faces = len(areas)

    print("\nField sizes:")
    print(f"    Surface faces = {n_faces}")
    print(f"    Pressure      = {len(p)}")

    if tau is not None:
        print(f"    Wall shear    = {len(tau)}")

    if len(p) != n_faces:
        raise RuntimeError(
            f"Pressure/surface mismatch: {len(p)} vs {n_faces}"
        )

    if tau is not None and len(tau) != n_faces:
        raise RuntimeError(
            f"Wall shear/surface mismatch: {len(tau)} vs {n_faces}"
        )

    # Forces.
    forces = compute_forces(
        p,
        tau,
        centres,
        Sf,
        normals,
        areas,
        AOA,
    )

    # Cp / Cf.
    coefficients = extract_surface_coefficients(
        p,
        tau,
        centres,
        normals,
        areas,
    )

    # Output.
    save_cp(coefficients, outdir)

    if tau is not None:
        save_cf(coefficients, outdir)
        save_raw_wall(coefficients, outdir)

    save_raw_surface(
        p,
        centres,
        areas,
        tau,
        outdir,
    )

    save_summary(
        forces,
        coefficients if tau is not None else None,
        time,
        AOA,
        outdir,
    )

    make_plot(
        coefficients,
        coefficients if tau is not None else None,
        outdir,
    )

    # Final results.
    print("\n" + "=" * 70)
    print("FINAL RESULTS")
    print("=" * 70)

    print(f"Cl       = {forces['Cl']:.8e}")
    print(f"Cd       = {forces['Cd']:.8e}")
    print(f"Cl,p     = {forces['Cl_pressure']:.8e}")
    print(f"Cd,p     = {forces['Cd_pressure']:.8e}")
    print(f"Cl,f     = {forces['Cl_viscous']:.8e}")
    print(f"Cd,f     = {forces['Cd_viscous']:.8e}")

    if tau is not None:
        print(f"Cf_avg   = {coefficients['Cf_area_average']:.8e}")

    print(f"Surface area = {forces['surface_area']:.8e} m2")
    print(f"Span         = {forces['span']:.8e} m")

    print(f"\nOutput:\n    {outdir}")
    print("\nDone.")


# ======================================================================
# ENTRY
# ======================================================================

if __name__ == "__main__":
    main()