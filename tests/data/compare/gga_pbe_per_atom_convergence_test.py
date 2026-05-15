"""Temporary test: per-atom GGA-PBE AE convergence (error vs Z), four panels.

Mimics ``atomSFE_with_boundary/scripts/visualize_gga_pbe_ae_fe_boundary_per_atom_convergence.py``,
but reads committed flat summaries under ``tests/data/summary/all_electron/gga_pbe/<sweep>/``
(``fe*_R*__*.json`` via ``summary_naming``).

Finest sweep case is the reference; each other case is one semilogy curve vs atomic number Z.

Run::

    python atomSFE/tests/data/compare/gga_pbe_per_atom_convergence_test.py
    python atomSFE/tests/data/compare/gga_pbe_per_atom_convergence_test.py --out path/to/figure.pdf
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np

_DATA_DIR = Path(__file__).resolve().parent.parent
if str(_DATA_DIR) not in sys.path:
    sys.path.insert(0, str(_DATA_DIR))
from summary_naming import glob_sweep_summaries, mesh_tag_from_summary_path

_COMPARE_DIR = Path(__file__).resolve().parent
_SUMMARY_DIR = _DATA_DIR / "summary"
_DEFAULT_GGA_PBE_ROOT = _SUMMARY_DIR / "all_electron" / "rscan"
_DEFAULT_OUT = _COMPARE_DIR / "gga_pbe_per_atom_convergence_test.pdf"

plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "stix",
        "font.size": 12,
    }
)


def _load_dataset_summary(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _parse_x_from_mesh_tag(mesh: str, mode: str) -> float:
    fe_txt, r_txt = mesh.split("_", 1)
    fe = float(fe_txt.replace("fe", ""))
    r = float(r_txt.replace("R", ""))
    return r if mode == "domain_radius_sweep" else fe


def _per_atom_metrics(payload: dict) -> dict[int, tuple[float, np.ndarray]]:
    out: dict[int, tuple[float, np.ndarray]] = {}
    for row in payload.get("config_summaries", []):
        z = row.get("atomic_number")
        if z is None:
            continue
        if not row.get("converged", True):
            continue
        e_tot = row.get("total_energy_ha")
        occ = row.get("occupied_eigenvalues_ha") or []
        if e_tot is None or len(occ) == 0:
            continue
        out[int(z)] = (float(e_tot), np.asarray(occ, dtype=float))
    return out


def _build_per_atom_curves(
    mode: str,
    gga_pbe_root: Path,
) -> tuple[list[tuple[float, np.ndarray, np.ndarray]], list[tuple[float, np.ndarray, np.ndarray]]]:
    """Return (energy_curves, eigen_curves); each item is (x_param, z_array, err_array)."""
    sweep_dir = gga_pbe_root / mode
    files = glob_sweep_summaries(sweep_dir)
    if not files:
        raise RuntimeError(f"No summary files found in {sweep_dir}")

    x_and_payload: list[tuple[float, dict]] = []
    for p in files:
        mesh = mesh_tag_from_summary_path(p)
        x = _parse_x_from_mesh_tag(mesh, mode)
        x_and_payload.append((x, _load_dataset_summary(p)))
    x_and_payload.sort(key=lambda t: t[0])

    x_ref, payload_ref = x_and_payload[-1]
    ref_map = _per_atom_metrics(payload_ref)

    energy_curves: list[tuple[float, np.ndarray, np.ndarray]] = []
    eigen_curves: list[tuple[float, np.ndarray, np.ndarray]] = []

    for x, payload in x_and_payload:
        if np.isclose(x, x_ref):
            continue
        cur_map = _per_atom_metrics(payload)
        shared = sorted(set(ref_map.keys()) & set(cur_map.keys()))
        if not shared:
            continue

        z_vals: list[int] = []
        e_errs: list[float] = []
        eig_errs: list[float] = []
        for z in shared:
            e_ref, eig_ref_all = ref_map[z]
            e_cur, eig_cur_all = cur_map[z]
            n = min(eig_ref_all.shape[0], eig_cur_all.shape[0])
            if n <= 0:
                continue
            z_vals.append(z)
            e_errs.append(abs(e_cur - e_ref))
            eig_errs.append(float(np.mean(np.abs(eig_cur_all[:n] - eig_ref_all[:n]))))

        if z_vals:
            z_arr = np.asarray(z_vals, dtype=float)
            energy_curves.append((float(x), z_arr, np.asarray(e_errs, dtype=float)))
            eigen_curves.append((float(x), z_arr, np.asarray(eig_errs, dtype=float)))

    return energy_curves, eigen_curves


def _plot_group(
    ax: plt.Axes,
    curves: list[tuple[float, np.ndarray, np.ndarray]],
    label_prefix: str,
    title: str,
) -> None:
    for x, z_arr, err_arr in curves:
        if label_prefix == "R":
            lbl = f"R={int(round(x))}"
        else:
            lbl = f"N_fe={int(round(x))}"
        ax.semilogy(z_arr, np.maximum(err_arr, 1e-20), marker="o", lw=1.3, ms=3.0, label=lbl)
    ax.set_title(title)
    ax.set_xlabel("Atomic number Z")
    ax.set_ylabel("Error (Ha)")
    ax.grid(True, which="major", axis="both", alpha=0.3, linestyle="-", linewidth=0.7)
    ax.grid(True, which="minor", axis="y", alpha=0.35, linestyle=":", linewidth=0.6)
    ax.minorticks_on()
    ax.legend(fontsize=8, ncol=2)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "Temporary: per-atom GGA-PBE AE convergence (4 panels: error vs Z per sweep case)."
        ),
    )
    ap.add_argument(
        "--gga-pbe-root",
        type=Path,
        default=_DEFAULT_GGA_PBE_ROOT,
        help="Path to gga_pbe summary root (domain_radius_sweep/, finite_element_sweep/).",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=_DEFAULT_OUT,
        help="Output PDF path (non-.pdf paths are forced to .pdf).",
    )
    args = ap.parse_args()

    root = args.gga_pbe_root.resolve()
    out_pdf = args.out.resolve()
    if out_pdf.suffix.lower() != ".pdf":
        out_pdf = out_pdf.with_suffix(".pdf")

    r_energy, r_eigen = _build_per_atom_curves("domain_radius_sweep", root)
    fe_energy, fe_eigen = _build_per_atom_curves("finite_element_sweep", root)

    fig, axes = plt.subplots(2, 2, figsize=(14, 9), constrained_layout=True)
    _plot_group(axes[0, 0], r_energy, "R", "Domain-radius sweep: Energy")
    _plot_group(axes[0, 1], r_eigen, "R", "Domain-radius sweep: Eigenvalues")
    _plot_group(axes[1, 0], fe_energy, "N_fe", "Finite-element sweep: Energy")
    _plot_group(axes[1, 1], fe_eigen, "N_fe", "Finite-element sweep: Eigenvalues")

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight", dpi=500)
    plt.close(fig)
    print(f"Wrote figure: {out_pdf}")


if __name__ == "__main__":
    main()