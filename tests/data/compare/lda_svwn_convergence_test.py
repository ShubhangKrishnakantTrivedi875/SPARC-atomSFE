"""LDA_SVWN sweep convergence test: max error vs finest reference (energy and mean occupied eigenvalue error).

Reads flat ``fe*_R*__*.json`` under ``summary/all_electron/lda_svwn/<sweep>/`` (from
``build_summary_from_out.py``). Reference for each sweep is the **finest** case (largest
:math:`R_{\\max}` for domain-radius sweep, largest :math:`N_{fe}` for FE sweep).

When aggregating max error over atomic numbers (main two-panel figure), **Z = 25 and Z = 51**
are omitted. ``--highlight-z`` plots that element alone (including 25 or 51 if requested).

Outputs PDF only (see ``--out``).

Run::

    python atomSFE/tests/data/compare/lda_svwn_convergence_test.py
    python atomSFE/tests/data/compare/lda_svwn_convergence_test.py --out path/to/figure.pdf
    python atomSFE/tests/data/compare/lda_svwn_convergence_test.py --highlight-z 26
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from matplotlib.ticker import LogFormatterMathtext, LogLocator, MultipleLocator

_DATA_DIR = Path(__file__).resolve().parent.parent
if str(_DATA_DIR) not in sys.path:
    sys.path.insert(0, str(_DATA_DIR))
from summary_naming import glob_sweep_summaries, mesh_tag_from_summary_path

_COMPARE_DIR = Path(__file__).resolve().parent
_SUMMARY_DIR = _DATA_DIR / "summary"
_DEFAULT_LDA_SVWN_ROOT = _SUMMARY_DIR / "all_electron" / "lda_svwn"
_DEFAULT_OUT_PDF = _COMPARE_DIR / "lda_svwn_convergence_test_summary.pdf"
_DEFAULT_OUT_HIGHLIGHT_PDF = _COMPARE_DIR / "lda_svwn_convergence_test_highlight.pdf"

plt.rcParams.update(
    {
        "font.family": "serif",
        "mathtext.fontset": "stix",
        "font.size": 15,
    }
)

AXIS_LABEL_FONTSIZE = 21
X_AXIS_LABEL_FONTSIZE = 24
TICK_LABEL_FONTSIZE = 15
LEGEND_FONTSIZE = 21

_COLOR_LINE_ENERGY = "C1"
_COLOR_LINE_EIGENVALUES = "C0"

_ATOMIC_NUMBERS_EXCLUDED_FROM_AGGREGATE_MAX: frozenset[int] = frozenset({25, 51})


def _load_dataset_summary(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def _parse_x_from_mesh_tag(mesh: str, mode: str) -> float:
    fe_txt, r_txt = mesh.split("_", 1)
    fe = float(fe_txt.replace("fe", ""))
    r = float(r_txt.replace("R", ""))
    return r if mode == "domain_radius_sweep" else fe


def _per_atom_metrics(payload: dict) -> dict[str, tuple[float, np.ndarray]]:
    out: dict[str, tuple[float, np.ndarray]] = {}
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
        out[str(int(z))] = (float(e_tot), np.asarray(occ, dtype=float))
    return out


def _build_curve(
    mode: str,
    lda_svwn_root: Path,
    atomic_number: int | None = None,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    sweep_dir = lda_svwn_root / mode
    files = glob_sweep_summaries(sweep_dir)
    if not files:
        raise RuntimeError(f"No summary files found in {sweep_dir}")

    x_and_payload: list[tuple[float, dict]] = []
    for p in files:
        x = _parse_x_from_mesh_tag(mesh_tag_from_summary_path(p), mode)
        x_and_payload.append((x, _load_dataset_summary(p)))
    x_and_payload.sort(key=lambda t: t[0])

    x_ref, payload_ref = x_and_payload[-1]
    ref_map = _per_atom_metrics(payload_ref)

    xs: list[float] = []
    y_energy: list[float] = []
    y_eigen: list[float] = []

    for x, payload in x_and_payload:
        if np.isclose(x, x_ref):
            continue
        cur_map = _per_atom_metrics(payload)
        shared = sorted(set(ref_map.keys()) & set(cur_map.keys()))
        if not shared:
            continue

        e_errs: list[float] = []
        eig_errs: list[float] = []
        for z in shared:
            zi = int(z)
            if atomic_number is None and zi in _ATOMIC_NUMBERS_EXCLUDED_FROM_AGGREGATE_MAX:
                continue
            if atomic_number is not None and zi != int(atomic_number):
                continue
            e_ref, eig_ref_all = ref_map[z]
            e_cur, eig_cur_all = cur_map[z]
            e_errs.append(abs(e_cur - e_ref))
            n = min(eig_ref_all.shape[0], eig_cur_all.shape[0])
            if n <= 0:
                continue
            per_atom_mean_eig_err = float(np.mean(np.abs(eig_cur_all[:n] - eig_ref_all[:n])))
            eig_errs.append(per_atom_mean_eig_err)

        if not e_errs or not eig_errs:
            continue
        xs.append(float(x))
        y_energy.append(float(np.max(np.asarray(e_errs, dtype=float))))
        y_eigen.append(float(np.max(np.asarray(eig_errs, dtype=float))))

    return np.asarray(xs), np.asarray(y_energy), np.asarray(y_eigen)


def _plot_convergence_panels(
    ax_l: plt.Axes,
    ax_r: plt.Axes,
    x_r: np.ndarray,
    y_r_energy: np.ndarray,
    y_r_eigen: np.ndarray,
    x_fe: np.ndarray,
    y_fe_energy: np.ndarray,
    y_fe_eigen: np.ndarray,
) -> None:
    ax_l.semilogy(
        x_r,
        np.maximum(y_r_energy, 1e-20),
        marker="o",
        lw=1.8,
        color=_COLOR_LINE_ENERGY,
        label="Energy",
    )
    ax_l.semilogy(
        x_r,
        np.maximum(y_r_eigen, 1e-20),
        marker="s",
        lw=1.8,
        color=_COLOR_LINE_EIGENVALUES,
        label="Eigenvalues",
    )
    ax_l.set_xlabel(r"$R_{max} \, \mathrm{(Bohr)}$", fontsize=X_AXIS_LABEL_FONTSIZE)
    ax_l.set_ylabel(r"Error (Ha)", fontsize=AXIS_LABEL_FONTSIZE)
    ax_l.yaxis.set_major_locator(LogLocator(base=10))
    ax_l.yaxis.set_minor_locator(LogLocator(base=10, subs=np.arange(2, 10)))
    ax_l.yaxis.set_major_formatter(LogFormatterMathtext())
    ax_l.minorticks_on()
    ax_l.tick_params(axis="both", which="major", labelsize=TICK_LABEL_FONTSIZE)
    ax_l.tick_params(axis="both", which="minor", labelsize=TICK_LABEL_FONTSIZE - 1)
    ax_l.grid(True, which="major", axis="both", alpha=0.32, linestyle="-", linewidth=0.75)
    ax_l.grid(True, which="minor", axis="y", alpha=0.40, linestyle=":", linewidth=0.70)
    ax_l.legend(fontsize=LEGEND_FONTSIZE)

    ax_r.semilogy(
        x_fe,
        np.maximum(y_fe_energy, 1e-20),
        marker="o",
        lw=1.8,
        color=_COLOR_LINE_ENERGY,
        label="Energy",
    )
    ax_r.semilogy(
        x_fe,
        np.maximum(y_fe_eigen, 1e-20),
        marker="s",
        lw=1.8,
        color=_COLOR_LINE_EIGENVALUES,
        label="Eigenvalues",
    )
    ax_r.set_xlabel(r"$N_{fe}$", fontsize=X_AXIS_LABEL_FONTSIZE)
    ax_r.set_ylabel(r"Error (Ha)", fontsize=AXIS_LABEL_FONTSIZE)
    ax_r.yaxis.set_major_locator(LogLocator(base=10, numticks=100))
    ax_r.yaxis.set_minor_locator(LogLocator(base=10, subs=tuple(np.arange(2, 10)), numticks=100))
    ax_r.yaxis.set_major_formatter(LogFormatterMathtext(base=10))
    ax_r.yaxis.set_minor_formatter(
        LogFormatterMathtext(base=10, labelOnlyBase=False, minor_thresholds=(np.inf, np.inf))
    )
    ax_r.minorticks_on()
    ax_r.tick_params(axis="both", which="major", labelsize=TICK_LABEL_FONTSIZE)
    ax_r.tick_params(axis="both", which="minor", labelsize=TICK_LABEL_FONTSIZE - 1)
    ax_r.tick_params(axis="y", which="minor", labelleft=False)
    if x_fe.size > 0:
        x_min = int(np.floor(float(np.min(x_fe))))
        x_max = int(np.ceil(float(np.max(x_fe))))
        ax_r.set_xlim(x_min - 0.5, x_max + 0.5)
        ax_r.xaxis.set_major_locator(MultipleLocator(1))
    ax_r.grid(True, which="major", axis="both", alpha=0.32, linestyle="-", linewidth=0.75)
    ax_r.grid(True, which="minor", axis="y", alpha=0.45, linestyle=":", linewidth=0.75)
    ax_r.legend(fontsize=LEGEND_FONTSIZE)


def main() -> None:
    ap = argparse.ArgumentParser(
        description=(
            "LDA_SVWN all-electron sweep convergence from summary JSON "
            "(max error vs finest case; Z=25 and Z=51 omitted from aggregate unless --highlight-z)."
        ),
    )
    ap.add_argument(
        "--lda-svwn-root",
        type=Path,
        default=_DEFAULT_LDA_SVWN_ROOT,
        help="Path to lda_svwn summary root (domain_radius_sweep/, finite_element_sweep/).",
    )
    ap.add_argument(
        "--out",
        type=Path,
        default=_DEFAULT_OUT_PDF,
        help="Output PDF path (aggregate: max error over Z, excluding 25 and 51).",
    )
    ap.add_argument(
        "--highlight-z",
        type=int,
        default=None,
        help="If set, write a second figure for this atomic number only (e.g. 26 for Fe).",
    )
    ap.add_argument(
        "--exclude-fe-x",
        type=float,
        nargs="*",
        default=(),
        help="Optional N_fe abscissa values to drop from the finite-element panel (e.g. 13.0).",
    )
    args = ap.parse_args()
    root = args.lda_svwn_root.resolve()
    out_pdf = args.out.resolve()
    if out_pdf.suffix.lower() != ".pdf":
        out_pdf = out_pdf.with_suffix(".pdf")
    exclude_fe = {float(v) for v in args.exclude_fe_x}

    def _curves_for(atomic_number: int | None) -> tuple[np.ndarray, ...]:
        x_r, y_r_energy, y_r_eigen = _build_curve("domain_radius_sweep", root, atomic_number=atomic_number)
        x_fe, y_fe_energy, y_fe_eigen = _build_curve("finite_element_sweep", root, atomic_number=atomic_number)
        if x_fe.size and exclude_fe:
            mask = np.ones(x_fe.shape[0], dtype=bool)
            for ex in exclude_fe:
                mask &= ~np.isclose(x_fe, ex)
            x_fe = x_fe[mask]
            y_fe_energy = y_fe_energy[mask]
            y_fe_eigen = y_fe_eigen[mask]
        fe_mask = x_fe > 1.0
        return (
            x_r,
            y_r_energy,
            y_r_eigen,
            x_fe[fe_mask],
            y_fe_energy[fe_mask],
            y_fe_eigen[fe_mask],
        )

    x_r, y_r_energy, y_r_eigen, x_fe, y_fe_energy, y_fe_eigen = _curves_for(None)

    fig, (ax_l, ax_r) = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
    _plot_convergence_panels(ax_l, ax_r, x_r, y_r_energy, y_r_eigen, x_fe, y_fe_energy, y_fe_eigen)

    out_pdf.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_pdf, format="pdf", bbox_inches="tight", dpi=600)
    plt.close(fig)
    print(f"Wrote figure: {out_pdf}")

    if args.highlight_z is not None:
        z = int(args.highlight_z)
        x_r_z, y_r_e_z, y_r_ev_z, x_fe_z, y_fe_e_z, y_fe_ev_z = _curves_for(z)
        if x_r_z.size == 0 and x_fe_z.size == 0:
            print(
                f"Warning: no data for Z={z} in summaries under {root}; skipped highlight figure.",
            )
        else:
            fig_z, (ax_l_z, ax_r_z) = plt.subplots(1, 2, figsize=(12, 4.8), constrained_layout=True)
            _plot_convergence_panels(
                ax_l_z, ax_r_z, x_r_z, y_r_e_z, y_r_ev_z, x_fe_z, y_fe_e_z, y_fe_ev_z
            )
            out_z = _DEFAULT_OUT_HIGHLIGHT_PDF
            if out_pdf != _DEFAULT_OUT_PDF:
                out_z = out_pdf.parent / f"{out_pdf.stem}_Z{z}.pdf"
            out_z = out_z.resolve()
            out_z.parent.mkdir(parents=True, exist_ok=True)
            fig_z.savefig(out_z, format="pdf", bbox_inches="tight", dpi=600)
            plt.close(fig_z)
            print(f"Wrote figure: {out_z}")


if __name__ == "__main__":
    main()
