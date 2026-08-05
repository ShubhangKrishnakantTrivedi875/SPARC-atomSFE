"""
Scaling report: phase-2 (N_fe / Lmax, Hg) and phase-3 (quadrature order q, Hg).

    python scaling_report.py

Writes to summary/tables/ :

    scaling_report.pdf   the presentable version
    scaling_report.md    the same tables as markdown
    scaling_fits.csv      fitted coefficients
    projections.csv       projected configurations
    time_ratio.csv         sc / nsc cost ratio, measured and projected

Pages, each self-contained (title + one plot OR one table, never both crowded onto
a page -- the previous version packed plot + table + footnote per page and text
collided):

    1  memory vs N_fe        (log-log line plot, both modes)
    2  time   vs N_fe        (log-log line plot, both modes)
    3  memory vs Lmax        (log-log line plot, both modes)
    4  time   vs Lmax        (log-log line plot, both modes)
    5  N_fe sweep table        (numbers behind pages 1-2)
    6  Lmax sweep table        (numbers behind pages 3-4)
    7  energy convergence table
    8  projections table  (N_fe=150/Lmax=80, N_fe=200/Lmax=100)
    9  sc / nsc cost ratio table
    10 memory vs q  (phase 3, Hg, full SCF)         log-log line plot
    11 time   vs q  (phase 3, Hg, full SCF)         log-log line plot
    12 phase-3 data table

Uses the project's categorical palette (blue = self-consistent RPA, orange =
RPA@DFT / non-self-consistent) -- slots 1-2, validated all-pairs in both color-vision
modes, so identity never rests on hue alone at just two series.

WHY THE MODELS LOOK LIKE THIS

Time.  The RPA is evaluated ONCE per outer cycle, and outer cycles do not depend on
grid size in a way that changes per-cycle cost, so:

    wall_full-SCF  =  outer_iters  x  channels  x  k n_quad^alpha

wall was regressed on channels at fixed n_quad first, as a check: the intercept is
consistent with zero at both n_fe=40 grid points tested, confirming that all
non-RPA cost (SCF outside the RPA build) is negligible next to the RPA itself.

n_quad = N_fe x q, but growing n_quad through N_fe also grows the FE problem
dimension, while growing it through q at fixed N_fe only densifies the
quadrature/projection matrices for an unchanged FE space.  These are NOT the same
operation and do not share one exponent -- phase 2 (grown via N_fe) and phase 3
(grown via q) are fitted and reported separately, never pooled into a single
n_quad^alpha.

Memory.  The n_quad x n_quad RPA kernels dominate (about 85% of the peak at the
phase-2 anchor), so for the N_fe sweep memory is fitted as A + B grid_points^2,
where grid_points = N_fe x p - 1 is the actual FE degrees of freedom -- the true
matrix-dimension variable, not n_quad = N_fe x q which mixes in the unrelated
quadrature order.  The Lmax sweep is fitted separately (affine in Lmax at fixed
grid points), since the two sweeps are orthogonal and a joint fit would be badly
conditioned on seven points.

The Lmax memory term is calibrated over Lmax 5-30 only, at a single grid-point
count, so how it behaves at Lmax=80-100 is genuinely uncertain.  Both bounding
assumptions are reported: `mem low` treats the Lmax term as additive and
grid-point-independent; `mem high` lets it scale with the grid-points-squared
term.  This is the largest uncertainty in the report and is stated as such
wherever a projection appears.

Energy convergence is reported as the difference between each row and the row
immediately before it in the same sweep -- no fitted or extrapolated reference.
The first row of a sweep has nothing to compare against, so its diff is blank.
"""

from __future__ import annotations

import json
import math
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from matplotlib.ticker import FuncFormatter, LogLocator, NullFormatter
import numpy as np
import pandas as pd

import cases as C

# phase3_quadrature also has a module named "cases" -- load it under a distinct
# name so it does not collide with (or shadow) this directory's own `cases` in
# sys.modules.
import importlib.util as _ilu
_p3_path = os.path.join(C.HERE, "..", "phase3_quadrature", "cases.py")
_spec = _ilu.spec_from_file_location("phase3_cases", _p3_path)
C3 = _ilu.module_from_spec(_spec)
_spec.loader.exec_module(C3)

OUT = os.path.join(C.HERE, "summary", "tables")
MIB = 1.048576
L_OCC = 3                       # Hg: 4f occupied, so channels = Lmax + L_OCC + 1
Q = C.FIXED["quadrature_point_number"]
P = C.FIXED["polynomial_order"]
L_REF = C.MESH_SWEEP[0][1]       # Lmax held fixed in the mesh sweep     (30)
N_REF = C.ANGULAR_SWEEP[0][0]    # N_fe held fixed in the angular sweep  (40)
GRID_REF = P * N_REF - 1         # grid points at the angular-sweep N_fe
PROJECTIONS = [(150, 80), (200, 100)]


def grid_points(n_fe):
    """
    FE degrees of freedom, N_fe elements of polynomial order p: p per element,
    minus 1 for the shared boundary condition.  This -- not n_quad = N_fe x q,
    which mixes in the unrelated quadrature order -- is what the RPA kernel
    matrices actually scale with when N_fe is the swept variable, so it is the
    fit variable for every N_fe-sweep model (pages 1-2).  The Lmax sweep and
    phase 3's q sweep are untouched by this: Lmax doesn't change n_fe, and
    phase 3 grows n_quad by varying q at fixed N_fe, which is a different
    operation with no grid_points analogue.
    """
    return P * n_fe - 1

MODE_LABEL = {"sc": "self-consistent RPA", "nsc": "non-self-consistent RPA@DFT"}
MODE_SHORT = {"sc": "scRPA", "nsc": "RPA@DFT"}

# palette.md slots 1 (blue) / 2 (orange) -- validated all-pairs, both CVD modes
COLOR = {"sc": "#2a78d6", "nsc": "#eb6834"}
INK        = "#0b0b0b"
INK_SOFT   = "#52514e"
MUTED      = "#898781"
GRID_LINE  = "#e1e0d9"
SURFACE    = "#fcfcfb"


def channels(l_max):
    return l_max + L_OCC + 1


def outer_count(outdir):
    """Outer-loop iterations, read from the log ('Outer iteration N' per cycle)."""
    p = os.path.join(outdir, "run.log")
    if not os.path.exists(p):
        return 0
    with open(p, errors="replace") as fh:
        return sum(1 for line in fh if line.startswith("Outer iteration"))


# ---------------------------------------------------------------- phase-2 data --

def load_phase2():
    rows = []
    for c in C.ALL_CASES:
        if c["mode"] != "scfull" and not (c["mode"] == "nsc" and c["code"] == "new"):
            continue
        p = os.path.join(c["outdir"], "result.json")
        if not os.path.exists(p):
            continue
        r = json.load(open(p))
        if not r.get("ok"):
            continue
        peak = (r.get("memory") or {}).get("peak_rss_tree_mb")
        mode = "sc" if c["mode"] == "scfull" else "nsc"
        rows.append(dict(
            mode=mode, n_fe=c["n_fe"], l_max=c["l_max"],
            n_quad=c["n_fe"] * Q, chan=channels(c["l_max"]),
            mem=(peak / MIB) if peak else np.nan,
            wall=r.get("wall_s"),
            outer=outer_count(c["outdir"]) or 1,
            energy=r.get("total_energy")))
    return pd.DataFrame(rows)


def load_phase3():
    rows = []
    for c in C3.ALL_CASES:
        if c["variant"] != "densep":
            continue
        p = os.path.join(c["outdir"], "result.json")
        if not os.path.exists(p):
            continue
        r = json.load(open(p))
        if not r.get("ok"):
            continue
        peak = (r.get("memory") or {}).get("peak_rss_tree_mb")
        rows.append(dict(
            mode=c["mode"], q=c["q"],
            n_quad=C3.FIXED["finite_element_number"] * c["q"],
            mem=(peak / MIB) if peak else np.nan,
            wall=r.get("wall_s"),
            outer=outer_count(c["outdir"]) or 1,
            energy=r.get("total_energy")))
    return pd.DataFrame(rows)


# ---------------------------------------------------------------- fits ----------

def fit_power(x, y):
    A = np.column_stack([np.ones(len(x)), np.log(x)])
    coef, *_ = np.linalg.lstsq(A, np.log(y), rcond=None)
    k, a = math.exp(coef[0]), coef[1]
    return k, a, float(np.max(np.abs(k * np.asarray(x) ** a - y) / y))


def fit_affine(x, y):
    A = np.column_stack([np.ones(len(x)), np.asarray(x, float)])
    coef, *_ = np.linalg.lstsq(A, np.asarray(y, float), rcond=None)
    c0, s = coef
    return c0, s, float(np.max(np.abs(c0 + s * np.asarray(x) - y) / y))


def fit_mode(d):
    mesh = d[d.l_max == L_REF].sort_values("n_fe")
    ang  = d[d.n_fe == N_REF].sort_values("l_max")
    gpts = grid_points(mesh.n_fe.values)

    mem_a, mem_b, mem_e = fit_affine(gpts ** 2.0, mesh.mem.values)
    mL_a, mL_b, mL_e    = fit_affine(ang.l_max.values, ang.mem.values)

    # per-(channel, outer-cycle) unit cost, as a function of grid points
    unit = mesh.wall.values / (mesh.chan.values * mesh.outer.values)
    t_k, t_alpha, t_e = fit_power(gpts, unit)

    anchor_mem = float(d[(d.n_fe == N_REF) & (d.l_max == L_REF)].mem.iloc[0])
    quad_share = mem_b * GRID_REF ** 2 / anchor_mem
    anchor_outer = int(d[(d.n_fe == N_REF) & (d.l_max == L_REF)].outer.iloc[0])

    return dict(mem_base=mem_a, mem_quad=mem_b, mem_err=mem_e,
                memL_base=mL_a, memL_slope=mL_b, memL_err=mL_e,
                quad_share=quad_share, outer_at_anchor=anchor_outer,
                t_k=t_k, t_alpha=t_alpha, t_err=t_e)


def project_memory(f, n_fe, l_max):
    gp = grid_points(n_fe)
    quad = f["mem_quad"] * gp ** 2
    low  = f["mem_base"] + quad + f["memL_slope"] * (l_max - L_REF)
    gamma = f["memL_slope"] / (f["mem_quad"] * GRID_REF ** 2)
    high = f["mem_base"] + quad * (1.0 + gamma * (l_max - L_REF))
    return min(low, high), max(low, high)


def project_time(f, n_fe, l_max, outer):
    return outer * channels(l_max) * f["t_k"] * grid_points(n_fe) ** f["t_alpha"]


def fit_phase3(d):
    """One power law per mode, wall/mem vs n_quad, both grown via q at fixed N_fe."""
    out = {}
    for m in ("sc", "nsc"):
        g = d[d["mode"] == m].sort_values("q")
        if g.empty:
            continue
        tk, ta, te = fit_power(g.n_quad.values, g.wall.values)
        mk, ma, me = fit_power(g.n_quad.values, g.mem.values)
        out[m] = dict(t_k=tk, t_alpha=ta, t_err=te, m_k=mk, m_alpha=ma, m_err=me,
                      outer=int(g.outer.iloc[0]), spread=float(g.energy.max() - g.energy.min()))
    return out


# ---------------------------------------------------------------- plot helpers --

def style_axes(ax, xlabel, ylabel, title, xticks):
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel(xlabel, fontsize=9.5, color=INK_SOFT)
    ax.set_ylabel(ylabel, fontsize=9.5, color=INK_SOFT)
    ax.set_title(title, fontsize=12.5, fontweight="bold", color=INK, loc="left", pad=10)
    ax.set_facecolor(SURFACE)
    ax.grid(True, which="major", color=GRID_LINE, linewidth=0.8, zorder=0)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    for spine in ("left", "bottom"):
        ax.spines[spine].set_color(MUTED)
    ax.tick_params(colors=INK_SOFT, labelsize=8.5)
    # x-ticks at exactly the sampled values, not matplotlib's default log decades --
    # the sweep only ever has 4-6 points, so decade ticks either miss them all or
    # place just one or two, neither of which lets a reader read a value off the axis
    ax.set_xticks(xticks)
    ax.set_xticks([], minor=True)
    ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: "%g" % v))
    ax.xaxis.set_minor_formatter(NullFormatter())


def _try_yticks(all_y, min_ratio=1.18):
    """
    Sampled y-values as tick locations, IF they clear a minimum log-spacing --
    below that, adjacent tick labels overlap on a compressed log axis.  Returns
    the tick list, or None to signal "fall back to per-point value labels".
    """
    ys = sorted(set(float(v) for v in all_y))
    for a, b in zip(ys, ys[1:]):
        if b / a < min_ratio:
            return None
    return ys


def line_plot(ax, series, xlabel, ylabel, title, note, value_fmt="%.0f"):
    """
    series = [(mode, x, y, fit_curve, fit_label), ...]

    fit_curve is a callable x -> y (the fitted scaling law) or None.  Measured
    points are markers only; the fit is a separate dashed curve in the same color,
    so "what was measured" and "what the model predicts" stay visually distinct
    even where they nearly overlap.

    Y-axis ticks are placed at the sampled values when they clear a minimum
    spacing; when two series' points sit too close together on the log scale for
    that to read cleanly, each point gets its own small value label instead (never
    both -- that doubles up the same number).
    """
    xticks = sorted(set(float(v) for s in series for v in s[1]))
    style_axes(ax, xlabel, ylabel, title, xticks)

    all_y = [v for s in series for v in s[2]]
    yticks = _try_yticks(all_y)
    if yticks is not None:
        ax.set_yticks(yticks)
        ax.set_yticks([], minor=True)
        # value_fmt, not a generic "%g" -- the data are floats with long tails
        # (e.g. 584.0439999...), and "%g" prints those tails as ugly 6-sig-fig
        # labels like "584.044" instead of the clean "584.0" every other number
        # in this report uses
        ax.yaxis.set_major_formatter(FuncFormatter(lambda v, _: value_fmt % v))
        ax.yaxis.set_minor_formatter(NullFormatter())
    else:
        ax.yaxis.set_major_formatter(NullFormatter())
        ax.yaxis.set_minor_formatter(NullFormatter())
        ax.tick_params(axis="y", length=0)

    for mode, x, y, fit_curve, fit_label in series:
        c = COLOR[mode]
        if fit_curve is not None:
            xs = np.geomspace(min(x), max(x), 200)
            ax.plot(xs, fit_curve(xs), color=c, linewidth=1.6, linestyle=(0, (5, 3)),
                    alpha=0.55, zorder=2)
        ax.plot(x, y, color=c, linewidth=0, marker="o", markersize=7,
                markerfacecolor=c, markeredgecolor=SURFACE, markeredgewidth=1.2,
                zorder=3, label=MODE_LABEL[mode])
        if yticks is None:
            # no y-ticks to read a value off of, so every point states its own,
            # alternating above/below so consecutive labels don't collide
            for i, (xi, yi) in enumerate(zip(x, y)):
                if xi == x[-1]:
                    continue  # last point already carries the fit-formula label
                dy = 9 if i % 2 == 0 else -12
                va = "bottom" if i % 2 == 0 else "top"
                ax.annotate(value_fmt % yi, xy=(xi, yi), xytext=(0, dy),
                            textcoords="offset points", fontsize=7.3, color=c,
                            ha="center", va=va)
        # direct label at the last point, anchored so it grows LEFT into the plot --
        # anchoring right and growing right runs the text off the figure edge, which
        # is exactly the overlap the first version of this report shipped with
        ax.annotate(fit_label, xy=(x[-1], y[-1]), xytext=(-8, 4),
                    textcoords="offset points", fontsize=8.2, color=c,
                    fontweight="bold", va="bottom", ha="right")
    handles = [plt.Line2D([], [], color=COLOR[m], marker="o", linestyle="none",
                          markersize=7, label=MODE_LABEL[m]) for m, *_ in
              {s[0]: s for s in series}.values()]
    dashed = plt.Line2D([], [], color=MUTED, linestyle=(0, (5, 3)), linewidth=1.6,
                        label="fitted scaling law")
    ax.legend(handles=handles + [dashed], loc="upper left", frameon=False,
              fontsize=9, labelcolor=INK_SOFT, handlelength=1.8, borderaxespad=0.3)
    if note:
        ax.text(0.0, -0.16, note, transform=ax.transAxes, fontsize=8.0,
                color=INK_SOFT, va="top", wrap=True)


def full_page(pdf, draw):
    fig = plt.figure(figsize=(9.5, 7.2))
    fig.patch.set_facecolor(SURFACE)
    # right margin is deliberately wide: the direct labels anchor at the last point
    # and grow left, but still need clearance so descenders/ascenders don't clip
    ax = fig.add_axes([0.12, 0.20, 0.78, 0.68])
    draw(ax, fig)
    pdf.savefig(fig, facecolor=SURFACE)
    plt.close(fig)


# ---------------------------------------------------------------- tables --------

INT_KEYS = ("N_fe", "Lmax", "n_quad", "grid points", "channels", "q", "outer")


def fmt(k, v):
    if isinstance(v, str):
        return v
    if v is None or (isinstance(v, float) and (pd.isna(v) or not np.isfinite(v))):
        return "-"
    if k in INT_KEYS:
        return "%d" % v
    if "ratio" in k:
        return "%.2fx" % v
    if "GiB" in k:
        return "%.1f" % v
    if "(h)" in k:
        return "%.1f" % v
    if "|E" in k:
        return "%.2e" % v
    if k == "dE vs previous (Ha)":
        return "%+.2e" % v
    if k == "E (Ha)":
        return "%.9f" % v
    if "MiB" in k or "(s)" in k:
        return "%.1f" % v
    return "%.3g" % v


def to_md(d):
    keys = list(d.columns)
    lines = ["| " + " | ".join(keys) + " |",
             "|" + "|".join(["---:"] * len(keys)) + "|"]
    for _, r in d.iterrows():
        lines.append("| " + " | ".join(fmt(k, r[k]) for k in keys) + " |")
    return "\n".join(lines)


def table_page(pdf, title, sub, d, highlight=(), extra_note=None):
    """A table gets its own page -- no plot sharing it, so nothing overlaps."""
    n_rows = len(d) + 1
    fig_h = min(10.5, max(3.4, 1.05 + 0.42 * n_rows + (0.55 if extra_note else 0)))
    fig = plt.figure(figsize=(11.0, fig_h))
    fig.patch.set_facecolor(SURFACE)
    fig.suptitle(title, fontsize=13.5, fontweight="bold", x=0.05, ha="left",
                y=1.0 - 0.35 / fig_h, color=INK)
    if sub:
        fig.text(0.05, 1.0 - 0.75 / fig_h, sub, fontsize=8.6, va="top",
                 color=INK_SOFT, wrap=True)
    top_frac = 1.0 - (1.35 if sub else 0.75) / fig_h
    bottom_frac = (0.65 if extra_note else 0.12) / fig_h
    ax = fig.add_axes([0.035, bottom_frac, 0.93, top_frac - bottom_frac])
    ax.axis("off")
    keys = list(d.columns)
    cells = [[fmt(k, r[k]) for k in keys] for _, r in d.iterrows()]
    t = ax.table(cellText=cells, colLabels=keys, cellLoc="right", loc="upper center")
    t.auto_set_font_size(False)
    t.set_fontsize(9.0)
    t.scale(1, 1.9)
    for j in range(len(keys)):
        c = t[0, j]
        c.set_facecolor("#e4e4e4")
        c.set_text_props(fontweight="bold", color=INK)
    for i in range(1, len(cells) + 1):
        for j, k in enumerate(keys):
            if i % 2 == 0:
                t[i, j].set_facecolor("#f7f7f7")
            t[i, j].set_text_props(color=INK)
            if k in highlight:
                t[i, j].set_text_props(fontweight="bold", color=INK)
    if "kind" in keys:
        j = keys.index("kind")
        for i in range(1, len(cells) + 1):
            if cells[i - 1][j] == "projected":
                for jj in range(len(keys)):
                    t[i, jj].set_facecolor("#fff1d6")
    for _, cell in t.get_celld().items():
        cell.set_linewidth(0.4)
        cell.set_edgecolor("#bbbbbb")
    if extra_note:
        fig.text(0.05, bottom_frac - 0.02, extra_note, fontsize=7.6,
                 va="top", color=INK_SOFT, wrap=True)
    pdf.savefig(fig, facecolor=SURFACE)
    plt.close(fig)


# ---------------------------------------------------------------- table builders

def mesh_table(d, mode):
    g = d[(d["mode"] == mode) & (d.l_max == L_REF)].sort_values("n_fe")
    out, prev = [], None
    for _, r in g.iterrows():
        gp = grid_points(r.n_fe)
        row = {"N_fe": int(r.n_fe), "grid points": int(gp), "outer": int(r.outer),
               "mem (MiB)": r.mem, "wall (s)": r.wall,
               "s/channel/outer": r.wall / (channels(r.l_max) * r.outer)}
        if prev is not None:
            b = math.log(gp / grid_points(prev.n_fe))
            row["mem exp"]  = math.log(r.mem / prev.mem) / b
            row["time exp"] = math.log(r.wall / prev.wall) / b
        out.append(row)
        prev = r
    return pd.DataFrame(out)


def ang_table(d, mode):
    g = d[(d["mode"] == mode) & (d.n_fe == N_REF)].sort_values("l_max")
    out, prev = [], None
    for _, r in g.iterrows():
        row = {"Lmax": int(r.l_max), "channels": int(r.chan), "outer": int(r.outer),
               "mem (MiB)": r.mem, "wall (s)": r.wall,
               "s/channel/outer": r.wall / (r.chan * r.outer)}
        if prev is not None:
            b = math.log(r.chan / prev.chan)
            row["mem exp"]  = math.log(r.mem / prev.mem) / b
            row["time exp"] = math.log(r.wall / prev.wall) / b
        out.append(row)
        prev = r
    return pd.DataFrame(out)


def sequential_diff_table(d, mode, axis):
    """
    Energy against its own immediately-preceding row -- no fitted/extrapolated
    reference.  The first row has nothing before it, so its diff is blank rather
    than compared against anything invented.
    """
    fixed_val = L_REF if axis == "n_fe" else N_REF
    fixed_col = "l_max" if axis == "n_fe" else "n_fe"
    g = d[(d["mode"] == mode) & (d[fixed_col] == fixed_val)].sort_values(axis)
    out, prev = [], None
    for _, r in g.iterrows():
        row = {"N_fe": int(r.n_fe), "Lmax": int(r.l_max), "E (Ha)": r.energy,
               "dE vs previous (Ha)": None if prev is None else r.energy - prev}
        out.append(row)
        prev = r.energy
    return pd.DataFrame(out)


def projection_table(fits):
    out = []
    for mode in ("sc", "nsc"):
        if mode not in fits:
            continue
        f = fits[mode]
        outer = f["outer_at_anchor"] if mode == "sc" else 1
        for n_fe, l_max in PROJECTIONS:
            lo, hi = project_memory(f, n_fe, l_max)
            t = project_time(f, n_fe, l_max, outer)
            out.append({"mode": MODE_SHORT[mode], "N_fe": n_fe, "Lmax": l_max,
                        "grid points": grid_points(n_fe), "channels": channels(l_max),
                        "outer": outer,
                        "mem low (GiB)": lo / 1024.0, "mem high (GiB)": hi / 1024.0,
                        "wall (h)": t / 3600.0})
    return pd.DataFrame(out)


def ratio_table(df, fits):
    out = []
    for (n_fe, l_max), g in df.groupby(["n_fe", "l_max"], sort=True):
        a, b = g[g["mode"] == "sc"], g[g["mode"] == "nsc"]
        if a.empty or b.empty:
            continue
        a, b = a.iloc[0], b.iloc[0]
        out.append({"N_fe": int(n_fe), "Lmax": int(l_max), "grid points": grid_points(n_fe),
                    "scRPA (s)": a.wall, "RPA@DFT (s)": b.wall, "ratio": a.wall / b.wall,
                    "scRPA MiB": a.mem, "RPA@DFT MiB": b.mem,
                    "mem ratio": a.mem / b.mem, "kind": "measured"})
    if "sc" in fits and "nsc" in fits:
        for n_fe, l_max in PROJECTIONS:
            ts = project_time(fits["sc"], n_fe, l_max, fits["sc"]["outer_at_anchor"])
            tn = project_time(fits["nsc"], n_fe, l_max, 1)
            ms = project_memory(fits["sc"], n_fe, l_max)
            mn = project_memory(fits["nsc"], n_fe, l_max)
            out.append({"N_fe": n_fe, "Lmax": l_max, "grid points": grid_points(n_fe),
                        "scRPA (s)": ts, "RPA@DFT (s)": tn, "ratio": ts / tn,
                        "scRPA MiB": ms[0], "RPA@DFT MiB": mn[0],
                        "mem ratio": ms[0] / mn[0], "kind": "projected"})
    return pd.DataFrame(out)


def phase3_table(d):
    out = []
    for _, r in d.sort_values(["mode", "q"]).iterrows():
        out.append({"mode": MODE_SHORT[r["mode"]], "q": int(r.q), "n_quad": int(r.n_quad),
                    "outer": int(r.outer), "mem (MiB)": r.mem, "wall (s)": r.wall,
                    "E (Ha)": r.energy})
    return pd.DataFrame(out)


# ---------------------------------------------------------------- main ----------

def main():
    os.makedirs(OUT, exist_ok=True)
    df = load_phase2()
    d3 = load_phase3()
    if df.empty:
        print("no completed phase-2 full-SCF / RPA@DFT cases")
        return

    fits = {m: fit_mode(df[df["mode"] == m]) for m in ("sc", "nsc")
            if not df[df["mode"] == m].empty}
    proj  = projection_table(fits)
    ratio = ratio_table(df, fits)
    fits3 = fit_phase3(d3) if not d3.empty else {}

    HEADER = ("Mercury (Z=80), all-electron, corrected SPARC-atomSFE, full SCF "
              "(outer loop converged).  Fixed: p=20, q=%d, domain=%.0f Bohr, omega=%d, "
              "polynomial mesh (concentration 2).  Channels = Lmax + 4 "
              "(Hg has 4f occupied, l_occ_max = 3)."
              % (Q, C.FIXED["domain_size"], C.FIXED["frequency_quadrature_point_number"]))

    md = ["# Scaling report -- phase 2 (N_fe, Lmax) and phase 3 (quadrature order)",
          "", HEADER, ""]

    with PdfPages(os.path.join(OUT, "scaling_report.pdf")) as pdf:

        # ---- cover ----
        fig = plt.figure(figsize=(9.5, 7.2))
        fig.patch.set_facecolor(SURFACE)
        fig.text(0.08, 0.80, "RPA scaling report", fontsize=22, fontweight="bold", color=INK)
        fig.text(0.08, 0.72, "Phase 2: N_fe / Lmax sweeps  --  Phase 3: quadrature order q",
                 fontsize=12, color=INK_SOFT)
        fig.text(0.08, 0.60, HEADER, fontsize=9.5, color=INK_SOFT, wrap=True)
        contents = ["1  Memory vs N_fe", "2  Wall time vs N_fe", "3  Memory vs Lmax",
                    "4  Wall time vs Lmax", "5  N_fe sweep, numeric",
                    "6  Lmax sweep, numeric",
                    "7  Total-energy convergence (successive differences, per sweep)",
                    "8  Projections (N_fe=150/Lmax=80, N_fe=200/Lmax=100)",
                    "9  Self-consistent / non-self-consistent cost ratio",
                    "10 Phase 3: memory vs quadrature order q  (Hg)",
                    "11 Phase 3: wall time vs quadrature order q  (Hg)",
                    "12 Phase 3: numeric data"]
        fig.text(0.08, 0.50, "Contents", fontsize=11, fontweight="bold", color=INK)
        for i, line in enumerate(contents):
            fig.text(0.09, 0.46 - 0.032 * i, line, fontsize=9.3, color=INK_SOFT)
        pdf.savefig(fig, facecolor=SURFACE)
        plt.close(fig)

        # ---- pages 1-2: N_fe ----
        mesh_n_fe = sorted({n for n, l in C.MESH_SWEEP})
        gpts_map = ", ".join("%d→%d" % (n, grid_points(n)) for n in mesh_n_fe)
        for ykey, ylabel, page_title, fmt in (
                ("mem", "peak memory (MiB)", "Memory vs N_fe", "%.0f"),
                ("wall", "wall time (s)", "Wall time vs N_fe", "%.1f")):
            series = []
            for m in fits:
                g = df[(df["mode"] == m) & (df.l_max == L_REF)].sort_values("n_fe")
                f = fits[m]
                if ykey == "mem":
                    lbl = "%s\nMiB ≈ %.0f + %.2e·(grid pts)²" % (
                        MODE_SHORT[m], f["mem_base"], f["mem_quad"])
                    curve = (lambda nfe, f=f:
                             f["mem_base"] + f["mem_quad"] * grid_points(nfe) ** 2)
                else:
                    outer = f["outer_at_anchor"] if m == "sc" else 1
                    lbl = "%s\nwall ∝ (grid pts)^%.2f  (%d outer cycle%s)" % (
                        MODE_SHORT[m], f["t_alpha"], outer, "" if outer == 1 else "s")
                    curve = (lambda nfe, f=f, outer=outer: outer * channels(L_REF) *
                             f["t_k"] * grid_points(nfe) ** f["t_alpha"])
                series.append((m, g.n_fe.values, g[ykey].values, curve, lbl))
            note = ("Lmax held at %d.  x-axis is N_fe; the fit is against grid points "
                    "(= N_fe x p - 1, the actual FE degrees of freedom -- p=%d here), "
                    "since that is what the RPA kernel matrices scale with, not the "
                    "quadrature order.  N_fe -> grid points: %s."
                    % (L_REF, P, gpts_map))
            full_page(pdf, lambda ax, fig, s=series, yl=ylabel, t=page_title, n=note, f=fmt:
                      line_plot(ax, s, "N_fe", yl, t, n, value_fmt=f))
            md += ["## %s" % page_title, "", note, ""]

        # ---- pages 3-4: Lmax ----
        for ykey, ylabel, page_title, fmt in (
                ("mem", "peak memory (MiB)", "Memory vs Lmax", "%.0f"),
                ("wall", "wall time (s)", "Wall time vs Lmax", "%.1f")):
            series = []
            for m in fits:
                g = df[(df["mode"] == m) & (df.n_fe == N_REF)].sort_values("l_max")
                f = fits[m]
                if ykey == "mem":
                    lbl = "%s\nMiB ≈ %.0f + %.2f·Lmax" % (
                        MODE_SHORT[m], f["memL_base"], f["memL_slope"])
                    curve = (lambda lm, f=f: f["memL_base"] + f["memL_slope"] * lm)
                else:
                    outer = f["outer_at_anchor"] if m == "sc" else 1
                    lbl = "%s\n%d outer cycle%s, linear in channels" % (
                        MODE_SHORT[m], outer, "" if outer == 1 else "s")
                    curve = (lambda lm, f=f, outer=outer:
                             outer * channels(lm) * f["t_k"] * GRID_REF ** f["t_alpha"])
                series.append((m, g.l_max.values, g[ykey].values, curve, lbl))
            note = ("N_fe held at %d (grid points = %d), so only the channel count "
                    "(Lmax + 4) moves.  x-axis is Lmax."
                    % (N_REF, GRID_REF))
            full_page(pdf, lambda ax, fig, s=series, yl=ylabel, t=page_title, n=note, f=fmt:
                      line_plot(ax, s, "Lmax", yl, t, n, value_fmt=f))
            md += ["## %s" % page_title, "", note, ""]

        # ---- page 5-6: numeric tables behind the plots ----
        for mode in fits:
            table_page(pdf, "N_fe sweep, numeric -- %s" % MODE_LABEL[mode], HEADER,
                      mesh_table(df, mode), highlight=("mem exp", "time exp"))
            md += ["### N_fe sweep -- %s" % MODE_LABEL[mode], "",
                   to_md(mesh_table(df, mode)), ""]
        for mode in fits:
            table_page(pdf, "Lmax sweep, numeric -- %s" % MODE_LABEL[mode], HEADER,
                      ang_table(df, mode), highlight=("mem exp", "time exp"))
            md += ["### Lmax sweep -- %s" % MODE_LABEL[mode], "",
                   to_md(ang_table(df, mode)), ""]

        # ---- page 7: energy convergence ----
        # Each row is compared to the row immediately before it in its own sweep --
        # no fitted or extrapolated reference.  The first row of each sweep has
        # nothing to compare against, so its diff is left blank rather than filled
        # in against an invented target.
        for mode in fits:
            n_fe_diff = sequential_diff_table(df, mode, "n_fe")
            l_max_diff = sequential_diff_table(df, mode, "l_max")
            note = ("`dE vs previous` is this row's energy minus the row above it, "
                    "within the same sweep (N_fe sweep at Lmax=%d; Lmax sweep at "
                    "N_fe=%d) -- not a comparison against any predicted or "
                    "extrapolated converged value." % (L_REF, N_REF))
            table_page(pdf, "Total-energy convergence, N_fe sweep -- %s" % MODE_LABEL[mode],
                      HEADER, n_fe_diff, highlight=("dE vs previous (Ha)",),
                      extra_note=note)
            table_page(pdf, "Total-energy convergence, Lmax sweep -- %s" % MODE_LABEL[mode],
                      HEADER, l_max_diff, highlight=("dE vs previous (Ha)",),
                      extra_note=note)
            md += ["## Energy convergence, N_fe sweep -- %s" % MODE_LABEL[mode], "",
                   note, "", to_md(n_fe_diff), "",
                   "## Energy convergence, Lmax sweep -- %s" % MODE_LABEL[mode], "",
                   note, "", to_md(l_max_diff), ""]

        # ---- page 8: projections ----
        note = ("`mem low` treats the Lmax term as additive and grid-point-independent; "
                "`mem high` lets it scale with the grid-points-squared term.  Only one "
                "grid-point count was sampled in the angular sweep, so which is correct "
                "is undetermined -- the true value lies between them.  This is the "
                "largest uncertainty in this report.  Wall time needs no such caveat: "
                "it is the measured outer-cycle count (scRPA: %d, from the "
                "N_fe=%d/Lmax=%d anchor; RPA@DFT: 1, by construction) times the fitted "
                "per-channel cost."
                % (fits.get("sc", {}).get("outer_at_anchor", 0), N_REF, L_REF))
        table_page(pdf, "Projections", HEADER, proj,
                  highlight=("mem low (GiB)", "mem high (GiB)", "wall (h)"),
                  extra_note=note)
        md += ["## Projections", "", note, "", to_md(proj), ""]

        # ---- page 9: cost ratio ----
        note = ""
        if "sc" in fits and "nsc" in fits:
            note = ("Ratio is scRPA / RPA@DFT wall time, using each mode's own measured "
                    "outer-cycle count -- no per-iteration normalisation is applied since "
                    "neither number was ever per-iteration.  The ratio grows with grid size "
                    "because scRPA's per-channel cost has the steeper exponent in grid "
                    "points (%.2f vs %.2f).  Shaded rows are projected."
                    % (fits["sc"]["t_alpha"], fits["nsc"]["t_alpha"]))
        table_page(pdf, "Self-consistent / non-self-consistent cost ratio", HEADER,
                  ratio, highlight=("ratio", "mem ratio"), extra_note=note)
        md += ["## Cost ratio", "", note, "", to_md(ratio), ""]

        # ---- pages 10-11: phase 3 (Hg, quadrature order q) ----
        if not d3.empty:
            H3 = ("Mercury (Z=80), all-electron, full SCF, N_fe=%d, p=%d, Lmax=%d, "
                  "omega=%d, dense basis order = p.  n_quad = N_fe x q grows through q "
                  "at FIXED N_fe here -- a different operation from phase 2's N_fe sweep, "
                  "so the exponents below are not comparable to pages 1-2."
                  % (C3.FIXED["finite_element_number"], C3.FIXED["polynomial_order"],
                     C3.FIXED["angular_momentum_cutoff"],
                     C3.FIXED["frequency_quadrature_point_number"]))
            for ykey, ylabel, page_title, fmt in (
                    ("mem", "peak memory (MiB)", "Phase 3 -- memory vs quadrature order q", "%.0f"),
                    ("wall", "wall time (s)", "Phase 3 -- wall time vs quadrature order q", "%.1f")):
                series = []
                for m in fits3:
                    g = d3[d3["mode"] == m].sort_values("q")
                    fk, fa = (fits3[m]["m_k"], fits3[m]["m_alpha"]) if ykey == "mem" \
                             else (fits3[m]["t_k"], fits3[m]["t_alpha"])
                    lbl = "%s\n%s ∝ n_quad^%.2f" % (
                        MODE_SHORT[m], "MiB" if ykey == "mem" else "s", fa)
                    curve = (lambda nq, k=fk, a=fa: k * nq ** a)
                    series.append((m, g.n_quad.values, g[ykey].values, curve, lbl))
                note = ("x-axis is n_quad = N_fe x q (N_fe=%d fixed).  scRPA: %d outer "
                        "cycles at every q; RPA@DFT: 1 (no outer loop)."
                        % (C3.FIXED["finite_element_number"],
                           fits3.get("sc", {}).get("outer", 0)))
                full_page(pdf, lambda ax, fig, s=series, yl=ylabel, t=page_title, n=note, f=fmt:
                          line_plot(ax, s, "n_quad", yl, t, n, value_fmt=f))
                md += ["## %s" % page_title, "", H3, "", note, ""]

            # ---- page 12: phase-3 numeric ----
            note3 = ("Energy spread across the whole q sweep: scRPA %.2e Ha, RPA@DFT "
                     "%.2e Ha -- q is converged at q=20 for both modes at this grid."
                     % (fits3.get("sc", {}).get("spread", float("nan")),
                        fits3.get("nsc", {}).get("spread", float("nan"))))
            table_page(pdf, "Phase 3 -- numeric data", H3, phase3_table(d3),
                      extra_note=note3)
            md += ["## Phase 3 data", "", note3, "", to_md(phase3_table(d3)), ""]

    # ---- csv ----
    pd.DataFrame([dict(mode=MODE_LABEL[m],
                       mem_base=fits[m]["mem_base"], mem_quad=fits[m]["mem_quad"],
                       memL_base=fits[m]["memL_base"], memL_slope=fits[m]["memL_slope"],
                       t_k=fits[m]["t_k"], t_alpha=fits[m]["t_alpha"],
                       outer_at_anchor=fits[m]["outer_at_anchor"])
                  for m in fits]).to_csv(os.path.join(OUT, "scaling_fits.csv"), index=False)
    proj.to_csv(os.path.join(OUT, "projections.csv"), index=False)
    ratio.to_csv(os.path.join(OUT, "time_ratio.csv"), index=False)
    with open(os.path.join(OUT, "scaling_report.md"), "w") as fh:
        fh.write("\n".join(md) + "\n")

    for n in ("scaling_report.pdf", "scaling_report.md", "scaling_fits.csv",
              "projections.csv", "time_ratio.csv"):
        print("wrote %s" % os.path.join(OUT, n))
    print()
    print(to_md(proj))
    print()
    print(to_md(ratio))


if __name__ == "__main__":
    main()
