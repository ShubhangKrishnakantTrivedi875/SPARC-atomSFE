"""
Clean phase-2 tables: peak memory and wall time, corrected vs original.

    python make_tables.py

Writes to summary/tables/ :

    phase2_memory_timing.csv    both sweeps, one row per configuration
    phase2_memory_timing.md     the same, as two markdown tables
    phase2_memory_timing.pdf    the same, typeset, one page

Eight columns per sweep -- the swept variable, the grid size, peak memory for
each code with the reduction, wall time for each code with the speedup:

    N_fe | n_quad | new (MiB) | orig (MiB) | reduction | new (s) | orig (s) | speedup

Peak memory is the sampled maximum resident set size, cross-checked against the
kernel high-water mark and SLURM MaxRSS (all three agree to ~1%).
"""

from __future__ import annotations

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd

import cases as C

OUT = os.path.join(C.HERE, "summary", "tables")
MIB = 1.048576

# The two comparisons drawn.  Each is (a, b, a_label, b_label): `a` is the column
# the reduction/speedup is measured *for*, `b` is the baseline it is measured
# *against*.  A (mode, code) pair identifies a run uniquely.
COMPARISONS = {
    "sc":  (("sc", "new"), ("sc", "orig"),  "new", "orig"),
    "nsc": (("nsc", "new"), ("sc", "new"),  "RPA@DFT", "scRPA"),
}


def head_for(a_label, b_label):
    return ["n_quad", "%s (MiB)" % a_label, "%s (MiB)" % b_label, "reduction",
            "%s (s)" % a_label, "%s (s)" % b_label, "speedup", "dE (Ha)"]


def load():
    rows = []
    for c in C.ALL_CASES:
        p = os.path.join(c["outdir"], "result.json")
        if not os.path.exists(p):
            continue
        r = json.load(open(p))
        if not r.get("ok"):
            continue
        peak = (r.get("memory") or {}).get("peak_rss_tree_mb")
        rows.append(dict(mode=c["mode"], code=c["code"], n_fe=c["n_fe"], l_max=c["l_max"],
                         n_quad=c["n_fe"] * C.FIXED["quadrature_point_number"],
                         peak_mib=(peak / MIB) if peak else None,
                         wall_s=r.get("wall_s"), energy=r.get("total_energy")))
    return pd.DataFrame(rows)


def sweep(df, pairs, axis, comparison):
    """One tidy row per configuration, ordered along `axis`."""
    (a_mode, a_code), (b_mode, b_code), a_label, b_label = COMPARISONS[comparison]
    head = head_for(a_label, b_label)
    out = []
    for n_fe, l_max in sorted(pairs, key=lambda k: k[0] if axis == "n_fe" else k[1]):
        g = df[(df.n_fe == n_fe) & (df.l_max == l_max)]
        a = g[(g["mode"] == a_mode) & (g.code == a_code)]
        b = g[(g["mode"] == b_mode) & (g.code == b_code)]
        if a.empty:
            continue
        a = a.iloc[0]
        row = {axis: (n_fe if axis == "n_fe" else l_max), "n_quad": int(a.n_quad),
               "%s (MiB)" % a_label: a.peak_mib, "%s (s)" % a_label: a.wall_s}
        if b.empty:
            row.update({"%s (MiB)" % b_label: None, "reduction": None,
                        "%s (s)" % b_label: None, "speedup": None, "dE (Ha)": None})
        else:
            b = b.iloc[0]
            row.update({"%s (MiB)" % b_label: b.peak_mib,
                        "reduction": 100.0 * (a.peak_mib - b.peak_mib) / b.peak_mib,
                        "%s (s)" % b_label: b.wall_s,
                        "speedup": b.wall_s / a.wall_s,
                        "dE (Ha)": a.energy - b.energy})
        out.append(row)
    cols = [axis] + head
    return pd.DataFrame(out)[cols] if out else pd.DataFrame(columns=cols)


def fmt(k, v):
    if v is None or pd.isna(v):
        return "-"
    if k in ("n_fe", "l_max", "n_quad"):
        return "%d" % v
    if k == "reduction":
        return "%+.1f%%" % v
    if k == "speedup":
        return "%.2fx" % v
    if k == "dE (Ha)":
        return "%+.2e" % v
    if k.endswith("(MiB)"):
        return "%,.0f".replace(",", "") % v if v < 1000 else "%.0f" % v
    return "%.1f" % v


def to_md(d, axis, label):
    keys = list(d.columns)
    head = [label] + keys[1:]
    lines = ["| " + " | ".join(head) + " |",
             "|" + "|".join(["---:"] * len(head)) + "|"]
    for _, r in d.iterrows():
        lines.append("| " + " | ".join(fmt(k, r[k]) for k in keys) + " |")
    return "\n".join(lines)


def pdf_table(ax, d, axis, label, title):
    ax.axis("off")
    ax.set_title(title, fontsize=11, fontweight="bold", loc="left", pad=14)
    keys = list(d.columns)
    cols = [label] + keys[1:]
    cells = [[fmt(k, r[k]) for k in keys] for _, r in d.iterrows()]
    t = ax.table(cellText=cells, colLabels=cols, cellLoc="right", loc="upper center")
    t.auto_set_font_size(False)
    t.set_fontsize(8.5)
    t.scale(1, 1.55)
    for j in range(len(cols)):                      # header styling
        c = t[0, j]
        c.set_facecolor("#e8e8e8")
        c.set_text_props(fontweight="bold")
    for i in range(1, len(cells) + 1):              # zebra + emphasis
        for j in range(len(cols)):
            if i % 2 == 0:
                t[i, j].set_facecolor("#f7f7f7")
            if cols[j] in ("reduction", "speedup"):
                t[i, j].set_text_props(fontweight="bold")
    for key, cell in t.get_celld().items():
        cell.set_linewidth(0.4)
        cell.set_edgecolor("#bbbbbb")


FIXED_SUB = ("Mercury (Z=80), all-electron, one outer cycle.  Fixed: p=20, q=%d, "
             "domain=%.0f Bohr, omega=%d, polynomial mesh (concentration 2)."
             % (C.FIXED["quadrature_point_number"], C.FIXED["domain_size"],
                C.FIXED["frequency_quadrature_point_number"]))

# comparison -> (file stem, page title, subtitle tail, footnote)
REPORTS = {
    "sc": ("phase2_memory_timing",
           "Phase 2  —  self-consistent RPA, peak memory and wall time",
           "Self-consistent-potential RPA (xc='RPA'), corrected code vs original.",
           "reduction = (new - orig)/orig on peak memory    speedup = orig/new on wall time"
           "    dE = new - orig on total energy"),
    "nsc": ("phase2_memory_timing_nonselfconsistent",
            "Phase 2  —  non-self-consistent RPA, peak memory and wall time",
            "E_c^RPA on frozen GGA_PBE orbitals (xc='RPA@DFT'), against self-consistent "
            "RPA in the corrected code.  No OEP equation is solved and no correlation "
            "energy density is built, so dE is physics, not error.",
            "reduction = (RPA@DFT - scRPA)/scRPA on peak memory    speedup = scRPA/RPA@DFT "
            "on wall time    dE = RPA@DFT - scRPA on total energy"),
}


def report(df, comparison):
    stem_name, title, sub_tail, foot = REPORTS[comparison]
    mesh = sweep(df, C.MESH_SWEEP,    "n_fe",  comparison)
    ang  = sweep(df, C.ANGULAR_SWEEP, "l_max", comparison)
    if mesh.empty and ang.empty:
        print("no completed %s cases, skipped" % comparison)
        return
    stem = os.path.join(OUT, stem_name)
    sub  = FIXED_SUB + "  " + sub_tail

    # ---- csv: every block, tagged ----
    tagged = []
    for d, axis, tag in ((mesh, "n_fe",  "mesh (Lmax=%d)" % C.MESH_SWEEP[0][1]),
                         (ang,  "l_max", "angular (N_fe=%d)" % C.ANGULAR_SWEEP[0][0])):
        if d.empty:
            continue
        d = d.copy().rename(columns={axis: "swept_value"})
        d.insert(0, "sweep", tag)
        tagged.append(d)
    pd.concat(tagged, ignore_index=True).to_csv(stem + ".csv", index=False)

    blocks = [(d, axis, label, t) for d, axis, label, t in (
        (mesh, "n_fe",  "N_fe", "Radial mesh sweep  (Lmax = %d)" % C.MESH_SWEEP[0][1]),
        (ang,  "l_max", "Lmax", "Angular sweep  (N_fe = %d)" % C.ANGULAR_SWEEP[0][0]),
    ) if not d.empty]

    # ---- markdown ----
    with open(stem + ".md", "w") as fh:
        fh.write("# %s\n\n%s\n\n" % (title.replace("  —  ", " -- "), sub))
        for d, axis, label, t in blocks:
            fh.write("## %s\n\n%s\n\n" % (t.replace("  ", " "), to_md(d, axis, label)))
        fh.write(foot + "\n")

    # ---- pdf ----
    fig = plt.figure(figsize=(11.0, 3.0 + 2.2 * len(blocks)))
    fig.suptitle(title, fontsize=13, fontweight="bold", x=0.06, ha="left", y=0.97)
    fig.text(0.06, 0.90, sub, fontsize=8, va="top", wrap=True)
    top, h = 0.80, 0.72 / len(blocks)
    for i, (d, axis, label, t) in enumerate(blocks):
        pdf_table(fig.add_axes([0.05, top - (i + 1) * h + 0.05, 0.90, h - 0.07]),
                  d, axis, label, t)
    fig.text(0.06, 0.035, foot + "\nPeak memory is the sampled maximum resident set "
             "size; kernel ru_maxrss and SLURM MaxRSS agree with it to ~1%.",
             fontsize=7.2, va="top", color="#444444")
    fig.savefig(stem + ".pdf", format="pdf")
    plt.close(fig)

    for e in ("csv", "md", "pdf"):
        print("wrote %s.%s" % (stem, e))
    for d, axis, label, _ in blocks:
        print()
        print(to_md(d, axis, label))


def main():
    os.makedirs(OUT, exist_ok=True)
    df = load()
    if df.empty:
        print("no completed cases under %s" % C.RESULTS_ROOT)
        return
    for comparison in C.MODE_ORDER:
        print("\n=== %s ===" % C.MODES[comparison]["tag"])
        report(df, comparison)


if __name__ == "__main__":
    main()
