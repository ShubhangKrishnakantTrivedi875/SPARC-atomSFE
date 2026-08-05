"""
Build the phase-3 quadrature-scaling tables.

    python collect.py            # tables to summary/ + terminal digest
    python collect.py --md       # also summary/tables.md

Writes to summary/ :

    all_runs.csv     one row per case
    scaling.csv      q, n_quad, peak memory, wall time, energy -- per RPA mode
    compare.csv      RPA@DFT vs self-consistent RPA at each q
    tables.md        markdown for the writeup

Both modes sweep the same q list on the same grid, so scaling.csv answers "is q
converged, and how do time and memory grow", and compare.csv answers "what does
dropping self-consistency buy".  dE in compare.csv is a physics difference between
two functionals, not an error.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess

import pandas as pd

import cases as C

SUMMARY = os.path.join(C.HERE, "summary")
MIB     = 1.048576


def sacct_maxrss(job_ids):
    want = sorted({str(j) for j in job_ids if j})
    if not want:
        return {}
    try:
        out = subprocess.run(
            ["sacct", "-j", ",".join(want), "--units=M", "--noheader", "-P",
             "--format=JobIDRaw,MaxRSS"], capture_output=True, text=True, timeout=120).stdout
    except (OSError, subprocess.SubprocessError):
        return {}
    got = {}
    for line in out.splitlines():
        f = line.split("|")
        if len(f) < 2 or "." not in f[0] or not f[1].endswith("M"):
            continue
        raw, mb = f[0].split(".")[0], float(f[1][:-1])
        got[raw] = max(got.get(raw, 0.0), mb)
    return got


def load():
    raw, missing = [], []
    for case in C.ALL_CASES:
        path = os.path.join(case["outdir"], "result.json")
        if not os.path.exists(path):
            missing.append((case["cid"], "not run"))
            continue
        with open(path) as fh:
            raw.append((case, json.load(fh)))
    acct = sacct_maxrss((r.get("slurm") or {}).get("job_id")
                        for _, r in raw if (r.get("slurm") or {}).get("array_task"))
    rows = []
    for case, r in raw:
        if not r.get("ok"):
            missing.append((case["cid"], "FAILED: %s" %
                            (r.get("error") or "exit %s" % r.get("exit_code"))))
        mem, s = r.get("memory", {}), r.get("slurm") or {}
        peak = mem.get("peak_rss_tree_mb")
        rows.append(dict(
            cid       = case["cid"],
            mode      = case["mode"],
            variant   = case["variant"],
            q         = case["q"],
            n_quad    = C.FIXED["finite_element_number"] * case["q"],
            ok        = bool(r.get("ok")),
            E         = r.get("total_energy"),
            eps_homo  = (r["eigenvalues"][-1] if r.get("eigenvalues") else None),
            peak_mib  = (peak / MIB) if peak else None,
            kernel_mb = mem.get("ru_maxrss_child_mb"),
            slurm_mb  = acct.get(str(s.get("job_id"))) if s.get("array_task") else None,
            wall_s    = r.get("wall_s"),
            solve_s   = r.get("t_solve_s"),
            scf_iters = r.get("scf_iterations"),
            converged = r.get("converged"),
        ))
    return pd.DataFrame(rows), missing


def scaling(df):
    """Per-mode table with the growth exponent between consecutive q."""
    out = []
    for m in [k for k in C.MODE_ORDER if k in set(df["mode"])]:
        d = df[(df["mode"] == m) & df.ok].sort_values("q")
        prev = None
        for _, r in d.iterrows():
            row = dict(mode=C.MODES[m]["tag"], q=int(r.q), n_quad=int(r.n_quad),
                       peak_mib=round(r.peak_mib, 1) if r.peak_mib else None,
                       wall_s=round(r.wall_s, 1) if r.wall_s else None,
                       scf_iters=r.scf_iters, E=r.E)
            if prev is not None and prev.peak_mib and r.peak_mib:
                import math
                lr = math.log(r.n_quad / prev.n_quad)
                row["mem_exponent"]  = round(math.log(r.peak_mib / prev.peak_mib) / lr, 2)
                if prev.wall_s and r.wall_s:
                    row["time_exponent"] = round(math.log(r.wall_s / prev.wall_s) / lr, 2)
            out.append(row)
            prev = r
    return pd.DataFrame(out)


def compare(df):
    """RPA@DFT vs self-consistent RPA at the q values where both ran."""
    out = []
    for q in sorted(df.q.unique()):
        a = df[(df["mode"] == "nsc") & (df.q == q) & df.ok]
        b = df[(df["mode"] == "sc")  & (df.q == q) & df.ok]
        if a.empty or b.empty:
            continue
        a, b = a.iloc[0], b.iloc[0]
        out.append(dict(
            q=int(q), n_quad=int(a.n_quad),
            nsc_mib=round(a.peak_mib, 1), sc_mib=round(b.peak_mib, 1),
            mem_rel_pct=round(100.0 * (a.peak_mib - b.peak_mib) / b.peak_mib, 1),
            nsc_s=round(a.wall_s, 1), sc_s=round(b.wall_s, 1),
            speedup=round(b.wall_s / a.wall_s, 2),
            nsc_E=a.E, sc_E=b.E, dE_Ha=a.E - b.E,
            d_eps_homo=(a.eps_homo - b.eps_homo)
                       if a.eps_homo is not None and b.eps_homo is not None else None))
    return pd.DataFrame(out)


def as_md(d, cols):
    cols = [(h, k) for h, k in cols if k in d.columns]
    lines = ["| " + " | ".join(h for h, _ in cols) + " |",
             "|" + "|".join(["---:"] * len(cols)) + "|"]
    for _, r in d.iterrows():
        cells = []
        for h, k in cols:
            v = r.get(k)
            if v is None or pd.isna(v):
                cells.append("-")
            elif k.endswith("_pct"):
                cells.append("%+.1f%%" % v)
            elif k in ("q", "n_quad", "scf_iters"):
                cells.append("%d" % v)
            elif k in ("E", "nsc_E", "sc_E", "dE_Ha", "d_eps_homo"):
                cells.append("%.9f" % v if abs(v) < 1 else "%.6f" % v)
            else:
                cells.append("%.1f" % v if isinstance(v, float) else str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--md", action="store_true")
    a = ap.parse_args()
    os.makedirs(SUMMARY, exist_ok=True)
    df, missing = load()
    if df.empty:
        print("no results yet under %s" % C.RESULTS_ROOT)
        return

    df.sort_values(["mode", "q"]).to_csv(
        os.path.join(SUMMARY, "all_runs.csv"), index=False)
    sc, cmp_ = scaling(df), compare(df)
    sc.to_csv(os.path.join(SUMMARY, "scaling.csv"), index=False)
    cmp_.to_csv(os.path.join(SUMMARY, "compare.csv"), index=False)

    print("Z=%d, all-electron, N_fe=%d, p=%d, Lmax=%d, omega=%d;  %d/%d cases done\n"
          % (C.Z, C.FIXED["finite_element_number"], C.FIXED["polynomial_order"],
             C.FIXED["angular_momentum_cutoff"],
             C.FIXED["frequency_quadrature_point_number"],
             len(df), len(C.ALL_CASES)))
    if not sc.empty:
        print("--- scaling with quadrature order ---")
        print(sc.to_string(index=False), "\n")
    if not cmp_.empty:
        print("--- RPA@DFT vs self-consistent RPA ---")
        print(cmp_.to_string(index=False), "\n")
    print("cross-check sampled / kernel / SLURM (MB):")
    print(df[["cid", "peak_mib", "kernel_mb", "slurm_mb"]].to_string(index=False))

    if a.md:
        with open(os.path.join(SUMMARY, "tables.md"), "w") as fh:
            fh.write("## Scaling with radial quadrature order "
                     "(Ar, all-electron, N_fe=%d, p=%d, "
                     "Lmax=%d, omega=%d)\n\n"
                     % (C.FIXED["finite_element_number"],
                        C.FIXED["polynomial_order"],
                        C.FIXED["angular_momentum_cutoff"],
                        C.FIXED["frequency_quadrature_point_number"]))
            fh.write(as_md(sc, [("mode", "mode"), ("q", "q"), ("n_quad", "n_quad"),
                                ("peak (MiB)", "peak_mib"), ("wall (s)", "wall_s"),
                                ("SCF its", "scf_iters"),
                                ("mem exp", "mem_exponent"),
                                ("time exp", "time_exponent")]) + "\n\n")
            if not cmp_.empty:
                fh.write("## RPA@DFT vs self-consistent RPA\n\n")
                fh.write(as_md(cmp_, [("q", "q"), ("n_quad", "n_quad"),
                                      ("RPA@DFT (MiB)", "nsc_mib"), ("scRPA (MiB)", "sc_mib"),
                                      ("rel", "mem_rel_pct"), ("RPA@DFT (s)", "nsc_s"),
                                      ("scRPA (s)", "sc_s"), ("speedup", "speedup"),
                                      ("dE (Ha)", "dE_Ha")]) + "\n")
                fh.write("\ndE is the physics difference between self-consistent RPA "
                         "and E_c^RPA on frozen GGA_PBE orbitals, not an error.\n")
        print("\nwrote %s" % os.path.join(SUMMARY, "tables.md"))

    with open(os.path.join(SUMMARY, "missing.txt"), "w") as fh:
        for cid, why in missing:
            fh.write("%-24s %s\n" % (cid, why))
    fails = [m for m in missing if m[1] != "not run"]
    print("\n%d not run, %d failed" % (len(missing) - len(fails), len(fails)))


if __name__ == "__main__":
    main()
