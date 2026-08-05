"""
Build the phase-2 memory and timing tables.

    python collect.py            # tables to summary/ + terminal digest
    python collect.py --md       # also emit summary/tables.md for the PR writeup

Writes to summary/ :

    all_runs.csv        one row per case
    mesh_sweep.csv      N_fe = 10, 15, 25   at Lmax = 20
    angular_sweep.csv   Lmax =  5, 10, 20   at N_fe = 25
    tables.md           the same two sweeps as markdown, ready to paste

Peak memory is the sampled maximum of resident set size over the whole run,
cross-checked against the kernel's own high-water mark (ru_maxrss) and, for array
tasks, against SLURM's MaxRSS.  Note SLURM polls at JobAcctGatherFrequency
(~30 s) so its figure understates short runs; the sampled one is primary.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess

import pandas as pd

import cases as C

SUMMARY = os.path.join(C.HERE, "summary")
MIB     = 1.048576          # MB (10^6) -> MiB (2^20) is /1.048576; psutil gives bytes


def sacct_maxrss(job_ids):
    want = sorted({str(j) for j in job_ids if j})
    if not want:
        return {}
    try:
        out = subprocess.run(
            ["sacct", "-j", ",".join(want), "--units=M", "--noheader", "-P",
             "--format=JobIDRaw,MaxRSS,Elapsed,State"],
            capture_output=True, text=True, timeout=120).stdout
    except (OSError, subprocess.SubprocessError):
        return {}
    got = {}
    for line in out.splitlines():
        f = line.split("|")
        if len(f) < 4 or "." not in f[0] or not f[1].endswith("M"):
            continue
        raw, mb = f[0].split(".")[0], float(f[1][:-1])
        got[raw] = max(got.get(raw, 0.0), mb)
    return got


def load():
    rows, missing = [], []
    raw = []
    for case in C.ALL_CASES:
        path = os.path.join(case["outdir"], "result.json")
        if not os.path.exists(path):
            missing.append((case["cid"], "not run"))
            continue
        with open(path) as fh:
            raw.append((case, json.load(fh)))

    acct = sacct_maxrss(
        (r.get("slurm") or {}).get("job_id")
        for _, r in raw if (r.get("slurm") or {}).get("array_task"))

    for case, r in raw:
        if not r.get("ok"):
            missing.append((case["cid"], "FAILED: %s" %
                            (r.get("error") or "exit %s" % r.get("exit_code"))))
        mem = r.get("memory", {})
        s   = r.get("slurm") or {}
        peak_mb = mem.get("peak_rss_tree_mb")
        rows.append(dict(
            cid          = case["cid"],
            code         = case["code"],
            n_fe         = case["n_fe"],
            l_max        = case["l_max"],
            n_quad       = case["n_fe"] * C.FIXED["quadrature_point_number"],
            ok           = bool(r.get("ok")),
            total_energy = r.get("total_energy"),
            peak_mb      = peak_mb,
            peak_mib     = (peak_mb / MIB) if peak_mb else None,
            kernel_mb    = mem.get("ru_maxrss_child_mb"),
            slurm_mb     = acct.get(str(s.get("job_id"))) if s.get("array_task") else None,
            wall_s       = r.get("wall_s"),
            solve_s      = r.get("t_solve_s"),
            n_states     = r.get("n_states"),
        ))
    return pd.DataFrame(rows), missing


def sweep(df, subset, axis):
    """One sweep table: new vs orig, with deltas, along `axis`."""
    d = df[df.apply(lambda r: (r.n_fe, r.l_max) in subset, axis=1)]
    if d.empty:
        return pd.DataFrame()
    out = []
    for key in sorted({(r.n_fe, r.l_max) for _, r in d.iterrows()},
                      key=lambda k: k[0] if axis == "n_fe" else k[1]):
        n_fe, l_max = key
        row = {"n_fe": n_fe, "l_max": l_max,
               "n_quad": n_fe * C.FIXED["quadrature_point_number"]}
        for code in ("new", "orig"):
            m = d[(d.n_fe == n_fe) & (d.l_max == l_max) & (d.code == code)]
            row["%s_mib" % code]  = round(m.peak_mib.iloc[0], 2) if len(m) else None
            row["%s_wall_s" % code] = round(m.wall_s.iloc[0], 1) if len(m) else None
            row["%s_E" % code]      = m.total_energy.iloc[0] if len(m) else None
        if row["new_mib"] and row["orig_mib"]:
            row["delta_mib"] = round(row["orig_mib"] - row["new_mib"], 2)
            row["mem_rel_%"] = round(100.0 * (row["new_mib"] - row["orig_mib"])
                                     / row["orig_mib"], 1)
        if row["new_wall_s"] and row["orig_wall_s"]:
            row["speedup"] = round(row["orig_wall_s"] / row["new_wall_s"], 2)
        if row["new_E"] is not None and row["orig_E"] is not None:
            row["dE_Ha"] = row["orig_E"] - row["new_E"]
        out.append(row)
    return pd.DataFrame(out)


def as_md(d, axis_key, axis_label):
    if d.empty:
        return "_no results for this sweep yet_"
    cols = [(axis_label, axis_key),
            ("n_quad", "n_quad"), ("new (MiB)", "new_mib"), ("orig (MiB)", "orig_mib"),
            ("Delta (MiB)", "delta_mib"), ("Rel.", "mem_rel_%"),
            ("new (s)", "new_wall_s"), ("orig (s)", "orig_wall_s"), ("speedup", "speedup")]
    cols = [(h, k) for h, k in cols if k in d.columns]
    lines = ["| " + " | ".join(h for h, _ in cols) + " |",
             "|" + "|".join(["---:"] * len(cols)) + "|"]
    for _, r in d.iterrows():
        cells = []
        for h, k in cols:
            v = r.get(k)
            if v is None or pd.isna(v):
                cells.append("-")
            elif k == "mem_rel_%":
                cells.append("%+.1f%%" % v)
            elif k == "delta_mib":
                cells.append("%+.2f" % v)
            elif k in ("n_fe", "l_max", "n_quad"):
                cells.append("%d" % v)
            elif isinstance(v, float):
                cells.append("%.2f" % v if "mib" in k else "%.1f" % v)
            else:
                cells.append(str(v))
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--md", action="store_true", help="also write summary/tables.md")
    a = ap.parse_args()

    os.makedirs(SUMMARY, exist_ok=True)
    df, missing = load()
    if df.empty:
        print("no results yet under %s" % C.RESULTS_ROOT)
        return

    df.sort_values(["l_max", "n_fe", "code"]).to_csv(
        os.path.join(SUMMARY, "all_runs.csv"), index=False)
    mesh = sweep(df, set(C.MESH_SWEEP), "n_fe")
    ang  = sweep(df, set(C.ANGULAR_SWEEP), "l_max")
    mesh.to_csv(os.path.join(SUMMARY, "mesh_sweep.csv"), index=False)
    ang.to_csv(os.path.join(SUMMARY, "angular_sweep.csv"), index=False)

    print("Z=%d, %s, all-electron, %d outer cycle;  %d/%d cases done\n"
          % (C.Z, C.XC, C.OUTER_CYCLES, len(df), len(C.ALL_CASES)))
    if not mesh.empty:
        print("--- mesh sweep (Lmax = 20) ---")
        print(mesh.to_string(index=False), "\n")
    if not ang.empty:
        print("--- angular sweep (N_fe = 25) ---")
        print(ang.to_string(index=False), "\n")
    print("cross-check, sampled vs kernel vs SLURM (MB):")
    print(df[["cid", "peak_mb", "kernel_mb", "slurm_mb"]].to_string(index=False))

    if a.md:
        with open(os.path.join(SUMMARY, "tables.md"), "w") as fh:
            fh.write("## Peak memory vs number of finite elements "
                     "(Hg, all-electron, Lmax = 20)\n\n")
            fh.write(as_md(mesh, "n_fe", "N_fe") + "\n\n")
            fh.write("## Peak memory vs angular-momentum cutoff "
                     "(Hg, all-electron, N_fe = 25)\n\n")
            fh.write(as_md(ang, "l_max", "Lmax") + "\n")
        print("\nwrote %s" % os.path.join(SUMMARY, "tables.md"))

    with open(os.path.join(SUMMARY, "missing.txt"), "w") as fh:
        for cid, why in missing:
            fh.write("%-30s %s\n" % (cid, why))
    fails = [m for m in missing if m[1] != "not run"]
    print("\n%d not run, %d failed" % (len(missing) - len(fails), len(fails)))


if __name__ == "__main__":
    main()
