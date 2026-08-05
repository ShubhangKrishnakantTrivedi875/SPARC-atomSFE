"""
Walk results/ and build the phase-1 comparison tables.

    python collect.py                      # all tables to summary/
    python collect.py --xc RPA --setting loose

Writes to summary/ :

    all_runs.csv          one row per case: energy, timing, memory, status
    energy_compare.csv    total energy per code, side by side, with differences
    eigen_compare.csv     per-state eigenvalues, new vs orig vs oneatom
    memory_compare.csv    peak RSS per code (sampled) + SLURM MaxRSS as a cross-check
    timing_compare.csv    wall time per code
    missing.txt           cases with no result.json yet, and failures

The energy and eigenvalue comparisons use 'new' as the reference column, since
the question phase 1 answers is whether the open-shell and pseudopotential
corrections moved anything relative to the other two codes.
"""

from __future__ import annotations

import argparse
import json
import os
import subprocess

import pandas as pd

import cases as C

SUMMARY = os.path.join(C.HERE, "summary")
CODES   = ["new", "orig", "oneatom"]
KEY     = ["Z", "electrons", "xc", "setting"]


# ---------------------------------------------------------------------------
def load_rows(filt):
    rows, missing, raw = [], [], []
    for case in C.ALL_CASES:
        if any(case[k] != v for k, v in filt.items() if v is not None):
            continue
        path = os.path.join(case["outdir"], "result.json")
        if not os.path.exists(path):
            missing.append((case["cid"], "not run"))
            continue
        with open(path) as fh:
            r = json.load(fh)
        raw.append((case, r))
        if not r.get("ok"):
            missing.append((case["cid"], "FAILED: %s" % (r.get("error") or
                                                         "exit %s" % r.get("exit_code"))))
    # Only array tasks get a SLURM figure.  An interactive run records the enclosing
    # session's job id, and sacct would then report that session's memory (1.3 GB of
    # dashboard) as if it were the case's.
    def acct_key(r):
        s = r.get("slurm") or {}
        return str(s["job_id"]) if s.get("job_id") and s.get("array_task") else None

    sacct = load_sacct(acct_key(r) for _, r in raw)

    for case, r in raw:
        mem = r.get("memory", {})
        acct = sacct.get(acct_key(r), {}) if acct_key(r) else {}
        rows.append(dict(
            cid            = case["cid"],
            code           = case["code"],
            Z              = case["Z"],
            electrons      = case["electrons"],
            xc             = case["xc"],
            setting        = case["setting"],
            ok             = bool(r.get("ok")),
            total_energy   = r.get("total_energy"),
            n_states       = r.get("n_states"),
            eps_homo       = (r["eigenvalues"][-1]
                              if r.get("eigenvalues") else None),
            converged      = r.get("converged"),
            scf_iters      = r.get("scf_iterations"),
            wall_s         = r.get("wall_s"),
            t_solve_s      = r.get("t_solve_s"),
            peak_rss_mb    = mem.get("peak_rss_tree_mb"),
            peak_rss_1p_mb = mem.get("peak_rss_single_mb"),
            slurm_maxrss_mb= acct.get("slurm_maxrss_mb"),
            slurm_elapsed  = acct.get("slurm_elapsed"),
            exit_code      = r.get("exit_code"),
            error          = r.get("error"),
        ))
    return pd.DataFrame(rows), missing


_SACCT_CACHE = {}


def load_sacct(job_ids):
    """
    MaxRSS and Elapsed as SLURM measured them, for a batch of job ids.

    One sacct call for the whole set, keyed by JobIDRaw so it matches the job_id
    each result.json recorded.  MaxRSS lives on the .batch step, not the parent.

    Caveat: SLURM samples at JobAcctGatherFrequency (typically 30 s), so for short
    cases MaxRSS UNDERSTATES the peak -- 137 MB against a sampled and kernel-confirmed
    198 MB on a 14 s case.  Treat it as a cross-check for the long runs, not as the
    primary figure.
    """
    want = sorted({str(j) for j in job_ids if j})
    todo = [j for j in want if j not in _SACCT_CACHE]
    for chunk in [todo[i:i + 400] for i in range(0, len(todo), 400)]:
        try:
            out = subprocess.run(
                ["sacct", "-j", ",".join(chunk), "--units=M", "--noheader", "-P",
                 "--format=JobIDRaw,MaxRSS,Elapsed,State"],
                capture_output=True, text=True, timeout=120).stdout
        except (OSError, subprocess.SubprocessError):
            out = ""
        for line in out.splitlines():
            f = line.split("|")
            if len(f) < 4 or "." not in f[0]:
                continue                          # parent row carries no MaxRSS
            raw, mem, elapsed, state = f[0].split(".")[0], f[1], f[2], f[3]
            if not mem.endswith("M"):
                continue
            mb = float(mem[:-1])
            cur = _SACCT_CACHE.get(raw)
            if cur is None or mb > cur["slurm_maxrss_mb"]:
                _SACCT_CACHE[raw] = {"slurm_maxrss_mb": mb,
                                     "slurm_elapsed": elapsed, "slurm_state": state}
        for j in chunk:
            _SACCT_CACHE.setdefault(j, {})
    return _SACCT_CACHE


def pivot(df, value, agg="first"):
    if df.empty:
        return pd.DataFrame()
    p = df.pivot_table(index=KEY, columns="code", values=value, aggfunc=agg)
    for c in CODES:
        if c not in p.columns:
            p[c] = float("nan")
    return p[CODES].reset_index()


def add_diffs(p, ref="new"):
    for c in CODES:
        if c != ref:
            p["d_%s_vs_%s" % (c, ref)] = p[c] - p[ref]
    return p


# ---------------------------------------------------------------------------
def eigen_table(filt):
    rows = []
    for case in C.ALL_CASES:
        if any(case[k] != v for k, v in filt.items() if v is not None):
            continue
        path = os.path.join(case["outdir"], "result.json")
        if not os.path.exists(path):
            continue
        with open(path) as fh:
            r = json.load(fh)
        eig = r.get("eigenvalues") or []
        nv, lv = r.get("n_values") or [], r.get("l_values") or []
        for i, e in enumerate(eig):
            rows.append(dict(
                Z=case["Z"], electrons=case["electrons"], xc=case["xc"],
                setting=case["setting"], code=case["code"], state=i,
                n=(nv[i] if i < len(nv) else None),
                l=(lv[i] if i < len(lv) else None),
                eigenvalue=e))
    if not rows:
        return pd.DataFrame()
    df = pd.DataFrame(rows)
    p = df.pivot_table(index=KEY + ["state", "n", "l"], columns="code",
                       values="eigenvalue", aggfunc="first")
    for c in CODES:
        if c not in p.columns:
            p[c] = float("nan")
    return add_diffs(p[CODES].reset_index())


# ---------------------------------------------------------------------------
def main():
    ap = argparse.ArgumentParser()
    for k in ("xc", "setting", "electrons"):
        ap.add_argument("--" + k, default=None)
    ap.add_argument("--Z", type=int, default=None)
    a = ap.parse_args()
    filt = dict(xc=a.xc, setting=a.setting, electrons=a.electrons, Z=a.Z)

    os.makedirs(SUMMARY, exist_ok=True)
    df, missing = load_rows(filt)

    if df.empty:
        print("no result.json found yet under %s" % C.RESULTS_ROOT)
    else:
        df.sort_values(["xc", "setting", "Z", "electrons", "code"]).to_csv(
            os.path.join(SUMMARY, "all_runs.csv"), index=False)

        add_diffs(pivot(df[df.ok], "total_energy")).to_csv(
            os.path.join(SUMMARY, "energy_compare.csv"), index=False)
        pivot(df, "wall_s").to_csv(
            os.path.join(SUMMARY, "timing_compare.csv"), index=False)
        mem = pivot(df, "peak_rss_mb")
        if not mem.empty:
            mem["ratio_orig_over_new"]    = mem["orig"] / mem["new"]
            mem["ratio_oneatom_over_new"] = mem["oneatom"] / mem["new"]
            # SLURM's own figure alongside, suffixed, so the two are never confused
            sl = pivot(df, "slurm_maxrss_mb")
            if not sl.empty:
                mem = mem.merge(sl.rename(columns={c: "slurm_" + c for c in CODES}),
                                on=KEY, how="left")
        mem.to_csv(os.path.join(SUMMARY, "memory_compare.csv"), index=False)

        et = eigen_table(filt)
        if not et.empty:
            et.to_csv(os.path.join(SUMMARY, "eigen_compare.csv"), index=False)

        # ---- terminal digest ----
        done = int(df.ok.sum())
        print("%d/%d cases have results, %d ok\n" % (len(df), len(C.ALL_CASES), done))
        ok = df[df.ok]
        if not ok.empty:
            print("--- total energy (Ha) ---")
            print(add_diffs(pivot(ok, "total_energy")).to_string(index=False))
            print("\n--- peak RSS (MB, whole process tree) ---")
            print(mem.to_string(index=False))
            print("\n--- wall time (s) ---")
            print(pivot(df, "wall_s").to_string(index=False))

    with open(os.path.join(SUMMARY, "missing.txt"), "w") as fh:
        for cid, why in missing:
            fh.write("%-58s %s\n" % (cid, why))
    fails = [m for m in missing if m[1] != "not run"]
    print("\n%d cases not run, %d failed -- see summary/missing.txt"
          % (len(missing) - len(fails), len(fails)))


if __name__ == "__main__":
    main()
