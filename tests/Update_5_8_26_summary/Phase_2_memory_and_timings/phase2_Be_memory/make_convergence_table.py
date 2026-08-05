"""
Energy and eigenvalue convergence vs N_fe, for one element's phase-2 study.

    python make_convergence_table.py

Writes convergence.csv / convergence.md to summary/.

Each row is compared to the row immediately before it (no fitted/extrapolated
reference, matching the phase-2 (Hg) scaling report's convention) -- the first
row of each mode has nothing to compare against, so its diff is blank.

Eigenvalue column differs by mode, since a single number is needed per row:

    self-consistent RPA ('sc')       HOMO eigenvalue, and its diff vs previous
    non-self-consistent RPA ('nsc')  max |eigenvalue diff| vs the PREVIOUS row,
                                      over every occupied state -- the strictest
                                      statement available, not just the HOMO
"""

from __future__ import annotations

import json
import os

import pandas as pd

import cases as C

OUT = os.path.join(C.HERE, "summary")


def load(mode):
    rows = []
    for c in sorted((c for c in C.ALL_CASES if c["mode"] == mode),
                    key=lambda c: c["n_fe"]):
        p = os.path.join(c["outdir"], "result.json")
        if not os.path.exists(p):
            continue
        r = json.load(open(p))
        if not r.get("ok"):
            continue
        eig = r.get("eigenvalues") or []
        rows.append(dict(n_fe=c["n_fe"], energy=r.get("total_energy"), eig=tuple(eig)))
    return rows


def table(mode, label):
    rows = load(mode)
    out, prev = [], None
    for r in rows:
        row = {"N_fe": r["n_fe"], "E (Ha)": r["energy"],
               "dE vs previous (Ha)": None if prev is None else r["energy"] - prev["energy"]}
        if mode == "sc":
            row["HOMO (Ha)"] = r["eig"][-1] if r["eig"] else None
            row["dHOMO vs previous (Ha)"] = (
                None if prev is None or not r["eig"] or not prev["eig"]
                else r["eig"][-1] - prev["eig"][-1])
        else:
            row["max |dE_eig| vs previous (Ha)"] = (
                None if prev is None or len(r["eig"]) != len(prev["eig"]) or not r["eig"]
                else max(abs(a - b) for a, b in zip(r["eig"], prev["eig"])))
        out.append(row)
        prev = r
    return pd.DataFrame(out)


def fmt(k, v):
    if v is None or (isinstance(v, float) and pd.isna(v)):
        return "-"
    if k == "N_fe":
        return "%d" % v
    if k == "E (Ha)" or k == "HOMO (Ha)":
        return "%.9f" % v
    return "%+.2e" % v


def to_md(d):
    keys = list(d.columns)
    lines = ["| " + " | ".join(keys) + " |",
             "|" + "|".join(["---:"] * len(keys)) + "|"]
    for _, r in d.iterrows():
        lines.append("| " + " | ".join(fmt(k, r[k]) for k in keys) + " |")
    return "\n".join(lines)


def main():
    os.makedirs(OUT, exist_ok=True)
    sc  = table("sc",  "self-consistent RPA")
    nsc = table("nsc", "non-self-consistent RPA@DFT")

    tagged = []
    for mode, d in (("self_consistent_rpa", sc), ("non_self_consistent_rpa", nsc)):
        t = d.copy()
        t.insert(0, "mode", mode)
        tagged.append(t)
    pd.concat(tagged, ignore_index=True).to_csv(
        os.path.join(OUT, "convergence.csv"), index=False)

    with open(os.path.join(OUT, "convergence.md"), "w") as fh:
        fh.write("# Energy and eigenvalue convergence vs N_fe (Z=%d)\n\n" % C.Z)
        fh.write("Lmax=%d, p=%d, q=%d, domain=%.0f Bohr, omega=%d, all-electron.  "
                 "Each row is compared to the row directly above it -- no fitted or "
                 "extrapolated reference.\n\n"
                 % (C.L_MAX, C.FIXED["polynomial_order"], C.FIXED["quadrature_point_number"],
                    C.FIXED["domain_size"], C.FIXED["frequency_quadrature_point_number"]))
        fh.write("## Self-consistent RPA\n\n%s\n\n" % to_md(sc))
        fh.write("## Non-self-consistent RPA@DFT\n\n%s\n\n" % to_md(nsc))
        fh.write("`max |dE_eig| vs previous` is the worst-case eigenvalue change over "
                 "ALL occupied states between consecutive N_fe, not just the HOMO.\n")

    print("wrote %s" % os.path.join(OUT, "convergence.csv"))
    print("wrote %s" % os.path.join(OUT, "convergence.md"))
    print()
    print("=== self-consistent RPA ===")
    print(to_md(sc))
    print()
    print("=== non-self-consistent RPA@DFT ===")
    print(to_md(nsc))


if __name__ == "__main__":
    main()
