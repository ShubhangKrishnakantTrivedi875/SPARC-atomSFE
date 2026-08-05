"""
Per-functional comparison tables: new vs original vs OneAtomFEM.

    python make_tables.py

Writes to summary/tables/ :

    RPA_resources.csv / .md    peak memory and wall time per code, self-consistent RPA
    <XC>_energy.csv / .md      total energy per code, and the two absolute differences
    <XC>_eigenvalues.csv / .md HOMO eigenvalue per code, and the MAX absolute deviation
                               taken over ALL occupied states of that case
    comparison_tables.md       everything in one file, ready to paste

Six columns throughout, as requested:

    element | new | original | oneatom | |new-original| | |new-oneatom|

For the eigenvalue tables columns 2-4 are the HOMO, since a single number is
needed per code, while columns 5-6 are the worst-case deviation over every
occupied state -- the strictest statement available.

Each functional gets four blocks, one per (electron treatment, setting), because
those are not comparable with each other.
"""

from __future__ import annotations

import json
import os

import pandas as pd

import cases as C

OUT = os.path.join(C.HERE, "summary", "tables")
CODES = ["new", "orig", "oneatom"]

ELEMENT = {2: "He", 7: "N", 25: "Mn", 42: "Mo", 79: "Au"}


MIB = 1.048576          # psutil reports MB (10^6); tables quote MiB (2^20)


def load():
    """One row per (case, code) with the total energy, HOMO, eigenvalues and resources."""
    rows = []
    for c in C.ALL_CASES:
        p = os.path.join(c["outdir"], "result.json")
        if not os.path.exists(p):
            continue
        r = json.load(open(p))
        if not r.get("ok"):
            continue
        eig = r.get("eigenvalues") or []
        mem = r.get("memory") or {}
        peak = mem.get("peak_rss_tree_mb")
        rows.append(dict(code=c["code"], Z=c["Z"], electrons=c["electrons"],
                         xc=c["xc"], setting=c["setting"],
                         energy=r.get("total_energy"),
                         homo=(eig[-1] if eig else None),
                         eig=tuple(eig),
                         peak_mib=(peak / MIB) if peak else None,
                         wall_s=r.get("wall_s")))
    return pd.DataFrame(rows)


def six_col(df, xc, electrons, setting, quantity):
    """
    The six-column block for one (xc, electrons, setting).

    quantity = 'energy'      -> columns 2-4 are total energies
    quantity = 'eigenvalues' -> columns 2-4 are the HOMO, 5-6 the max deviation
                                over all occupied states
    """
    out = []
    for Z in C.Z_LIST:
        sel = {}
        for code in CODES:
            m = df[(df.code == code) & (df.Z == Z) & (df.electrons == electrons)
                   & (df.xc == xc) & (df.setting == setting)]
            sel[code] = m.iloc[0] if len(m) else None

        row = {"element": "%s (Z=%d)" % (ELEMENT.get(Z, "?"), Z)}
        if quantity == "energy":
            for code in CODES:
                row[code] = sel[code].energy if sel[code] is not None else None
            base = row["new"]
            for other, lab in (("orig", "|new-original|"), ("oneatom", "|new-oneatom|")):
                row[lab] = (abs(base - row[other])
                            if base is not None and row[other] is not None else None)
        else:
            for code in CODES:
                row[code] = sel[code].homo if sel[code] is not None else None
            ref = sel["new"]
            for other, lab in (("orig", "|new-original|"), ("oneatom", "|new-oneatom|")):
                o = sel[other]
                if ref is None or o is None or len(ref.eig) != len(o.eig) or not ref.eig:
                    row[lab] = None
                else:
                    row[lab] = max(abs(a - b) for a, b in zip(ref.eig, o.eig))
        out.append(row)
    cols = ["element"] + CODES + ["|new-original|", "|new-oneatom|"]
    return pd.DataFrame(out)[cols].rename(columns={"orig": "original"})


def resources_block(df, electrons, setting):
    """
    Peak memory and wall time for self-consistent RPA, one row per element.

    Peak memory is the sampled maximum resident set size of the whole run; for
    OneAtomFEM that is the sum over its 8 MPI ranks, so it is not a per-process
    figure.  Ratios are original/new, i.e. how much more the original needed.
    """
    out = []
    for Z in C.Z_LIST:
        row = {"element": "%s (Z=%d)" % (ELEMENT.get(Z, "?"), Z)}
        got = {}
        for code in CODES:
            m = df[(df.code == code) & (df.Z == Z) & (df.electrons == electrons)
                   & (df.xc == "RPA") & (df.setting == setting)]
            got[code] = m.iloc[0] if len(m) else None
            key = "original" if code == "orig" else code
            row["%s_MiB" % key] = got[code].peak_mib if got[code] is not None else None
            row["%s_s" % key]   = got[code].wall_s   if got[code] is not None else None
        if got["new"] is not None and got["orig"] is not None:
            row["mem_orig/new"]  = got["orig"].peak_mib / got["new"].peak_mib
            row["time_orig/new"] = got["orig"].wall_s / got["new"].wall_s
        else:
            row["mem_orig/new"] = row["time_orig/new"] = None
        out.append(row)
    cols = ["element", "new_MiB", "original_MiB", "oneatom_MiB",
            "new_s", "original_s", "oneatom_s", "mem_orig/new", "time_orig/new"]
    return pd.DataFrame(out)[cols]


def resources_md(d):
    head = ["element", "new (MiB)", "original (MiB)", "oneatom (MiB)",
            "new (s)", "original (s)", "oneatom (s)", "mem orig/new", "time orig/new"]
    keys = list(d.columns)
    lines = ["| " + " | ".join(head) + " |",
             "|:---|" + "|".join(["---:"] * (len(head) - 1)) + "|"]
    for _, r in d.iterrows():
        cells = [str(r["element"])]
        for k in keys[1:]:
            v = r[k]
            if pd.isna(v):
                cells.append("-")
            elif k.endswith("/new"):
                cells.append("%.1fx" % v)
            elif k.endswith("_MiB"):
                cells.append("%.0f" % v)
            else:
                cells.append("%.0f" % v)
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def to_md(d, value_fmt="%.8f", diff_fmt="%.2e"):
    head = ["element", "new", "original", "oneatom", "\\|new-original\\|", "\\|new-oneatom\\|"]
    keys = ["element", "new", "original", "oneatom", "|new-original|", "|new-oneatom|"]
    lines = ["| " + " | ".join(head) + " |",
             "|:---|" + "|".join(["---:"] * 5) + "|"]
    for _, r in d.iterrows():
        cells = [str(r["element"])]
        for k in keys[1:4]:
            cells.append("-" if pd.isna(r[k]) else value_fmt % r[k])
        for k in keys[4:]:
            cells.append("-" if pd.isna(r[k]) else diff_fmt % r[k])
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main():
    os.makedirs(OUT, exist_ok=True)
    df = load()
    if df.empty:
        print("no completed cases found")
        return

    everything = ["# Phase-1 comparison tables",
                  "",
                  "`new` = corrected SPARC-atomSFE, `original` = SPARC-atomSFE before the "
                  "corrections, `oneatom` = OneAtomFEM reference.",
                  "",
                  "Eigenvalue tables: columns 2-4 are the HOMO; columns 5-6 are the "
                  "**maximum** absolute deviation over **all** occupied states of that case.",
                  ""]

    for xc in ("LDA", "HF", "EXX", "RPA"):
        for quantity, unit in (("energy", "total energy, Ha"),
                               ("eigenvalues", "eigenvalues, Ha")):
            frames = []
            everything += ["", "## %s -- %s" % (xc, unit), ""]
            for electrons in ("ae", "psp"):
                for setting in ("loose", "tight"):
                    blk = six_col(df, xc, electrons, setting, quantity)
                    if blk[CODES[0]].isna().all():
                        continue
                    label = "%s, %s" % ("all-electron" if electrons == "ae"
                                        else "pseudopotential", setting)
                    everything += ["### %s" % label, "", to_md(blk), ""]
                    blk.insert(0, "setting", setting)
                    blk.insert(0, "electrons", electrons)
                    frames.append(blk)
            if frames:
                all_blk = pd.concat(frames, ignore_index=True)
                stem = os.path.join(OUT, "%s_%s" % (xc, quantity))
                all_blk.to_csv(stem + ".csv", index=False)
                with open(stem + ".md", "w") as fh:
                    fh.write("## %s -- %s\n\n" % (xc, unit))
                    for (el, st), g in all_blk.groupby(["electrons", "setting"], sort=False):
                        fh.write("### %s, %s\n\n%s\n\n"
                                 % ("all-electron" if el == "ae" else "pseudopotential",
                                    st, to_md(g.drop(columns=["electrons", "setting"]))))
                print("wrote %s.csv / .md" % stem)

    # ---- self-consistent RPA resources ----
    frames = []
    everything += ["", "## RPA -- peak memory and wall time", "",
                   "Peak memory is the sampled maximum resident set size over the whole "
                   "run.  For `oneatom` that is the sum across its 8 MPI ranks, so it is "
                   "not comparable per-process with the two thread-parallel codes.", ""]
    for electrons in ("ae", "psp"):
        for setting in ("loose", "tight"):
            blk = resources_block(df, electrons, setting)
            if blk["new_MiB"].isna().all():
                continue
            label = "%s, %s" % ("all-electron" if electrons == "ae"
                                else "pseudopotential", setting)
            everything += ["### %s" % label, "", resources_md(blk), ""]
            blk.insert(0, "setting", setting)
            blk.insert(0, "electrons", electrons)
            frames.append(blk)
    if frames:
        allr = pd.concat(frames, ignore_index=True)
        stem = os.path.join(OUT, "RPA_resources")
        allr.to_csv(stem + ".csv", index=False)
        with open(stem + ".md", "w") as fh:
            fh.write("## RPA -- peak memory and wall time\n\n")
            for (el, st), g in allr.groupby(["electrons", "setting"], sort=False):
                fh.write("### %s, %s\n\n%s\n\n"
                         % ("all-electron" if el == "ae" else "pseudopotential", st,
                            resources_md(g.drop(columns=["electrons", "setting"]))))
        print("wrote %s.csv / .md" % stem)

    with open(os.path.join(OUT, "comparison_tables.md"), "w") as fh:
        fh.write("\n".join(everything) + "\n")
    print("wrote %s" % os.path.join(OUT, "comparison_tables.md"))


if __name__ == "__main__":
    main()
