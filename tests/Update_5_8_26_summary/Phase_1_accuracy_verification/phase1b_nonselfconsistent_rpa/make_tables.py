"""
Phase-1b comparison tables: non-self-consistent RPA accuracy, new SPARC-atomSFE
vs OneAtomFEM_rpa_at_dft.

    python make_tables.py

Writes to summary/ :

    eigenvalues.csv / .md      HOMO + max |deviation| over ALL occupied states
    rpa_correlation.csv / .md  E_c^RPA component + |deviation| -- the actual
                               quantity under test
    total_energy.csv / .md     total energy + |deviation| -- see the caveat below,
                               this is NOT expected to agree
    comparison_tables.md       everything in one file, ready to paste

Only two codes here (no 'orig' -- the original SPARC-atomSFE has no 'RPA@DFT'),
so four columns throughout: element | new | oneatom | |new-oneatom|.  Eigenvalue
columns 2-3 are the HOMO; column 4 is the worst-case deviation over every
occupied state.

CAVEAT ON TOTAL ENERGY, load-bearing for reading this table:
SPARC's 'RPA@DFT' reports E_x^HF + E_c^RPA in place of E_x^GGA + E_c^GGA (see
cases.py's module docstring).  OneAtomFEM's emulation keeps E_x^GGA + E_c^GGA and
adds E_c^RPA on top (its 'Exact exchange (HF)' term is zeroed by alpha_x=0).  The
two codes are therefore computing two DIFFERENT, both internally-consistent
definitions of "RPA@DFT total energy" -- the |new-oneatom| column is expected to be
large and growing with Z, and is NOT a measure of RPA error.  Eigenvalues and
RPA_correlation are the columns that actually validate the RPA implementation,
since both codes converge the same GGA_PBE orbitals and evaluate a genuine
many-body RPA correlation energy on them.
"""

from __future__ import annotations

import json
import os

import pandas as pd

import cases as C

OUT = os.path.join(C.HERE, "summary")
CODES = ["new", "oneatom"]
ELEMENT = {2: "He", 7: "N", 25: "Mn", 42: "Mo", 79: "Au"}


def load():
    rows = []
    for c in C.ALL_CASES:
        p = os.path.join(c["outdir"], "result.json")
        if not os.path.exists(p):
            continue
        r = json.load(open(p))
        if not r.get("ok"):
            continue
        eig = r.get("eigenvalues") or []
        comps = r.get("energy_components") or {}
        rows.append(dict(code=c["code"], Z=c["Z"], setting=c["setting"],
                         energy=r.get("total_energy"),
                         homo=(eig[-1] if eig else None), eig=tuple(eig),
                         rpa_corr=comps.get("rpa_correlation")))
    return pd.DataFrame(rows)


def four_col(df, setting, quantity):
    out = []
    for Z in C.Z_LIST:
        sel = {}
        for code in CODES:
            m = df[(df.code == code) & (df.Z == Z) & (df.setting == setting)]
            sel[code] = m.iloc[0] if len(m) else None
        row = {"element": "%s (Z=%d)" % (ELEMENT.get(Z, "?"), Z)}
        if quantity == "energy":
            for code in CODES:
                row[code] = sel[code].energy if sel[code] is not None else None
            row["|new-oneatom|"] = (abs(row["new"] - row["oneatom"])
                                    if row["new"] is not None and row["oneatom"] is not None
                                    else None)
        elif quantity == "rpa_corr":
            for code in CODES:
                row[code] = sel[code].rpa_corr if sel[code] is not None else None
            row["|new-oneatom|"] = (abs(row["new"] - row["oneatom"])
                                    if row["new"] is not None and row["oneatom"] is not None
                                    else None)
        else:  # eigenvalues
            for code in CODES:
                row[code] = sel[code].homo if sel[code] is not None else None
            n, o = sel["new"], sel["oneatom"]
            if n is not None and o is not None and len(n.eig) == len(o.eig) and n.eig:
                row["max |dE| (occ)"] = max(abs(a - b) for a, b in zip(n.eig, o.eig))
            else:
                row["max |dE| (occ)"] = None
        out.append(row)
    cols = ["element", "new", "oneatom"] + \
           (["max |dE| (occ)"] if quantity == "eigenvalues" else ["|new-oneatom|"])
    return pd.DataFrame(out)[cols]


def to_md(d, value_fmt="%.9f", diff_fmt="%.2e"):
    keys = list(d.columns)
    lines = ["| " + " | ".join(keys) + " |",
             "|:---|" + "|".join(["---:"] * (len(keys) - 1)) + "|"]
    for _, r in d.iterrows():
        cells = [str(r["element"])]
        for k in keys[1:3]:
            v = r[k]
            cells.append("-" if pd.isna(v) else value_fmt % v)
        v = r[keys[3]]
        cells.append("-" if pd.isna(v) else diff_fmt % v)
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def main():
    os.makedirs(OUT, exist_ok=True)
    df = load()
    if df.empty:
        print("no completed cases found")
        return

    everything = ["# Phase-1b comparison tables",
                  "",
                  "Non-self-consistent RPA accuracy: `new` = corrected SPARC-atomSFE "
                  "(xc='RPA@DFT', ground_state_functional='GGA_PBE'), `oneatom` = "
                  "OneAtomFEM_rpa_at_dft emulation (xc='RPA', double_hybrid_flag=1, "
                  "alpha_x=alpha_c=0, OEP_method skipped).  All-electron only.",
                  "",
                  "**Total energy is not expected to agree** -- the two codes report "
                  "structurally different quantities (see this file's module "
                  "docstring).  Eigenvalues and RPA_correlation are the columns that "
                  "actually validate the RPA implementation.", ""]

    for quantity, unit, stem, vfmt in (
            ("eigenvalues", "HOMO eigenvalue, Ha", "eigenvalues", "%.9f"),
            ("rpa_corr", "E_c^RPA, Ha", "rpa_correlation", "%.9f"),
            ("energy", "total energy, Ha (see caveat)", "total_energy", "%.6f")):
        frames = []
        everything += ["", "## %s" % unit, ""]
        for setting in ("loose", "tight"):
            blk = four_col(df, setting, quantity)
            if blk["new"].isna().all():
                continue
            everything += ["### %s" % setting, "", to_md(blk, value_fmt=vfmt), ""]
            blk.insert(0, "setting", setting)
            frames.append(blk)
        if frames:
            all_blk = pd.concat(frames, ignore_index=True)
            path = os.path.join(OUT, stem)
            all_blk.to_csv(path + ".csv", index=False)
            with open(path + ".md", "w") as fh:
                fh.write("## %s\n\n" % unit)
                for setting, g in all_blk.groupby("setting", sort=False):
                    fh.write("### %s\n\n%s\n\n" %
                            (setting, to_md(g.drop(columns=["setting"]), value_fmt=vfmt)))
            print("wrote %s.csv / .md" % path)

    with open(os.path.join(OUT, "comparison_tables.md"), "w") as fh:
        fh.write("\n".join(everything) + "\n")
    print("wrote %s" % os.path.join(OUT, "comparison_tables.md"))


if __name__ == "__main__":
    main()
