"""
Phase-1 validation matrix.

Three codes x five atoms x two parameter settings x four functionals x {AE, PSP}.
Each case carries a stable cid and a hierarchical output directory

    results/<code>/Z<zz>/<ae|psp>/<xc>/<setting>/

so cases can be added or re-run individually without disturbing the rest.

The two SPARC codes share one API and are driven in-process by a child script.
OneAtomFEM is driven through environment overrides and mpirun, and its
functional names differ -- see XC_NAME.
"""

from __future__ import annotations

import os

HERE         = os.path.dirname(os.path.abspath(__file__))
RESULTS_ROOT = os.path.join(HERE, "results")
LOG_ROOT     = os.path.join(HERE, "logs")

# ---------------------------------------------------------------------------
#  Codes under test
# ---------------------------------------------------------------------------
CODE_ROOT = {
    "new":     "/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE",
    "orig":    "/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE_original/SPARC-atomSFE",
    "oneatom": "/storage/home/hcoda1/5/strivedi44/r-phanish6-0/OneAtomFEM",
}
SPARC_CODES   = ("new", "orig")      # same API, run through _child_sparc.py
ONEATOM_CODES = ("oneatom",)         # env overrides + mpirun

# ---------------------------------------------------------------------------
#  Functional names.  'LDA' is PW92 in both codes; the OEP-based ones are
#  spelled differently.
# ---------------------------------------------------------------------------
XC_NAME = {
    "LDA": {"new": "LDA_PW", "orig": "LDA_PW", "oneatom": "LDA_SPW"},
    "HF":  {"new": "HF",     "orig": "HF",     "oneatom": "HFx"},
    "EXX": {"new": "EXX",    "orig": "EXX",    "oneatom": "OEPx"},
    "RPA": {"new": "RPA",    "orig": "RPA",    "oneatom": "RPA"},
}
XC_LIST  = ("LDA", "HF", "EXX", "RPA")
USES_OEP = {"LDA": False, "HF": False, "EXX": True, "RPA": True}
# omega / Lmax only mean anything for the OEP-based functionals
USES_RPA_GRID = {"LDA": False, "HF": False, "EXX": True, "RPA": True}

# ---------------------------------------------------------------------------
#  Systems and parameter settings
# ---------------------------------------------------------------------------
Z_LIST = (2, 7, 25, 42, 79)

SETTINGS = {
    "loose": dict(finite_element_number=15, polynomial_order=20,
                  quadrature_point_number=45, domain_size=13.0,
                  frequency_quadrature_point_number=8, angular_momentum_cutoff=4),
    "tight": dict(finite_element_number=30, polynomial_order=20,
                  quadrature_point_number=45, domain_size=13.0,
                  frequency_quadrature_point_number=8, angular_momentum_cutoff=25),
}

COMMON = dict(
    mesh_type          = "polynomial",
    mesh_concentration = 2,
    scf_tolerance      = 1e-8,
)

ELECTRONS = ("ae", "psp")            # all-electron / pseudopotential


# ---------------------------------------------------------------------------
#  Cost model -- only used to order the array so cheap cases land first and to
#  pick a per-case time limit.  Rough, deliberately.
# ---------------------------------------------------------------------------
def est_cost(case) -> float:
    s   = SETTINGS[case["setting"]]
    n_q = s["finite_element_number"] * s["quadrature_point_number"]
    c   = (n_q / 500.0) ** 3                                  # dense linear algebra
    if case["xc"] == "RPA":
        c *= s["frequency_quadrature_point_number"] * (s["angular_momentum_cutoff"] + 4)
    elif case["xc"] == "EXX":
        c *= 8
    elif case["xc"] == "HF":
        c *= 4
    c *= 1.0 + case["Z"] / 20.0                               # more occupied states
    return c


def outdir(case) -> str:
    return os.path.join(RESULTS_ROOT, case["code"], "Z%02d" % case["Z"],
                        case["electrons"], case["xc"], case["setting"])


def cid(case) -> str:
    return "{code}__Z{Z:02d}__{electrons}__{xc}__{setting}".format(**case)


def _build():
    out = []
    for code in CODE_ROOT:
        for Z in Z_LIST:
            for electrons in ELECTRONS:
                for xc in XC_LIST:
                    for setting in SETTINGS:
                        c = dict(code=code, Z=Z, electrons=electrons,
                                 xc=xc, setting=setting)
                        c["cid"]    = cid(c)
                        c["outdir"] = outdir(c)
                        c["cost"]   = est_cost(c)
                        out.append(c)
    return out


ALL_CASES = _build()


def cases_for(codes) -> list:
    """Cases for the given code names, cheapest first (stable ordering)."""
    sel = [c for c in ALL_CASES if c["code"] in codes]
    sel.sort(key=lambda c: (c["cost"], c["cid"]))
    return sel


def sparc_cases()   -> list: return cases_for(SPARC_CODES)
def oneatom_cases() -> list: return cases_for(ONEATOM_CODES)


def by_cid(cid_str):
    for c in ALL_CASES:
        if c["cid"] == cid_str:
            return c
    raise KeyError(cid_str)


def params_for(case) -> dict:
    """Flat parameter dict for one case, in SPARC-atomSFE naming."""
    s = dict(SETTINGS[case["setting"]])
    p = dict(COMMON)
    p.update(
        atomic_number           = case["Z"],
        xc_functional           = XC_NAME[case["xc"]][case["code"]],
        all_electron_flag       = (case["electrons"] == "ae"),
        use_oep                 = USES_OEP[case["xc"]],
        finite_element_number   = s["finite_element_number"],
        polynomial_order        = s["polynomial_order"],
        quadrature_point_number = s["quadrature_point_number"],
        domain_size             = s["domain_size"],
    )
    if USES_RPA_GRID[case["xc"]]:
        p["frequency_quadrature_point_number"] = s["frequency_quadrature_point_number"]
        p["angular_momentum_cutoff"]           = s["angular_momentum_cutoff"]
    return p


if __name__ == "__main__":
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else "all"
    sel = {"sparc": sparc_cases, "oneatom": oneatom_cases,
           "all": lambda: cases_for(CODE_ROOT)}[which]()
    print("# %d cases (%s), cheapest first" % (len(sel), which))
    print("# %-5s %-56s %10s" % ("idx", "cid", "est_cost"))
    for i, c in enumerate(sel):
        print("  %-5d %-56s %10.1f" % (i, c["cid"], c["cost"]))
