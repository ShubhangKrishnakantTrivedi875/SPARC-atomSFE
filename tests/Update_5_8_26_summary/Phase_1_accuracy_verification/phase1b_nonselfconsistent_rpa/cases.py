"""
Phase-1b: non-self-consistent RPA accuracy, new SPARC-atomSFE vs OneAtomFEM.

Same 5 systems and same two settings as phase1_validation, all-electron only:

    xc = 'RPA@DFT', ground_state_functional = 'GGA_PBE'   (new SPARC-atomSFE)
    xc = 'RPA', double_hybrid_flag=1, alpha_x=0, alpha_c=0,
        OEP_method call skipped                            (OneAtomFEM_rpa_at_dft)

OneAtomFEM has no 'RPA@DFT' functional; OneAtomFEM_rpa_at_dft is a dedicated copy
that emulates it by hardcoding double_hybrid_flag/alpha_x/alpha_c and skipping the
(otherwise needless) OEP_method call -- see that copy's module docstring diff.

Both codes then converge plain GGA_PBE orbitals and evaluate the RPA correlation
energy once on them, so eigenvalues should agree closely between codes.  Total
energy is NOT expected to agree: SPARC's 'RPA@DFT' reports E_x^HF + E_c^RPA in
place of E_x^GGA + E_c^GGA, while OneAtomFEM's emulation keeps E_x^GGA + E_c^GGA
and adds E_c^RPA on top (its 'Exact exchange (HF)' term is zeroed by alpha_x=0).
Compare eigenvalues and the RPA_correlation component directly; total energy carries
this known, structural offset on the OneAtomFEM side.

Output tree:
    results/<code>/Z<zz>/<setting>/
"""

from __future__ import annotations

import os

HERE         = os.path.dirname(os.path.abspath(__file__))
RESULTS_ROOT = os.path.join(HERE, "results")

CODE_ROOT = {
    "new":     "/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE",
    "oneatom": "/storage/home/hcoda1/5/strivedi44/r-phanish6-0/OneAtomFEM_rpa_at_dft",
}
SPARC_CODES   = ("new",)
ONEATOM_CODES = ("oneatom",)

XC_NAME  = {"RPA@DFT": {"new": "RPA@DFT", "oneatom": "RPA"}}
XC_LIST  = ("RPA@DFT",)
USES_OEP = {"RPA@DFT": False}          # 'new' solves no OEP equation; harness compat
USES_RPA_GRID = {"RPA@DFT": True}      # omega / Lmax still apply (RPA correlation)

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

ELECTRONS = ("ae",)                    # all-electron only, per this study


def est_cost(case) -> float:
    s   = SETTINGS[case["setting"]]
    n_q = s["finite_element_number"] * s["quadrature_point_number"]
    c   = (n_q / 500.0) ** 3 * s["frequency_quadrature_point_number"] * \
          (s["angular_momentum_cutoff"] + 4)
    c *= 1.0 + case["Z"] / 20.0
    return c


def cid(case) -> str:
    return "{code}__Z{Z:02d}__ae__RPA@DFT__{setting}".format(**case)


def outdir(case) -> str:
    return os.path.join(RESULTS_ROOT, case["code"], "Z%02d" % case["Z"],
                        case["setting"])


def params_for(case) -> dict:
    p = dict(COMMON)
    p.update(SETTINGS[case["setting"]])
    p.update(
        xc_functional     = XC_NAME["RPA@DFT"][case["code"]],
        all_electron_flag = True,
        atomic_number     = case["Z"],
    )
    if case["code"] == "new":
        p["ground_state_functional"] = "GGA_PBE"
        p["use_oep"] = False
    return p


def _build():
    out = []
    for code in CODE_ROOT:
        for Z in Z_LIST:
            for setting in SETTINGS:
                c = dict(code=code, Z=Z, xc="RPA@DFT", electrons="ae", setting=setting)
                c["cid"]    = cid(c)
                c["outdir"] = outdir(c)
                c["cost"]   = est_cost(c)
                out.append(c)
    out.sort(key=lambda c: (c["cost"], c["cid"]))
    return out


ALL_CASES = _build()


def cases_for(codes) -> list:
    return [c for c in ALL_CASES if c["code"] in codes]


def sparc_cases()   -> list: return cases_for(SPARC_CODES)
def oneatom_cases() -> list: return cases_for(ONEATOM_CODES)


def by_cid(cid_str):
    for c in ALL_CASES:
        if c["cid"] == cid_str:
            return c
    raise KeyError(cid_str)


if __name__ == "__main__":
    print("# %d cases, cheapest first  (all-electron, RPA@DFT / GGA_PBE)"
          % len(ALL_CASES))
    print("# %-4s %-30s %6s %8s" % ("idx", "cid", "Z", "est_cost"))
    for i, c in enumerate(ALL_CASES):
        print("  %-4d %-30s %6d %8.1f" % (i, c["cid"], c["Z"], c["cost"]))
