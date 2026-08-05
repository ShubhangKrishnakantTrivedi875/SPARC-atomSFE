"""
Beryllium N_fe convergence: RPA peak memory, wall time, and total energy vs
radial mesh, self-consistent RPA and RPA@DFT.

Same fixed parameters as the Mercury phase-2 study (p=20, q=45, domain=13 Bohr,
polynomial mesh concentration 2, omega=24), Lmax held at 30, all-electron.  Only
N_fe is swept: 5, 10, 15, 20, 25, 30.  No Lmax sweep and no capped 1-outer-cycle
profiling mode here -- this is a convergence study, so self-consistent RPA is run
to full outer-loop convergence (outer_cycles left uncapped), matching the
'self_consistent_rpa_full' mode from the Mercury study, just not sharing its name
since there is no capped variant alongside it in this study.

Two RPA modes, corrected code only (no 'orig' or 'oneatom' comparison here):

    sc   self-consistent RPA       xc='RPA',     full outer-loop convergence
    nsc  non-self-consistent RPA@DFT  xc='RPA@DFT', GGA_PBE orbitals

Output tree:
    results/<mode_tag>/Nfe<nn>/
"""

from __future__ import annotations

import os

HERE         = os.path.dirname(os.path.abspath(__file__))
RESULTS_ROOT = os.path.join(HERE, "results")

CODE_ROOT = {
    "new": "/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE",
}
SPARC_CODES   = ("new",)
ONEATOM_CODES = ()
USES_OEP = {"RPA": True, "RPA@DFT": False}      # required by run_case.py's n_ranks logic

Z  = 4                              # beryllium
L_MAX = 30

MODES = {
    "sc": dict(tag="self_consistent_rpa",
               xc_functional="RPA", use_oep=True, ground_state_functional=None,
               outer_cycles=None),  # uncapped -- converges the outer loop
    "nsc": dict(tag="non_self_consistent_rpa",
                xc_functional="RPA@DFT", use_oep=False,
                ground_state_functional="GGA_PBE", outer_cycles=1),
}
MODE_ORDER = ["sc", "nsc"]

FIXED = dict(
    all_electron_flag        = True,
    polynomial_order          = 20,
    quadrature_point_number   = 45,
    domain_size                = 13.0,
    mesh_type                  = "polynomial",
    mesh_concentration         = 2,
    scf_tolerance               = 1e-8,
    frequency_quadrature_point_number = 24,
)

N_FE_SWEEP = [5, 10, 15, 20, 25, 30]


def cid(case) -> str:
    return "{mode}__Nfe{n_fe:02d}".format(**case)


def outdir(case) -> str:
    return os.path.join(RESULTS_ROOT, MODES[case["mode"]]["tag"], "Nfe%02d" % case["n_fe"])


def est_cost(case) -> float:
    n_q = case["n_fe"] * FIXED["quadrature_point_number"]
    return (n_q / 500.0) ** 3 * (L_MAX + 4)


def _build():
    out = []
    for mode in MODE_ORDER:
        block = []
        for n_fe in N_FE_SWEEP:
            c = dict(mode=mode, code="new", n_fe=n_fe, l_max=L_MAX, Z=Z,
                     xc=MODES[mode]["xc_functional"], electrons="ae",
                     setting="Nfe%02d" % n_fe)
            c["cid"]    = cid(c)
            c["outdir"] = outdir(c)
            c["cost"]   = est_cost(c)
            block.append(c)
        block.sort(key=lambda c: (c["cost"], c["cid"]))
        out += block
    return out


ALL_CASES = _build()


def cases_for(codes) -> list:
    return [c for c in ALL_CASES if c["code"] in codes]


def sparc_cases()   -> list: return cases_for(("new",))
def oneatom_cases() -> list: return []


def by_cid(cid_str):
    for c in ALL_CASES:
        if c["cid"] == cid_str:
            return c
    raise KeyError(cid_str)


def params_for(case) -> dict:
    m = MODES[case["mode"]]
    p = dict(FIXED)
    p.update(
        xc_functional            = m["xc_functional"],
        use_oep                  = m["use_oep"],
        ground_state_functional  = m["ground_state_functional"],
        atomic_number            = case["Z"],
        finite_element_number    = case["n_fe"],
        angular_momentum_cutoff  = case["l_max"],
        max_scf_iterations_outer = m["outer_cycles"],
    )
    return {k: v for k, v in p.items() if v is not None}


def indices_for(mode) -> list:
    return [i for i, c in enumerate(ALL_CASES) if c["mode"] == mode]


if __name__ == "__main__":
    print("# %d cases, mode-block then cheapest-first  (Z=%d, all-electron, Lmax=%d)"
          % (len(ALL_CASES), Z, L_MAX))
    print("# %-4s %-28s %-9s %6s %8s %10s" %
          ("idx", "cid", "xc", "N_fe", "n_quad", "est_cost"))
    for i, c in enumerate(ALL_CASES):
        print("  %-4d %-28s %-9s %6d %8d %10.1f"
              % (i, c["cid"], c["xc"], c["n_fe"],
                 c["n_fe"] * FIXED["quadrature_point_number"], c["cost"]))
    for mode in MODE_ORDER:
        ix = indices_for(mode)
        print("# %-24s --array=%d-%d" % (MODES[mode]["tag"], ix[0], ix[-1]))
