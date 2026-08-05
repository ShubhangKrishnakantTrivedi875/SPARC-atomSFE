"""
Phase-2 study: RPA peak memory and wall time, corrected code vs original.

One system (mercury, all-electron, RPA) profiled along two axes:

    radial mesh   N_fe = 10, 15, 25, 40    at Lmax = 30
    angular reach Lmax =  5, 10, 20, 30    at N_fe = 40

so seven distinct configurations, sharing the (N_fe=40, Lmax=30) anchor.  Results
under the earlier (N_fe=25, Lmax=20) anchor are still on disk but are no longer
part of the study.  Everything else is held fixed:

    Z = 80 (Hg), all-electron, XC = RPA
    polynomial_order = 20, quadrature_point_number = 45, domain_size = 13 Bohr
    polynomial mesh, concentration 2, omega = 4

Three RPA modes are profiled over the same seven configurations:

    self_consistent_rpa       xc='RPA',     outer capped at 1 -- one RPA evaluation
    self_consistent_rpa_full  xc='RPA',     outer converged   -- the real cost
    non_self_consistent_rpa   xc='RPA@DFT', E_c^RPA on frozen GGA_PBE orbitals

The capped mode exists because it isolates the cost and the peak memory of a single
RPA evaluation; its total energy is not converged and must not be quoted.  The full
mode is what a real self-consistent RPA calculation costs, and dividing the two gives
the outer-iteration count directly.  Peak memory is the same in both, since only one
RPA build is ever live.

Only the corrected code implements 'RPA@DFT', and only it is worth running to outer
convergence, so both of those blocks are `new`-only: 14 capped + 7 full + 7
non-self-consistent runs.

Output tree:
    results/<mode>/<code>/Nfe<nn>/Lmax<nn>/

Cases are sorted mode-block first, then cheapest-first inside a block, so the
self-consistent indices are stable as modes are added.
"""

from __future__ import annotations

import os

HERE         = os.path.dirname(os.path.abspath(__file__))
RESULTS_ROOT = os.path.join(HERE, "results")

CODE_ROOT = {
    "new":  "/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE",
    "orig": "/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE_original/SPARC-atomSFE",
}
SPARC_CODES   = ("new", "orig")
ONEATOM_CODES = ()                 # not part of this study
USES_OEP      = {"RPA": True}      # kept for run_case.py compatibility

Z  = 80                            # mercury
XC = "RPA"

# `tag` names the results/ subtree, `codes` is which codes can run the mode, and the
# rest goes straight into the solver keywords.  'RPA@DFT' takes the
# orbitals from a converged GGA_PBE ground state and solves no OEP equation, hence
# use_oep=False -- the driver asserts exactly that.
MODES = {
    # one outer cycle: isolates the cost and the peak memory of a SINGLE RPA
    # evaluation.  Not a converged calculation -- its total energy is meaningless.
    "sc": dict(tag="self_consistent_rpa", codes=("new", "orig"),
               xc_functional="RPA", use_oep=True, ground_state_functional=None,
               outer_cycles=1),
    # outer loop left to converge: the real cost of a self-consistent RPA run.
    # Peak memory is unchanged (still one RPA build live at a time), so this mode
    # exists for the timing and the outer-iteration count.
    "scfull": dict(tag="self_consistent_rpa_full", codes=("new",),
                   xc_functional="RPA", use_oep=True, ground_state_functional=None,
                   outer_cycles=None),
    # no outer loop at all: E_c^RPA once, on frozen GGA_PBE orbitals
    "nsc": dict(tag="non_self_consistent_rpa", codes=("new",),
                xc_functional="RPA@DFT", use_oep=False,
                ground_state_functional="GGA_PBE", outer_cycles=1),
}
MODE_ORDER = ["sc", "scfull", "nsc"]

FIXED = dict(
    all_electron_flag        = True,
    polynomial_order         = 20,
    quadrature_point_number  = 45,
    domain_size              = 13.0,
    mesh_type                = "polynomial",
    mesh_concentration       = 2,
    scf_tolerance            = 1e-8,
    frequency_quadrature_point_number = 4,
)

#            (N_fe, Lmax)   -- the anchor (40, 30) belongs to both sweeps
MESH_SWEEP    = [(10, 30), (15, 30), (25, 30), (40, 30)]   # every N_fe at Lmax = 30
ANGULAR_SWEEP = [(40, 5), (40, 10), (40, 20), (40, 30)]    # every Lmax at N_fe = 40

CONFIGS = []
for n_fe, l_max in MESH_SWEEP + ANGULAR_SWEEP:
    if (n_fe, l_max) not in CONFIGS:
        CONFIGS.append((n_fe, l_max))


def cid(case) -> str:
    return "{mode}__{code}__Nfe{n_fe:02d}__Lmax{l_max:02d}".format(**case)


def outdir(case) -> str:
    return os.path.join(RESULTS_ROOT, MODES[case["mode"]]["tag"], case["code"],
                        "Nfe%02d" % case["n_fe"], "Lmax%02d" % case["l_max"])


def est_cost(case) -> float:
    n_q = case["n_fe"] * FIXED["quadrature_point_number"]
    return (n_q / 500.0) ** 3 * (case["l_max"] + 4)


def _build():
    out = []
    for mode in MODE_ORDER:
        block = []
        for code in MODES[mode]["codes"]:
            for n_fe, l_max in CONFIGS:
                c = dict(mode=mode, code=code, n_fe=n_fe, l_max=l_max, Z=Z,
                         xc=MODES[mode]["xc_functional"], electrons="ae",
                         setting="Nfe%02d_Lmax%02d" % (n_fe, l_max))
                c["cid"]    = cid(c)
                c["outdir"] = outdir(c)
                c["cost"]   = est_cost(c)
                block.append(c)
        # cheapest-first inside a mode block; blocks stay contiguous so adding a
        # mode never renumbers the ones before it
        block.sort(key=lambda c: (c["cost"], c["cid"]))
        out += block
    return out


ALL_CASES = _build()


def cases_for(codes) -> list:
    return [c for c in ALL_CASES if c["code"] in codes]


def sparc_cases()   -> list: return cases_for(SPARC_CODES)
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
    # ground_state_functional is a 'RPA@DFT'-only keyword and does not exist on the
    # original code at all, so it is dropped rather than passed as None
    return {k: v for k, v in p.items() if v is not None}


def indices_for(mode=None, codes=None) -> list:
    """Array indices matching a mode and/or a set of codes, for --array=..."""
    return [i for i, c in enumerate(ALL_CASES)
            if (mode is None or c["mode"] == mode)
            and (codes is None or c["code"] in codes)]


if __name__ == "__main__":
    print("# %d cases, mode-block then cheapest-first  (Z=%d, all-electron)"
          % (len(ALL_CASES), Z))
    print("# %-4s %-34s %-9s %6s %6s %8s %10s" %
          ("idx", "cid", "xc", "N_fe", "Lmax", "n_quad", "est_cost"))
    for i, c in enumerate(ALL_CASES):
        print("  %-4d %-34s %-9s %6d %6d %8d %10.1f"
              % (i, c["cid"], c["xc"], c["n_fe"], c["l_max"],
                 c["n_fe"] * FIXED["quadrature_point_number"], c["cost"]))
    for mode in MODE_ORDER:
        print("# %-24s --array=%s"
              % (MODES[mode]["tag"],
                 ",".join(str(i) for i in indices_for(mode, ("new",)))))
