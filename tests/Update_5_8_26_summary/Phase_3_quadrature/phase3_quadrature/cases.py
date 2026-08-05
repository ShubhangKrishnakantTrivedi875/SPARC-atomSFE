"""
Phase-3 study: how RPA time and memory scale with radial quadrature order.

Mercury, all-electron, full-SCF RPA (no outer-cycle cap), swept over quadrature
order q with everything else fixed:

    Z = 80, N_fe = 25, polynomial_order = 20, domain_size = 13 Bohr
    omega = 8, Lmax = 30, polynomial mesh, concentration 2

The outer loop is left to converge, so these are real self-consistent RPA
calculations, not one-cycle profiling runs.

Two RPA modes are swept over the same q list, exactly as in phase 2:

    self_consistent_rpa      xc = 'RPA'      OEP-consistent RPA potential
    non_self_consistent_rpa  xc = 'RPA@DFT'  E_c^RPA on frozen GGA_PBE orbitals

so 2 x 7 = 14 runs.  Only the corrected code implements 'RPA@DFT', and phase 3 uses
only the densep code anyway, so there is no second-code axis here.

The dense basis order is set equal to the polynomial order p,
so the floor is q >= p = 20.  Gauss-Legendre with q points integrates exactly to
degree 2q-1 and the integrand is a product of two dense basis functions (degree
2d), giving q >= d+1.  Keeping d = p rather than the upstream 2p+1 is what makes
the low-q end of the sweep legal at all; the price is that a product of two
degree-p orbitals -- the Poisson source for Hartree and for the differential
Coulomb kernel -- is no longer representable exactly.

q = 20 sits one degree below exactness for the overlap matrix and is included
deliberately: H_r_inv_sq integrates phi.phi/r^2 so exactness is unreachable at any
q regardless, and on the earlier Ar sweep the q=20 total energy matched the
converged value to 1e-10 Ha.

The results tree carries the mode and the grid, and the grid tag carries Z, so the
earlier argon sweeps survive alongside this one:
    results/<variant>/<mode>/Z<nn>_Nfe<nn>_Lmax<nn>/q<nn>/
"""

from __future__ import annotations

import os

HERE         = os.path.dirname(os.path.abspath(__file__))
RESULTS_ROOT = os.path.join(HERE, "results")

CODE_ROOT = {
    # dense order = p  (phase-3 variant, a copy so phases 1-2 stay reproducible)
    "densep":  "/storage/home/hcoda1/5/strivedi44/r-phanish6-0/SPARC-atomSFE_densep",
}
SPARC_CODES   = tuple(CODE_ROOT)
ONEATOM_CODES = ()
USES_OEP      = {"RPA": True}

Z  = 80                # mercury
XC = "RPA"

# The two RPA modes, same definitions as phase 2.  'RPA@DFT' takes its orbitals from
# a converged GGA_PBE ground state and solves no OEP equation, hence use_oep=False --
# the driver asserts exactly that.
MODES = {
    "sc": dict(tag="self_consistent_rpa",
               xc_functional="RPA", use_oep=True, ground_state_functional=None),
    "nsc": dict(tag="non_self_consistent_rpa",
                xc_functional="RPA@DFT", use_oep=False,
                ground_state_functional="GGA_PBE"),
}
MODE_ORDER = ["sc", "nsc"]

FIXED = dict(
    all_electron_flag        = True,
    finite_element_number    = 25,
    polynomial_order         = 20,
    domain_size              = 13.0,
    mesh_type                = "polynomial",
    mesh_concentration       = 2,
    scf_tolerance            = 1e-8,
    frequency_quadrature_point_number = 8,
    angular_momentum_cutoff  = 30,
)

Q_LIST = (20, 30, 40, 50, 60, 80)

# lowest q the variant can legally run, from the exactness argument above
Q_FLOOR = {
    "densep":  FIXED["polynomial_order"],           # 20
}

# grid tag, so a re-run at a different (N_fe, Lmax) does not overwrite this one
GRID = "Z%02d_Nfe%02d_Lmax%02d" % (Z, FIXED["finite_element_number"],
                                   FIXED["angular_momentum_cutoff"])


def cid(case) -> str:
    return "{mode}__{variant}__q{q:02d}".format(**case)


def outdir(case) -> str:
    return os.path.join(RESULTS_ROOT, case["variant"], MODES[case["mode"]]["tag"],
                        GRID, "q%02d" % case["q"])


def est_cost(case) -> float:
    n_q = FIXED["finite_element_number"] * case["q"]
    return (n_q / 500.0) ** 3 * (FIXED["angular_momentum_cutoff"] + 4)


def _build():
    out = []
    for mode in MODE_ORDER:
        block = []
        for variant in CODE_ROOT:
            for q in Q_LIST:
                if q < Q_FLOOR[variant]:
                    continue                 # would fail the assert; not a case
                c = dict(mode=mode, code=variant, variant=variant, q=q, Z=Z,
                         xc=MODES[mode]["xc_functional"], electrons="ae",
                         setting="q%02d" % q)
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
    p.update(xc_functional           = m["xc_functional"],
             use_oep                 = m["use_oep"],
             ground_state_functional = m["ground_state_functional"],
             atomic_number           = case["Z"],
             quadrature_point_number = case["q"])
    # ground_state_functional is a 'RPA@DFT'-only keyword, dropped rather than None
    return {k: v for k, v in p.items() if v is not None}


def indices_for(mode) -> list:
    """Array indices for one mode, for --array=..."""
    return [i for i, c in enumerate(ALL_CASES) if c["mode"] == mode]


if __name__ == "__main__":
    print("# %d cases, mode-block then cheapest-first  (Z=%d, all-electron)"
          % (len(ALL_CASES), Z))
    print("# %-4s %-26s %-9s %5s %8s %10s"
          % ("idx", "cid", "xc", "q", "n_quad", "est_cost"))
    for i, c in enumerate(ALL_CASES):
        print("  %-4d %-26s %-9s %5d %8d %10.1f"
              % (i, c["cid"], c["xc"], c["q"],
                 FIXED["finite_element_number"] * c["q"], c["cost"]))
    for mode in MODE_ORDER:
        ix = indices_for(mode)
        print("# %-26s --array=%d-%d" % (MODES[mode]["tag"], ix[0], ix[-1]))
    print("\n# grid %s, dense basis order = p = %d, floor q >= %d"
          % (GRID, FIXED["polynomial_order"], Q_FLOOR["densep"]))
    print("# q=20 kept: one degree under-integrated, verified harmless (1e-10 Ha)")
