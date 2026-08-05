# Snapshot: one-outer-cycle phase 2

Frozen copy of the phase-2 tables and the scaling report as they stood before the
full-SCF (`self_consistent_rpa_full`) runs were added.  Everything here was produced
from these two modes only:

    self_consistent_rpa      xc='RPA',     max_scf_iterations_outer = 1
    non_self_consistent_rpa  xc='RPA@DFT', GGA_PBE orbitals

so every self-consistent wall time in these files is the cost of ONE RPA evaluation,
not of a converged self-consistent RPA calculation.  Peak memory is unaffected by the
outer cap -- only one RPA build is ever live -- so the memory numbers and the memory
fits here remain the current best values.

Retained deliberately: the one-cycle timings are the cleanest measurement of the
per-evaluation cost, which is what the time exponents (n_quad^2.34 scRPA,
n_quad^2.16 RPA@DFT) and the channel-linearity check are fitted to.  The
regenerated top-level files add the full-SCF mode on top of this.
