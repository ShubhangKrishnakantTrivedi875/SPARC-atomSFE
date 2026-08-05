# Phase 2 -- non-self-consistent RPA, peak memory and wall time

Mercury (Z=80), all-electron, one outer cycle.  Fixed: p=20, q=45, domain=13 Bohr, omega=4, polynomial mesh (concentration 2).  E_c^RPA on frozen GGA_PBE orbitals (xc='RPA@DFT'), against self-consistent RPA in the corrected code.  No OEP equation is solved and no correlation energy density is built, so dE is physics, not error.

## Radial mesh sweep (Lmax = 30)

| N_fe | n_quad | RPA@DFT (MiB) | scRPA (MiB) | reduction | RPA@DFT (s) | scRPA (s) | speedup | dE (Ha) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 450 | 241 | 263 | -8.5% | 6.6 | 17.8 | 2.70x | -3.39e+00 |
| 15 | 675 | 333 | 375 | -11.1% | 13.2 | 45.7 | 3.46x | -3.39e+00 |
| 25 | 1125 | 586 | 650 | -9.9% | 40.5 | 149.4 | 3.69x | -3.39e+00 |
| 40 | 1800 | 1194 | 1390 | -14.1% | 129.3 | 460.8 | 3.56x | -3.39e+00 |

## Angular sweep (N_fe = 40)

| Lmax | n_quad | RPA@DFT (MiB) | scRPA (MiB) | reduction | RPA@DFT (s) | scRPA (s) | speedup | dE (Ha) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 1800 | 916 | 1138 | -19.5% | 35.0 | 107.2 | 3.06x | -3.39e+00 |
| 10 | 1800 | 967 | 1154 | -16.2% | 54.0 | 177.7 | 3.29x | -3.39e+00 |
| 20 | 1800 | 1094 | 1258 | -13.1% | 91.4 | 314.4 | 3.44x | -3.39e+00 |
| 30 | 1800 | 1194 | 1390 | -14.1% | 129.3 | 460.8 | 3.56x | -3.39e+00 |

reduction = (RPA@DFT - scRPA)/scRPA on peak memory    speedup = scRPA/RPA@DFT on wall time    dE = RPA@DFT - scRPA on total energy
