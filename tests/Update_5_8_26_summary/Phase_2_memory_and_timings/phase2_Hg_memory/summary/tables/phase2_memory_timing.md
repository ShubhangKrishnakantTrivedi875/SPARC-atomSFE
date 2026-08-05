# Phase 2 -- self-consistent RPA, peak memory and wall time

Mercury (Z=80), all-electron, one outer cycle.  Fixed: p=20, q=45, domain=13 Bohr, omega=4, polynomial mesh (concentration 2).  Self-consistent-potential RPA (xc='RPA'), corrected code vs original.

## Radial mesh sweep (Lmax = 30)

| N_fe | n_quad | new (MiB) | orig (MiB) | reduction | new (s) | orig (s) | speedup | dE (Ha) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 450 | 263 | - | - | 17.8 | - | - | - |
| 15 | 675 | 375 | - | - | 45.7 | - | - | - |
| 25 | 1125 | 650 | - | - | 149.4 | - | - | - |
| 40 | 1800 | 1390 | - | - | 460.8 | - | - | - |

## Angular sweep (N_fe = 40)

| Lmax | n_quad | new (MiB) | orig (MiB) | reduction | new (s) | orig (s) | speedup | dE (Ha) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 1800 | 1138 | - | - | 107.2 | - | - | - |
| 10 | 1800 | 1154 | - | - | 177.7 | - | - | - |
| 20 | 1800 | 1258 | - | - | 314.4 | - | - | - |
| 30 | 1800 | 1390 | - | - | 460.8 | - | - | - |

reduction = (new - orig)/orig on peak memory    speedup = orig/new on wall time    dE = new - orig on total energy
