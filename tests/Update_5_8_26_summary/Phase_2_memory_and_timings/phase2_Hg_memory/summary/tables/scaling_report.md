# Scaling report -- phase 2 (N_fe, Lmax) and phase 3 (quadrature order)

Mercury (Z=80), all-electron, corrected SPARC-atomSFE, full SCF (outer loop converged).  Fixed: p=20, q=45, domain=13 Bohr, omega=4, polynomial mesh (concentration 2).  Channels = Lmax + 4 (Hg has 4f occupied, l_occ_max = 3).

## Memory vs N_fe

Lmax held at 30.  x-axis is N_fe; the fit is against grid points (= N_fe x p - 1, the actual FE degrees of freedom -- p=20 here), since that is what the RPA kernel matrices scale with, not the quadrature order.  N_fe -> grid points: 10→199, 15→299, 25→499, 40→799.

## Wall time vs N_fe

Lmax held at 30.  x-axis is N_fe; the fit is against grid points (= N_fe x p - 1, the actual FE degrees of freedom -- p=20 here), since that is what the RPA kernel matrices scale with, not the quadrature order.  N_fe -> grid points: 10→199, 15→299, 25→499, 40→799.

## Memory vs Lmax

N_fe held at 40 (grid points = 799), so only the channel count (Lmax + 4) moves.  x-axis is Lmax.

## Wall time vs Lmax

N_fe held at 40 (grid points = 799), so only the channel count (Lmax + 4) moves.  x-axis is Lmax.

### N_fe sweep -- self-consistent RPA

| N_fe | grid points | outer | mem (MiB) | wall (s) | s/channel/outer | mem exp | time exp |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 199 | 5 | 295.1 | 55.1 | 0.324 | - | - |
| 15 | 299 | 5 | 450.0 | 165.4 | 0.973 | 1.04 | 2.7 |
| 25 | 499 | 5 | 737.7 | 584.0 | 3.44 | 0.965 | 2.46 |
| 40 | 799 | 5 | 1747.9 | 1769.3 | 10.4 | 1.83 | 2.35 |

### N_fe sweep -- non-self-consistent RPA@DFT

| N_fe | grid points | outer | mem (MiB) | wall (s) | s/channel/outer | mem exp | time exp |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 199 | 1 | 240.6 | 6.6 | 0.194 | - | - |
| 15 | 299 | 1 | 333.0 | 13.2 | 0.388 | 0.798 | 1.7 |
| 25 | 499 | 1 | 586.0 | 40.5 | 1.19 | 1.1 | 2.19 |
| 40 | 799 | 1 | 1193.9 | 129.3 | 3.8 | 1.51 | 2.47 |

### Lmax sweep -- self-consistent RPA

| Lmax | channels | outer | mem (MiB) | wall (s) | s/channel/outer | mem exp | time exp |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 9 | 4 | 1150.0 | 371.2 | 10.3 | - | - |
| 10 | 14 | 5 | 1188.9 | 741.6 | 10.6 | 0.0751 | 1.57 |
| 20 | 24 | 5 | 1402.3 | 1344.8 | 11.2 | 0.306 | 1.1 |
| 30 | 34 | 5 | 1747.9 | 1769.3 | 10.4 | 0.633 | 0.788 |

### Lmax sweep -- non-self-consistent RPA@DFT

| Lmax | channels | outer | mem (MiB) | wall (s) | s/channel/outer | mem exp | time exp |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 5 | 9 | 1 | 916.3 | 35.0 | 3.89 | - | - |
| 10 | 14 | 1 | 966.5 | 54.0 | 3.86 | 0.121 | 0.979 |
| 20 | 24 | 1 | 1093.9 | 91.4 | 3.81 | 0.23 | 0.977 |
| 30 | 34 | 1 | 1193.9 | 129.3 | 3.8 | 0.251 | 0.995 |

## Energy convergence, N_fe sweep -- self-consistent RPA

`dE vs previous` is this row's energy minus the row above it, within the same sweep (N_fe sweep at Lmax=30; Lmax sweep at N_fe=40) -- not a comparison against any predicted or extrapolated converged value.

| N_fe | Lmax | E (Ha) | dE vs previous (Ha) |
|---:|---:|---:|---:|
| 10 | 30 | -18413.639227975 | - |
| 15 | 30 | -18413.641135948 | -1.91e-03 |
| 25 | 30 | -18413.641200805 | -6.49e-05 |
| 40 | 30 | -18413.641206868 | -6.06e-06 |

## Energy convergence, Lmax sweep -- self-consistent RPA

`dE vs previous` is this row's energy minus the row above it, within the same sweep (N_fe sweep at Lmax=30; Lmax sweep at N_fe=40) -- not a comparison against any predicted or extrapolated converged value.

| N_fe | Lmax | E (Ha) | dE vs previous (Ha) |
|---:|---:|---:|---:|
| 40 | 5 | -18413.171562837 | - |
| 40 | 10 | -18413.619296343 | -4.48e-01 |
| 40 | 20 | -18413.640585551 | -2.13e-02 |
| 40 | 30 | -18413.641206868 | -6.21e-04 |

## Energy convergence, N_fe sweep -- non-self-consistent RPA@DFT

`dE vs previous` is this row's energy minus the row above it, within the same sweep (N_fe sweep at Lmax=30; Lmax sweep at N_fe=40) -- not a comparison against any predicted or extrapolated converged value.

| N_fe | Lmax | E (Ha) | dE vs previous (Ha) |
|---:|---:|---:|---:|
| 10 | 30 | -18413.631984602 | - |
| 15 | 30 | -18413.631983599 | +1.00e-06 |
| 25 | 30 | -18413.631981613 | +1.99e-06 |
| 40 | 30 | -18413.631969526 | +1.21e-05 |

## Energy convergence, Lmax sweep -- non-self-consistent RPA@DFT

`dE vs previous` is this row's energy minus the row above it, within the same sweep (N_fe sweep at Lmax=30; Lmax sweep at N_fe=40) -- not a comparison against any predicted or extrapolated converged value.

| N_fe | Lmax | E (Ha) | dE vs previous (Ha) |
|---:|---:|---:|---:|
| 40 | 5 | -18413.162880405 | - |
| 40 | 10 | -18413.610108664 | -4.47e-01 |
| 40 | 20 | -18413.631348768 | -2.12e-02 |
| 40 | 30 | -18413.631969526 | -6.21e-04 |

## Projections

`mem low` treats the Lmax term as additive and grid-point-independent; `mem high` lets it scale with the grid-points-squared term.  Only one grid-point count was sampled in the angular sweep, so which is correct is undetermined -- the true value lies between them.  This is the largest uncertainty in this report.  Wall time needs no such caveat: it is the measured outer-cycle count (scRPA: 5, from the N_fe=40/Lmax=30 anchor; RPA@DFT: 1, by construction) times the fitted per-channel cost.

| mode | N_fe | Lmax | grid points | channels | outer | mem low (GiB) | mem high (GiB) | wall (h) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| scRPA | 150 | 80 | 2999 | 84 | 5 | 22.5 | 37.9 | 33.8 |
| scRPA | 200 | 100 | 3999 | 104 | 5 | 39.4 | 79.2 | 85.6 |
| RPA@DFT | 150 | 80 | 2999 | 84 | 1 | 14.6 | 21.9 | 1.4 |
| RPA@DFT | 200 | 100 | 3999 | 104 | 1 | 25.7 | 44.2 | 3.2 |

## Cost ratio

Ratio is scRPA / RPA@DFT wall time, using each mode's own measured outer-cycle count -- no per-iteration normalisation is applied since neither number was ever per-iteration.  The ratio grows with grid size because scRPA's per-channel cost has the steeper exponent in grid points (2.49 vs 2.15).  Shaded rows are projected.

| N_fe | Lmax | grid points | scRPA (s) | RPA@DFT (s) | ratio | scRPA MiB | RPA@DFT MiB | mem ratio | kind |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 30 | 199 | 55.1 | 6.6 | 8.34x | 295.1 | 240.6 | 1.23x | measured |
| 15 | 30 | 299 | 165.4 | 13.2 | 12.53x | 450.0 | 333.0 | 1.35x | measured |
| 25 | 30 | 499 | 584.0 | 40.5 | 14.43x | 737.7 | 586.0 | 1.26x | measured |
| 40 | 5 | 799 | 371.2 | 35.0 | 10.60x | 1150.0 | 916.3 | 1.26x | measured |
| 40 | 10 | 799 | 741.6 | 54.0 | 13.74x | 1188.9 | 966.5 | 1.23x | measured |
| 40 | 20 | 799 | 1344.8 | 91.4 | 14.72x | 1402.3 | 1093.9 | 1.28x | measured |
| 40 | 30 | 799 | 1769.3 | 129.3 | 13.69x | 1747.9 | 1193.9 | 1.46x | measured |
| 150 | 80 | 2999 | 121618.0 | 5081.5 | 23.93x | 23013.4 | 14979.6 | 1.54x | projected |
| 200 | 100 | 3999 | 308232.0 | 11677.9 | 26.39x | 40308.5 | 26276.1 | 1.53x | projected |

## Phase 3 -- memory vs quadrature order q

Mercury (Z=80), all-electron, full SCF, N_fe=25, p=20, Lmax=30, omega=8, dense basis order = p.  n_quad = N_fe x q grows through q at FIXED N_fe here -- a different operation from phase 2's N_fe sweep, so the exponents below are not comparable to pages 1-2.

x-axis is n_quad = N_fe x q (N_fe=25 fixed).  scRPA: 4 outer cycles at every q; RPA@DFT: 1 (no outer loop).

## Phase 3 -- wall time vs quadrature order q

Mercury (Z=80), all-electron, full SCF, N_fe=25, p=20, Lmax=30, omega=8, dense basis order = p.  n_quad = N_fe x q grows through q at FIXED N_fe here -- a different operation from phase 2's N_fe sweep, so the exponents below are not comparable to pages 1-2.

x-axis is n_quad = N_fe x q (N_fe=25 fixed).  scRPA: 4 outer cycles at every q; RPA@DFT: 1 (no outer loop).

## Phase 3 data

Energy spread across the whole q sweep: scRPA 8.31e-09 Ha, RPA@DFT 2.25e-06 Ha -- q is converged at q=20 for both modes at this grid.

| mode | q | n_quad | outer | mem (MiB) | wall (s) | E (Ha) |
|---:|---:|---:|---:|---:|---:|---:|
| RPA@DFT | 20 | 500 | 1 | 360.7 | 13.6 | -18414.910538018 |
| RPA@DFT | 30 | 750 | 1 | 471.0 | 20.9 | -18414.910536919 |
| RPA@DFT | 40 | 1000 | 1 | 584.6 | 33.8 | -18414.910537648 |
| RPA@DFT | 50 | 1250 | 1 | 733.9 | 49.2 | -18414.910538091 |
| RPA@DFT | 60 | 1500 | 1 | 882.0 | 66.2 | -18414.910537448 |
| RPA@DFT | 80 | 2000 | 1 | 1250.6 | 116.7 | -18414.910535844 |
| scRPA | 20 | 500 | 4 | 582.9 | 164.4 | -18414.920617675 |
| scRPA | 30 | 750 | 4 | 646.2 | 280.7 | -18414.920617674 |
| scRPA | 40 | 1000 | 4 | 837.6 | 403.4 | -18414.920617679 |
| scRPA | 50 | 1250 | 4 | 1066.8 | 570.1 | -18414.920617671 |
| scRPA | 60 | 1500 | 4 | 1308.9 | 768.0 | -18414.920617671 |
| scRPA | 80 | 2000 | 4 | 1818.5 | 1362.6 | -18414.920617674 |

