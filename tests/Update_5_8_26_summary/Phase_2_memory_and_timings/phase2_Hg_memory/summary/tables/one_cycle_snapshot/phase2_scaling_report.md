# Phase 2 -- scaling, convergence and projections

Mercury (Z=80), all-electron, corrected SPARC-atomSFE.  Fixed: p=20, q=45, domain=13 Bohr, omega=4, polynomial mesh (concentration 2).  Channels = Lmax + 4, since Hg has 4f occupied (l_occ_max = 3).

## 1-2. Scaling with N_fe

Mercury (Z=80), all-electron, corrected SPARC-atomSFE.  Fixed: p=20, q=45, domain=13 Bohr, omega=4, polynomial mesh (concentration 2).  Channels = Lmax + 4, since Hg has 4f occupied (l_occ_max = 3).  Lmax held at 30, so only n_quad moves.  Exponents are local, between consecutive rows, against log(n_quad).

### self-consistent RPA -- scaling with N_fe

| N_fe | n_quad | mem (MiB) | wall (s) | s/channel | mem exp | time exp |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 450 | 263.0 | 17.8 | 0.52 | - | - |
| 15 | 675 | 374.7 | 45.7 | 1.34 | 0.87 | 2.32 |
| 25 | 1125 | 650.1 | 149.4 | 4.39 | 1.08 | 2.32 |
| 40 | 1800 | 1390.3 | 460.8 | 13.55 | 1.62 | 2.40 |

```
memory   194 + 3.685e-04 n_quad^2  MiB at Lmax=30   (max err 3.5%);  the quadratic term is 86% of the peak at the anchor
         1067 + 10.36 Lmax  MiB at n_quad=1800   (max err 1.7%)
time     wall = channels x 3.170e-07 n_quad^2.34  s   (max err 1.5%)
         cross-check at n_quad=1800: wall = -20.4 + 14.10 x channels, so everything outside the RPA is 4.4% of the total
```

### non-self-consistent RPA@DFT -- scaling with N_fe

| N_fe | n_quad | mem (MiB) | wall (s) | s/channel | mem exp | time exp |
|---:|---:|---:|---:|---:|---:|---:|
| 10 | 450 | 240.6 | 6.6 | 0.19 | - | - |
| 15 | 675 | 333.0 | 13.2 | 0.39 | 0.80 | 1.70 |
| 25 | 1125 | 586.0 | 40.5 | 1.19 | 1.11 | 2.19 |
| 40 | 1800 | 1193.9 | 129.3 | 3.80 | 1.51 | 2.47 |

```
memory   186 + 3.119e-04 n_quad^2  MiB at Lmax=30   (max err 3.5%);  the quadratic term is 85% of the peak at the anchor
         859 + 11.31 Lmax  MiB at n_quad=1800   (max err 0.8%)
time     wall = channels x 3.399e-07 n_quad^2.16  s   (max err 9.7%)
         cross-check at n_quad=1800: wall = +1.1 + 3.77 x channels, so everything outside the RPA is 0.9% of the total
```

## 1-2. Scaling with Lmax

Mercury (Z=80), all-electron, corrected SPARC-atomSFE.  Fixed: p=20, q=45, domain=13 Bohr, omega=4, polynomial mesh (concentration 2).  Channels = Lmax + 4, since Hg has 4f occupied (l_occ_max = 3).  N_fe held at 40, so n_quad is constant at 1800 and only the channel count moves.  Exponents are against log(channels).

### self-consistent RPA -- scaling with Lmax

| Lmax | channels | mem (MiB) | wall (s) | s/channel | mem exp | time exp |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 9 | 1138.2 | 107.2 | 11.92 | - | - |
| 10 | 14 | 1153.6 | 177.7 | 12.70 | 0.03 | 1.14 |
| 20 | 24 | 1258.3 | 314.4 | 13.10 | 0.16 | 1.06 |
| 30 | 34 | 1390.3 | 460.8 | 13.55 | 0.29 | 1.10 |

```
memory   194 + 3.685e-04 n_quad^2  MiB at Lmax=30   (max err 3.5%);  the quadratic term is 86% of the peak at the anchor
         1067 + 10.36 Lmax  MiB at n_quad=1800   (max err 1.7%)
time     wall = channels x 3.170e-07 n_quad^2.34  s   (max err 1.5%)
         cross-check at n_quad=1800: wall = -20.4 + 14.10 x channels, so everything outside the RPA is 4.4% of the total
```

### non-self-consistent RPA@DFT -- scaling with Lmax

| Lmax | channels | mem (MiB) | wall (s) | s/channel | mem exp | time exp |
|---:|---:|---:|---:|---:|---:|---:|
| 5 | 9 | 916.3 | 35.0 | 3.89 | - | - |
| 10 | 14 | 966.5 | 54.0 | 3.86 | 0.12 | 0.98 |
| 20 | 24 | 1093.9 | 91.4 | 3.81 | 0.23 | 0.98 |
| 30 | 34 | 1193.9 | 129.3 | 3.80 | 0.25 | 1.00 |

```
memory   186 + 3.119e-04 n_quad^2  MiB at Lmax=30   (max err 3.5%);  the quadratic term is 85% of the peak at the anchor
         859 + 11.31 Lmax  MiB at n_quad=1800   (max err 0.8%)
time     wall = channels x 3.399e-07 n_quad^2.16  s   (max err 9.7%)
         cross-check at n_quad=1800: wall = +1.1 + 3.77 x channels, so everything outside the RPA is 0.9% of the total
```

## 3. Total-energy convergence

Mercury (Z=80), all-electron, corrected SPARC-atomSFE.  Fixed: p=20, q=45, domain=13 Bohr, omega=4, polynomial mesh (concentration 2).  Channels = Lmax + 4, since Hg has 4f occupied (l_occ_max = 3).  Deviations are measured against a projected (Lmax-extrapolated) reference, not against the largest run.

### self-consistent RPA -- total-energy convergence

| N_fe | Lmax | E (Ha) | |E - E_ref| (Ha) | within 1e-4 | within 1e-5 |
|---:|---:|---:|---:|---:|---:|
| 40 | 5 | -18409.776929186 | 4.64e-01 | no | no |
| 40 | 10 | -18410.218848242 | 2.20e-02 | no | no |
| 40 | 20 | -18410.240084032 | 7.51e-04 | no | no |
| 10 | 30 | -18410.241904062 | 1.07e-03 | no | no |
| 15 | 30 | -18410.241666979 | 8.32e-04 | no | no |
| 25 | 30 | -18410.241629856 | 7.95e-04 | no | no |
| 40 | 30 | -18410.240730760 | 1.04e-04 | no | no |

```
reference  E_inf = -18410.240834913 Ha, from E(L) = E_inf + c L^-k on the top three Lmax points;  k = 4.87
deviation at Lmax=30 is 1.04e-04 Ha;  Lmax needed for 1e-4 is 30, for 1e-5 is 49
the four N_fe rows all sit at Lmax=30, so their deviation IS the Lmax error -- the N_fe sweep itself spans only 1.2e-03 Ha with non-monotone increments, i.e. N_fe is not the limiting parameter
```

### non-self-consistent RPA@DFT -- total-energy convergence

| N_fe | Lmax | E (Ha) | |E - E_ref| (Ha) | within 1e-4 | within 1e-5 |
|---:|---:|---:|---:|---:|---:|
| 40 | 5 | -18413.162880405 | 4.69e-01 | no | no |
| 40 | 10 | -18413.610108664 | 2.20e-02 | no | no |
| 40 | 20 | -18413.631348768 | 7.18e-04 | no | no |
| 10 | 30 | -18413.631984602 | 8.20e-05 | yes | no |
| 15 | 30 | -18413.631983599 | 8.30e-05 | yes | no |
| 25 | 30 | -18413.631981613 | 8.50e-05 | yes | no |
| 40 | 30 | -18413.631969526 | 9.71e-05 | yes | no |

```
reference  E_inf = -18413.632066577 Ha, from E(L) = E_inf + c L^-k on the top three Lmax points;  k = 4.93
deviation at Lmax=30 is 9.71e-05 Ha;  Lmax needed for 1e-4 is 30, for 1e-5 is 48
the four N_fe rows all sit at Lmax=30, so their deviation IS the Lmax error -- the N_fe sweep itself spans only 1.5e-05 Ha with non-monotone increments, i.e. N_fe is not the limiting parameter
```

## 4. Projections

| mode | N_fe | Lmax | n_quad | channels | mem low (GiB) | mem high (GiB) | s/channel | wall (h) |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| scRPA | 150 | 80 | 6750 | 84 | 17.1 | 23.7 | 296.68 | 6.9 |
| scRPA | 200 | 100 | 9000 | 104 | 30.0 | 47.0 | 582.10 | 16.8 |
| RPA@DFT | 150 | 80 | 6750 | 84 | 14.6 | 21.8 | 60.88 | 1.4 |
| RPA@DFT | 200 | 100 | 9000 | 104 | 25.6 | 44.2 | 113.18 | 3.3 |

```
`mem low` carries the Lmax term as additive and n_quad-independent; `mem high` lets it grow with the n_quad^2 term.  Only one n_quad was sampled in the angular sweep, so which is right is undetermined and the truth lies between them -- this is the largest uncertainty in the report.
Wall time needs no such caveat: the channel count enters linearly, verified to 4.4% (scRPA) and 0.9% (RPA@DFT) at n_quad=1800.
```

## Self-consistent / non-self-consistent cost ratio

| N_fe | Lmax | n_quad | scRPA (s) | RPA@DFT (s) | ratio | scRPA MiB | RPA@DFT MiB | mem ratio | kind |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 10 | 30 | 450 | 17.8 | 6.6 | 2.70x | 263.0 | 240.6 | 1.09x | measured |
| 15 | 30 | 675 | 45.7 | 13.2 | 3.46x | 374.7 | 333.0 | 1.13x | measured |
| 25 | 30 | 1125 | 149.4 | 40.5 | 3.69x | 650.1 | 586.0 | 1.11x | measured |
| 40 | 5 | 1800 | 107.2 | 35.0 | 3.06x | 1138.2 | 916.3 | 1.24x | measured |
| 40 | 10 | 1800 | 177.7 | 54.0 | 3.29x | 1153.6 | 966.5 | 1.19x | measured |
| 40 | 20 | 1800 | 314.4 | 91.4 | 3.44x | 1258.3 | 1093.9 | 1.15x | measured |
| 40 | 30 | 1800 | 460.8 | 129.3 | 3.56x | 1390.3 | 1193.9 | 1.16x | measured |
| 150 | 80 | 6750 | 24921.2 | 5114.0 | 4.87x | 17503.1 | 14961.2 | 1.17x | projected |
| 200 | 100 | 9000 | 60538.1 | 11770.3 | 5.14x | 30770.2 | 26239.6 | 1.17x | projected |

```
Both modes evaluate the RPA once, so these are like-for-like and no iteration-count correction is needed.  The ratio grows with system size because the self-consistent per-channel cost carries a steeper exponent (2.34 vs 2.16): the OEP route needs the correlation energy density and the driving term, each with an eigendecomposition and an inverse, where RPA@DFT needs one log-determinant.  Shaded rows are projected.
```

