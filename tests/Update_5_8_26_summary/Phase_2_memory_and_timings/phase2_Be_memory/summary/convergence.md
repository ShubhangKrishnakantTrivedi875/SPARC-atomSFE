# Energy and eigenvalue convergence vs N_fe (Z=4)

Lmax=30, p=20, q=45, domain=13 Bohr, omega=24, all-electron.  Each row is compared to the row directly above it -- no fitted or extrapolated reference.

## Self-consistent RPA

| N_fe | E (Ha) | dE vs previous (Ha) | HOMO (Ha) | dHOMO vs previous (Ha) |
|---:|---:|---:|---:|---:|
| 5 | -14.754500884 | - | -0.357223967 | - |
| 10 | -14.754506449 | -5.56e-06 | -0.355639836 | +1.58e-03 |
| 15 | -14.754506545 | -9.63e-08 | -0.355763421 | -1.24e-04 |
| 20 | -14.754506552 | -6.96e-09 | -0.355704309 | +5.91e-05 |
| 25 | -14.754506553 | -4.87e-10 | -0.355776240 | -7.19e-05 |
| 30 | -14.754506553 | +3.08e-10 | -0.355805204 | -2.90e-05 |

## Non-self-consistent RPA@DFT

| N_fe | E (Ha) | dE vs previous (Ha) | max |dE_eig| vs previous (Ha) |
|---:|---:|---:|---:|
| 5 | -14.752234543 | - | - |
| 10 | -14.752236013 | -1.47e-06 | +2.61e-06 |
| 15 | -14.752236098 | -8.52e-08 | +8.75e-10 |
| 20 | -14.752236105 | -6.60e-09 | +1.95e-09 |
| 25 | -14.752236104 | +5.79e-10 | +6.73e-08 |
| 30 | -14.752236104 | +1.68e-10 | +2.29e-08 |

`max |dE_eig| vs previous` is the worst-case eigenvalue change over ALL occupied states between consecutive N_fe, not just the HOMO.
