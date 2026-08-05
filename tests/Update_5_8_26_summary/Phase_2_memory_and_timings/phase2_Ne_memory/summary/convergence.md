# Energy and eigenvalue convergence vs N_fe (Z=10)

Lmax=30, p=20, q=45, domain=13 Bohr, omega=24, all-electron.  Each row is compared to the row directly above it -- no fitted or extrapolated reference.

## Self-consistent RPA

| N_fe | E (Ha) | dE vs previous (Ha) | HOMO (Ha) | dHOMO vs previous (Ha) |
|---:|---:|---:|---:|---:|
| 5 | -129.146982455 | - | -0.812575848 | - |
| 10 | -129.147020587 | -3.81e-05 | -0.796110103 | +1.65e-02 |
| 15 | -129.147022056 | -1.47e-06 | -0.796147226 | -3.71e-05 |
| 20 | -129.147022128 | -7.16e-08 | -0.796108647 | +3.86e-05 |
| 25 | -129.147022137 | -9.05e-09 | -0.796156105 | -4.75e-05 |
| 30 | -129.147022138 | -1.49e-09 | -0.796176919 | -2.08e-05 |

## Non-self-consistent RPA@DFT

| N_fe | E (Ha) | dE vs previous (Ha) | max |dE_eig| vs previous (Ha) |
|---:|---:|---:|---:|
| 5 | -129.144730393 | - | - |
| 10 | -129.144758619 | -2.82e-05 | +2.84e-07 |
| 15 | -129.144759332 | -7.12e-07 | +7.01e-07 |
| 20 | -129.144759399 | -6.69e-08 | +6.50e-08 |
| 25 | -129.144759369 | +2.97e-08 | +2.65e-07 |
| 30 | -129.144759346 | +2.32e-08 | +1.84e-07 |

`max |dE_eig| vs previous` is the worst-case eigenvalue change over ALL occupied states between consecutive N_fe, not just the HOMO.
