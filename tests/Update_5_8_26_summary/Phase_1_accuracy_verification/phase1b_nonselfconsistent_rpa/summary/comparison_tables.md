# Phase-1b comparison tables

Non-self-consistent RPA accuracy: `new` = corrected SPARC-atomSFE (xc='RPA@DFT', ground_state_functional='GGA_PBE'), `oneatom` = OneAtomFEM_rpa_at_dft emulation (xc='RPA', double_hybrid_flag=1, alpha_x=alpha_c=0, OEP_method skipped).  All-electron only.

**Total energy is not expected to agree** -- the two codes report structurally different quantities (see this file's module docstring).  Eigenvalues and RPA_correlation are the columns that actually validate the RPA implementation.


## HOMO eigenvalue, Ha

### loose

| element | new | oneatom | max |dE| (occ) |
|:---|---:|---:|---:|
| He (Z=2) | -0.579290748 | -0.579290748 | 1.71e-10 |
| N (Z=7) | -0.260724179 | -0.260724177 | 9.38e-09 |
| Mn (Z=25) | -0.184768742 | -0.184767802 | 5.44e-06 |
| Mo (Z=42) | -0.139281037 | -0.139285005 | 3.44e-05 |
| Au (Z=79) | -0.151377102 | -0.151344750 | 7.96e-05 |

### tight

| element | new | oneatom | max |dE| (occ) |
|:---|---:|---:|---:|
| He (Z=2) | -0.579290753 | -0.579290748 | 4.62e-09 |
| N (Z=7) | -0.260723970 | -0.260724193 | 2.67e-07 |
| Mn (Z=25) | -0.184767594 | -0.184767713 | 2.13e-06 |
| Mo (Z=42) | -0.139286760 | -0.139285613 | 5.25e-06 |
| Au (Z=79) | -0.151332685 | -0.151345664 | 2.97e-05 |


## E_c^RPA, Ha

### loose

| element | new | oneatom | |new-oneatom| |
|:---|---:|---:|---:|
| He (Z=2) | -0.083551555 | -0.083551555 | 3.28e-12 |
| N (Z=7) | -0.336934666 | -0.336934666 | 6.75e-11 |
| Mn (Z=25) | -1.525020797 | -1.525021436 | 6.39e-07 |
| Mo (Z=42) | -2.506810539 | -2.506813082 | 2.54e-06 |
| Au (Z=79) | -4.479239836 | -4.479226801 | 1.30e-05 |

### tight

| element | new | oneatom | |new-oneatom| |
|:---|---:|---:|---:|
| He (Z=2) | -0.084389678 | -0.084389678 | 1.52e-11 |
| N (Z=7) | -0.346794563 | -0.346794562 | 1.28e-09 |
| Mn (Z=25) | -1.652090140 | -1.652090219 | 7.87e-08 |
| Mo (Z=42) | -2.872240061 | -2.872240049 | 1.25e-08 |
| Au (Z=79) | -5.926641658 | -5.926636628 | 5.03e-06 |


## total energy, Ha (see caveat)

### loose

| element | new | oneatom | |new-oneatom| |
|:---|---:|---:|---:|
| He (Z=2) | -2.943693 | -2.976486 | 3.28e-02 |
| N (Z=7) | -54.171343 | -54.757932 | 5.87e-01 |
| Mn (Z=25) | -1150.171930 | -1152.076623 | 1.90e+00 |
| Mo (Z=42) | -3977.211527 | -3979.471959 | 2.26e+00 |
| Au (Z=79) | -17869.748206 | -17873.615316 | 3.87e+00 |

### tight

| element | new | oneatom | |new-oneatom| |
|:---|---:|---:|---:|
| He (Z=2) | -2.944531 | -2.977325 | 3.28e-02 |
| N (Z=7) | -54.181203 | -54.767792 | 5.87e-01 |
| Mn (Z=25) | -1150.298997 | -1152.203692 | 1.90e+00 |
| Mo (Z=42) | -3977.576944 | -3979.837386 | 2.26e+00 |
| Au (Z=79) | -17871.195607 | -17875.062726 | 3.87e+00 |

